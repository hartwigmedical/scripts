#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use File::Slurp;
use JSON::XS;
use File::Find::Rule;
use 5.010.000;

my $DATETIME = localtime;
my $SCRIPT = basename $0;
my $LIMS_JSON = '/data/ops/lims/prod/lims.json';
my $NONREPORTABLE_TSV = '/data/ops/lims/prod/samples_without_pdf_report.tsv';
my @PDF_BUCKETS = ( 'gs://patient-reporter-final-prod-1', 'gs://patient-reporter-final-prod-1/old_cpct_reports' );
my $PDF_BUCKETS_STRING = join( ", ", sort @PDF_BUCKETS );

my $DELIM = "\t";
my $NA_CHAR = "NA";
my $TYPE;
my $MUST_MATCH_EXACT;

my %OUT_FIELDS_PER_TYPE = (
  'samples'     => [qw(submission sample_id sample_name arrival_date label analysis_type project_name lab_sop_versions lab_status)],
  'submissions' => [qw(submission project_name project_type analysis_type sample_count has_lab_finished)],
  'contact_groups' => [qw(group_id report_contact_email data_contact_email lab_contact_email client_contact_email)],
  'unreported'  => [qw(sample_name lab_status lab_ref_status sample_id ref_sample_id sampling_date arrival_date api_info)],
  'dateTable'   => [qw(sample_name sampling_date arrival_date)],
);
my $available_types = join(", ", sort keys %OUT_FIELDS_PER_TYPE);
my $last_update_epoch = (stat($LIMS_JSON))[9]; 
my $last_update = strftime( "%Y-%m-%d %H:%M:%S", localtime($last_update_epoch));

my $HELP_TEXT =<<HELP;

  Description
    Parses LIMS JSON file and prints information about the
    requested samples / submissions to screen. Uses regex
    matching of filters by default.

  Updated: $last_update
    
  Usage
    $SCRIPT -type samples
    $SCRIPT -type submissions
    $SCRIPT -type contact_groups
    $SCRIPT -type samples -filter "submission=HMFreg"
    $SCRIPT -type samples -filter "sample_id=FR11111111"
    $SCRIPT -type samples -filter "sample_name=CPCT01010001T" -exact
    $SCRIPT -type unreported 
      (prints a list of biopsies without pdf report)

  Available query types:
     $available_types
    
  Selection options
    -exact      (match filters with "eq" instead of regex)
    -rna        (include rna samples in output)
    -plasma     (include plasma samples in output)
    
  Input/Output options:
    -delim  <s> (output delim)
    -lims   <s> (use alternate input json file)
    -json       (output objects in json format)
    
  Files/locations that are used
    LIMS input json file: $LIMS_JSON
    Patient report sources: $PDF_BUCKETS_STRING
HELP

die $HELP_TEXT . "\n" if scalar @ARGV == 0;

my %opt = ();
GetOptions (
    "type=s"    => \$TYPE,
    "exact"     => \$MUST_MATCH_EXACT,
    "delim=s"   => \$DELIM,
    "filter=s@" => \$opt{ filters },
    "rna"       => \$opt{ include_rna },
    "plasma"    => \$opt{ include_plasma },
    "lims=s"    => \$opt{ lims_input },
    "json"      => \$opt{ json },
    "debug"     => \$opt{ debug },
    "help|h"    => \$opt{ help },
) or die "Error in command line arguments\n";
die $HELP_TEXT . "\n" if $opt{'help'};
die $HELP_TEXT . "\n[ERROR] Please provide type with -type\n" unless $TYPE;
die $HELP_TEXT . "\n[ERROR] Type ($TYPE) is not supported\n" unless exists $OUT_FIELDS_PER_TYPE{ $TYPE };

my $nonreportable = readNonReportableSamples( $NONREPORTABLE_TSV );
my $lims = readJson( $opt{ lims_input } or $LIMS_JSON );
my $out_fields = $OUT_FIELDS_PER_TYPE{ $TYPE };

if ( $TYPE eq 'unreported' ){
    my $samples = getUnreportedBiopsies( $lims, $nonreportable, \%opt );
    printObjectInfo( $samples, $out_fields, 'lab_status', \%opt );
}
elsif ( $TYPE eq 'dateTable' ){
    my $samples = getDateTable( $lims, \%opt );
    printObjectInfo( $samples, $out_fields, 'sample_name', \%opt );
}
elsif( $TYPE eq 'samples' ){
    my $samples = filterSamples( $lims, \%opt  );
    printObjectInfo( $samples, $out_fields, 'sample_name', \%opt );
}
elsif( $TYPE eq 'submissions' ){
    my $submissions = filterByCategory( $lims, \%opt, $TYPE );
    printObjectInfo( $submissions, $out_fields, 'submission', \%opt );
}
elsif( $TYPE eq 'contact_groups' ){
    my $contact_groups = filterByCategory( $lims, \%opt, $TYPE );
    printObjectInfo( $contact_groups, $out_fields, 'group_id', \%opt );
}

sub filterByCategory{
    my ($lims, $opt, $category) = @_;
    my $objects = $lims->{ $category };
    my @selected_objects = ();
    foreach my $id ( keys %{$objects} ){
        push @selected_objects, $objects->{ $id };
    }
    my $filtered_objects = filterObjects( \@selected_objects, $opt );
    
    return $filtered_objects;
}

sub filterSamples{
    my ($lims, $opt) = @_;
    my $samples = $lims->{ 'samples' };
    
    my @selected_objects = ();    
    foreach my $sample_id ( keys %{$samples} ){
        my $sample = $samples->{ $sample_id };
        next if ( not $opt{ include_rna } ) and ($sample->{analysis_type} =~ /rna/i);
        next if ( not $opt{ include_plasma } ) and ($sample->{analysis_type} =~ /plasma/i);
        push @selected_objects, $sample;
    }
    my $filtered_objects = filterObjects( \@selected_objects, $opt );
    
    return $filtered_objects;
}

sub filterObjects{
    my ($objects, $opt) = @_;
    my $filters = $opt{ 'filters' };
    my $jsonOut = $opt{ 'json' };    

    my @out = ();
    my %filter_counts = ();
    
    foreach my $obj ( @$objects ){
        my $do_skip_object = applyFiltersOnObject( $obj, $filters, \%filter_counts );
        push( @out, $obj ) unless $do_skip_object;        
    }
    
    
    foreach my $filter ( keys %filter_counts ){
        my $count = $filter_counts{ $filter };
        unless ( $jsonOut ){
            say "## Filter: $count filtered away by filter \"$filter\"";
        }
    }
    
    return \@out;
}

sub applyFiltersOnObject{
    my ($object, $filters, $counts) = @_;
        
    foreach my $filter ( @$filters ){
        my $negate_search;
        my $field;
        my $match_str;

        if ($filter =~ /\!\=/){
            $negate_search = 1;
            ($field, $match_str) = split( "!=", $filter, 2 );
        }
        elsif ($filter =~ /=/){
            $negate_search = 0;
            ($field, $match_str) = split( "=", $filter, 2 );
        }
        else{
            die "[ERROR] Incorrect filter format ($filter). Requires a = or != in string.\n"
        }
                    
        if ( not exists $object->{ $field } ){
            $counts->{ $field.'=KeyNotExists' }++;
            return(1);
        }
        elsif ( not defined $object->{ $field } ){
            $counts->{ $field.'=KeyNotDefined' }++;
            return(1);
        }
        else{
            my $exact_match = $object->{ $field } eq $match_str;
            my $regex_match = $object->{ $field } =~ m/$match_str/i;

            if ( $negate_search ){
                $exact_match = !$exact_match;
                $regex_match = !$regex_match;
            }
            
            if ( $MUST_MATCH_EXACT ){
                if ( not $exact_match ){
                    $counts->{ $filter }++;
                    return(1);
                }
            }
            elsif( not $regex_match ){
                $counts->{ $filter }++;
                return(1);
            }
        }
    }
    
    ## all filters applied and still here so return OK response
    return 0;
}

sub readJson{
    my ($json_file) = @_;
    my $json_txt = read_file( $json_file );
    my $json_obj = decode_json( $json_txt );
    return( $json_obj );
}

sub readNonReportableSamples{
    my ($tsv_file) = @_;
    my %non_reportable_samples = ();
    open IN, "<", $tsv_file or die "[ERROR] Unable to open file ($tsv_file): $!\n";
    while ( my $sample = <IN> ){
        chomp $sample;
        next if $sample =~ /^#/;
        if ( $sample =~ /\s/ ){
            die "[ERROR] Whitespace detected in sample \"$sample\" (file: $tsv_file)\n";
        }
        $non_reportable_samples{ $sample } = 1;
    }
    close IN;
    return( \%non_reportable_samples );
}

sub printObjectInfo{
    my ($objects, $fields, $sortkey, $opt) = @_;
    my $jsonOut = $opt{ 'json' };
    my $count = scalar @$objects;
    
    unless ( $jsonOut ){
        say "## Script: $SCRIPT";
        say "## DateTime: $DATETIME";
        say "## InputJson: $LIMS_JSON";
    }
    
    @$objects = sort { 
        $a->{$sortkey} cmp $b->{$sortkey} or
        $a->{"sample_name"} cmp $b->{"sample_name"}
    } @$objects;
    
    if ( $jsonOut ){
        my $json_obj = JSON::XS->new->allow_nonref;
        my $json_txt = $json_obj->pretty->encode( $objects );
        print $json_txt;
    }
    else{
        say "## ObjectCount: $count";
        say "#".join( $DELIM, @$fields);
        foreach my $obj ( @$objects ){
            my @out = ();
            push @out, getValueByKey( $obj, $_) foreach @$fields;
            say join $DELIM, @out;
        }
    }
}

sub getDateTable{
    my ($lims) = @_;
    my @out_samples = ();
    my $samples = $lims->{ 'samples' };
    
    foreach my $sample_id ( keys %{$samples} ){

        my $sample = $samples->{ $sample_id };
        my $sample_name = getValueByKey( $sample, 'sample_name' );
        my $submission = getValueByKey( $sample, 'submission' );
        
        ## skip if not production biopsy
        next unless defined $sample->{'analysis_type'};
        next unless $sample->{'analysis_type'} =~ /^(Somatic_T|Somatic_R)$/;
        next unless defined $sample->{'submission'};
        next if $submission =~ /^HMFreg\d+$/ and $sample_name !~ /^CORE/;
        
        push @out_samples, $sample;
    }
    return( \@out_samples );
}


sub getUnreportedBiopsies{
    my ($lims, $nonreportable, $opt) = @_;
    my @out_samples = ();
    my $samples = $lims->{ 'samples' };
    
    my %pdfs = ();
    my @pdf_paths = ();

    foreach my $pdf_bucket ( @PDF_BUCKETS ){
        push( @pdf_paths, glob(`gsutil ls "$pdf_bucket/*.pdf"`) );
    }
    foreach my $pdf_path ( @pdf_paths ){
        my ($pdf_sample) = split( /[\.\_]/, basename( $pdf_path ) );
        $pdfs{ $pdf_sample } = 1;
    }
    
    foreach my $sample_id ( keys %{$samples} ){
        
        my $sample = $samples->{ $sample_id };
        my $sample_name = getValueByKey( $sample, 'sample_name' );
        my $submission = getValueByKey( $sample, 'submission' );
        
        ## skip if sample is on non reportable list
        next if exists $nonreportable->{ $sample_name };

        ## skip if report already exists
        next if exists $pdfs{ $sample_name };
        
        ## skip unless is biopsy
        next unless defined $sample->{'analysis_type'};
        next unless $sample->{'analysis_type'} eq 'Somatic_T';

        ## skip if sample is too old
        next unless defined $sample->{'arrival_date'};
        next if $sample->{'arrival_date'} =~ /(2015|2016|2017|2018|2019|2020)/;
                
        ## skip if not production biopsy
        next unless defined $sample->{'label'};
        next unless $sample->{'label'} =~ /^(CPCT|DRUP|WIDE|CORE|ACTN|SHRP|GAYA)$/;
        
        ## skip T0 biopsies (these should always have a T as well)
        next if $sample_name =~ /T0$/;
        
        ## skip if lab not ready yet
        next unless defined $sample->{'lab_status'};
        if ( $submission eq 'HMFregDRUP' ){
            ## DRUP will only be sequenced and reported when cohort is full
            next unless $sample->{'lab_status'} =~ /finished|failed/i;
        }
        else{
            next unless $sample->{'lab_status'} =~ /finished|failed|storage/i;
        }
        
        ## find and add R status
        $sample->{ 'ref_status' } = $NA_CHAR;
        if ( exists $sample->{ 'ref_sample_id' } ){
            my $ref_sample_id = $sample->{ 'ref_sample_id' };
            if ( exists $samples->{ $ref_sample_id } ){
                my $ref_sample = $samples->{ $ref_sample_id };
                my $ref_status = getValueByKey( $ref_sample, 'lab_status' );
                $sample->{ lab_ref_status } = $ref_status;
            }
        }
        
        ## ok: ready to report with API info added
        addApiInfoToSampleBySampleId($sample, $sample_id);
        push @out_samples, $sample;
    }
    
    return( \@out_samples );
}

sub addApiInfoToSampleBySampleId{
    my ($sample, $sample_id) = @_;

    my $api_cmd = "hmf_api_get 'runs?barcode=${sample_id}'";
    my $api_txt = `$api_cmd`;
    my $runs = decode_json( $api_txt );
    
    ## init new fields in case no runs found
    $sample->{api_info} = $NA_CHAR;

    foreach my $run ( @$runs ){
        my $ini = getValueByKey( $run, 'ini' );
        my $set_name = getValueByKey( $run, 'set', 'name' );
        my $bucket = getValueByKey( $run, 'bucket' );
        my $status = getValueByKey( $run, 'status' );
        next if $ini =~ /rerun/i;
        next if $bucket =~ /research-pipeline-output-prod/i;
        $sample->{api_info} = "$set_name ($status)";
    }
    return($sample);
}

sub getValueByKey{
    my ($obj, $key, $key2) = @_;
    my $out = $NA_CHAR;
    if ( defined $key2 and defined $obj->{ $key }{ $key2 }){
        $out = $obj->{ $key }{ $key2 };
    }elsif( defined $key and defined $obj->{ $key }){
        $out = $obj->{ $key };
    }
    return($out);
}
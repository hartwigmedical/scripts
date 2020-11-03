#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use File::Basename;
use Getopt::Long;
use File::Slurp;
use JSON::XS;
use Text::CSV;
use File::Copy;
use Email::Valid;
use 5.01000;

use constant EMPTY => q{ };
use constant NACHAR => 'NA';

## These fields are copied from ancestor sample
## in case they are not set at DNA/RNA sample
use constant FIELDS_TO_COPY_FROM_TISSUE_ANCESTOR => qw(
    report_germline report_germline_level report_viral report_pgx ptum
);

## These fields will be set to json boolean for json output
use constant BOOLEAN_FIELDS => qw(
    shallowseq report_germline report_viral report_pgx
    add_to_database add_to_datarequest
);

## Setup help msg
my $SCRIPT  = `basename $0`; chomp( $SCRIPT );
my $HELP_TEXT = <<"HELP";

  Description
    Parses LIMS text files (derived from excel and MS Access) and 
    writes to JSON output.
    
  Usage
    $SCRIPT -lims_dir /path/to/in/dir/ -out_json /path/to/out.json

  Required params:
    -lims_dir <str>  Path to dir with input files (eg /data/ops/lims/prod)
    -out_json <str>  Path to output json (eg /data/tmp/lims.json)
    
HELP

## Get input and setup all paths
my %opt = ();
GetOptions (
    "lims_dir=s" => \$opt{ lims_dir },
    "out_json=s" => \$opt{ out_json },
    "debug"      => \$opt{ debug },
    "help|h"     => \$opt{ help },
) or die("Error in command line arguments\n");

die $HELP_TEXT if $opt{ help };
die $HELP_TEXT unless $opt{ lims_dir };
die $HELP_TEXT unless $opt{ out_json };

my $CNTR_TSV = '/data/ops/lims/prod/center2entity.tsv';
my $LIMS_DIR = $opt{lims_dir};
my $JSON_OUT = $opt{out_json};
my $LATEST_DIR = $LIMS_DIR . "/lab_files/latest";

my $ACCESS_SAMPLES_CSV = $LATEST_DIR . '/access_samples.csv';
my $ACCESS_ACTIONS_CSV = $LATEST_DIR . '/access_actions.csv';
my $ACCESS_REGISTRATIONS_CSV = $LATEST_DIR . '/access_registration.csv';
my $FOR_001_SUBM_TSV = $LATEST_DIR . '/for001_submissions.tsv';
my $FOR_001_CONT_TSV = $LATEST_DIR . '/for001_contacts.tsv';
my $FOR_001_SAMP_TSV = $LATEST_DIR . '/for001_samples.tsv';
my $FOR_002_PROC_TSV = $LATEST_DIR . '/for002_processing.tsv';

## Files from previous years
my $SUBM_TSV_2019 = $LATEST_DIR . '/2019_for001_submissions.tsv';
my $SAMP_TSV_2019 = $LATEST_DIR . '/2019_for001_samples.tsv';
my $PROC_TSV_2019 = $LATEST_DIR . '/2019_for001_processing.tsv';
my $SUBM_TSV_2018 = $LATEST_DIR . '/2018_subm';
my $SAMP_TSV_2018 = $LATEST_DIR . '/2018_samp';
my $PROC_TSV_2018 = $LATEST_DIR . '/2018_proc';
my $PROC_TSV_2017 = $LATEST_DIR . '/2017_proc';
my $LIMS_JSN_2017 = $LATEST_DIR . '/2017_lims.json'; # excel LIMS pre-2018

my @ALL_INPUT_FILES = ( 
    $ACCESS_SAMPLES_CSV, $ACCESS_ACTIONS_CSV, $ACCESS_REGISTRATIONS_CSV,
    $FOR_001_SUBM_TSV, $FOR_001_SAMP_TSV, $FOR_002_PROC_TSV,
    $SUBM_TSV_2019, $SAMP_TSV_2019, $PROC_TSV_2019,
    $SUBM_TSV_2018, $SAMP_TSV_2018, $PROC_TSV_2018,
    $PROC_TSV_2017, $LIMS_JSN_2017,
    $CNTR_TSV
);

## Some final checks

foreach ( $LIMS_DIR ){
    die "[ERROR] Input dir does not exist ($_)\n" unless -e $_;
    die "[ERROR] Input dir is not a directory ($_)\n" unless -d $_;
}
foreach ( @ALL_INPUT_FILES ){
    die "[ERROR] Ihput file does not exist ($_)\n" unless -f $_;
}
foreach ( $JSON_OUT ){
    die "[ERROR] Output file exists and is not writable ($_)\n" if ( -f $_ and not -w $_ );
}
    
## ---------- 
## MAIN
## ----------
say "[INFO] START with \"$SCRIPT\"";

my $name_dict = getFieldNameTranslations();
my $cntr_dict = parseDictFile( $CNTR_TSV, 'center2centername' );

my $proc_objs = {}; # will contain objects from InProcess sheet
my $subm_objs = {}; # will contain objects from Received-Samples shipments sheet
my $cont_objs = {}; # will contain objects from Received-Samples contact sheet
my $samp_objs = {}; # will contain objects from Received-Samples samples sheet
my $cpct_objs = {}; # will contain objects from MS Access LIMS samples table
my $acti_objs = {}; # will contain objects from MS Access LIMS actions table
my $regi_objs = {}; # will contain objects from MS Access LIMS registration table
my $lims_objs = {}; # will contain all sample objects

$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2017, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2018, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2019, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $FOR_002_PROC_TSV, "\t" );

$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_2018'}, 'submission', 0, $SUBM_TSV_2018, "\t" );
$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_2019'}, 'submission', 0, $SUBM_TSV_2019, "\t" );

$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_CURR'}, 'submission', 0, $FOR_001_SUBM_TSV, "\t" );
$cont_objs = parseTsvCsv( $cont_objs, $name_dict->{'CONT_CURR'}, 'group_id',   1, $FOR_001_CONT_TSV, "\t" );

$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_2018'}, 'sample_id',  1, $SAMP_TSV_2018, "\t" );
$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_CURR'}, 'sample_id',  1, $SAMP_TSV_2019, "\t" );
$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_CURR'}, 'sample_id',  1, $FOR_001_SAMP_TSV, "\t" );

$cpct_objs = parseTsvCsv( $cpct_objs, $name_dict->{'CPCT_CURR'}, 'sample_id',  1, $ACCESS_SAMPLES_CSV, "," );
$acti_objs = parseTsvCsv( $acti_objs, $name_dict->{'ACTI_CURR'}, 'action_id',  1, $ACCESS_ACTIONS_CSV, "," );
$regi_objs = parseTsvCsv( $regi_objs, $name_dict->{'REGI_CURR'}, 'registration_id',  1, $ACCESS_REGISTRATIONS_CSV, "," );

checkContactInfo( $cont_objs );

$subm_objs = addContactInfoToSubmissions( $subm_objs, $cont_objs );
$lims_objs = addExcelSamplesToSamples( $lims_objs, $samp_objs, $subm_objs );
$lims_objs = addAccessSamplesToSamples( $lims_objs, $cpct_objs, $subm_objs, $cntr_dict );
$lims_objs = addLabSopStringToSamples( $lims_objs, $proc_objs );
$lims_objs = addIsoAndPrepExperimentIdsToSamples( $lims_objs, $regi_objs, $acti_objs );

checkDrupStage3Info( $subm_objs, $lims_objs );

printLimsToJson( $lims_objs, $subm_objs, $cont_objs, $JSON_OUT );

say "[INFO] DONE with \"$SCRIPT\"";

## ---------- 
## /MAIN
## ----------



## ---------- 
## SUBs
## ----------
sub parseTsvCsv{
    my ($objects, $fields, $store_field_name, $should_be_unique, $file, $sep) = @_;
    my $csv = Text::CSV->new({ binary => 1, auto_diag => 1, sep_char => $sep });
    my %store = %$objects;
    
    say "[INFO]   Parsing input file $file";
    open IN, "<", $file or die "Unable to open file ($file): $!\n";
    my $header_line = <IN>; chomp($header_line);
    die "[ERROR] Cannot parse line ($header_line)\n" unless $csv->parse($header_line);
    my @header_fields = $csv->fields();
    my %fields_map = map { $_ => 1 } @header_fields;
    
    ## Checking header content
    my $header_misses_field = 0;
    foreach my $field (keys %$fields) {
        if ( not exists $fields_map{ $field } ){
            warn "[WARN] Missing header field ($field) in file ($file)\n";
            $header_misses_field = 1;
        }
    }
    if ( $header_misses_field ){
        print Dumper \%fields_map and die "[ERROR] Header incomplete ($file)\n";
    }
    
    ## Header OK: continue reading in all data lines
    while ( <IN> ){
        chomp;
        die "[ERROR] Cannot parse line ($_)\n" unless $csv->parse($_);
        my @values = $csv->fields();
        
        my %raw_object = ();
        foreach my $field ( @header_fields ){
            my $next_value = shift @values;
            $next_value = NACHAR if not defined $next_value;
            $raw_object{ $field } = $next_value;
        }
        
        my $obj = selectAndRenameFields( \%raw_object, $fields );
        my $key = $obj->{ $store_field_name } || NACHAR;
        my $source = $obj->{ 'sample_source' } || NACHAR;
        my $name = $obj->{ 'sample_name' } || "SampleNameUnknown";
        
        if ( $source eq 'BLOOD' ){
            $key = $name;
        }
        
        next if isSkipValue( $key );
        my $reason_not_to_store = checkKeyToStore( \%store, $key );
        if ( $should_be_unique and $reason_not_to_store ){
            warn "[WARN] SKIPPING sample (name: $name) for reason: $reason_not_to_store\n" and next;
        }
        
        ## Checks OK: fix some fields and store object
        fixDateFields( $obj );
        fixIntegerFields( $obj );
        fixBooleanFields( $obj );
        $store{ $key } = $obj;
    }
    close IN;
    
    return \%store;
}

sub selectAndRenameFields{
    my ($obj_in, $fields) = @_;
    my %obj_out = ();
    foreach my $key ( keys %$obj_in ){
        if ( defined $fields->{ $key } ){
            $obj_out{ $fields->{$key} } = $obj_in->{ $key };
        }
    }
    return \%obj_out;
}

sub readJson{
    my ($json_file) = @_;
    say "[INFO]   Parsing input json file $json_file";
    my $json_txt = read_file( $json_file );
    my $json_obj = decode_json( $json_txt );
    return( $json_obj );
}

sub printLimsToJson{
    my ($samples, $submissions, $contact_groups, $lims_file) = @_;
    my $samp_count = scalar keys %$samples;
    my $subm_count = scalar keys %$submissions;
    my $cont_count = scalar keys %$contact_groups;

    my %lims = ( 'samples' => $samples, 'submissions' => $submissions, 'contact_groups' => $contact_groups );
    my $coder = JSON::XS->new->utf8->canonical;
    my $lims_txt = $coder->encode(\%lims);
    
    say "[INFO]   Writing output to $lims_file ($cont_count contact groups, $subm_count submissions and $samp_count samples)";
    open my $lims_json_fh, '>', $lims_file or die "Unable to open output file ($lims_file): $!\n";
        print $lims_json_fh $lims_txt;
    close $lims_json_fh;
}

sub addLabSopStringToSamples{
    my ($samples, $inprocess) = @_;
    my %store = %$samples;
    my $sop_field_name = 'lab_sop_versions';
    foreach my $id ( keys %store ){
        if ( exists $inprocess->{ $id } ){
            ## format: PREP(\d+)V(\d+)-QC(\d+)V(\d+)-SEQ(\d+)V(\d+)
            $store{ $id }{ $sop_field_name } = $inprocess->{ $id }{ $sop_field_name };
        }
        elsif ( defined $samples->{ $id }{ $sop_field_name } ){
            ## keep whatever is present
        }
        else{
            ## fallback to NA default
            $store{ $id }{ $sop_field_name } = NACHAR;
        }
    }
    return \%store;
}

sub addIsoAndPrepExperimentIdsToSamples{
    my ($samples, $registrations, $actions) = @_;
    my %store = %$samples;

    say "[INFO]   Adding isolation and prep dates to samples";

    my %conversion = (
        7  => { 'new_name' => 'isolation_start', lims_name => 'Compose blood isolation experiment'},
        8  => { 'new_name' => 'isolation_start', lims_name => 'Compose tissue isolation experiment'},
        11 => { 'new_name' => 'isolation_end', lims_name => 'Finished blood isolation experiment'},
        12 => { 'new_name' => 'isolation_end', lims_name => 'Finished tissue isolation experiment'},
    );

    my %experiment_dates = ();
    while (my ($registration_id, $obj) = each %$registrations){
        my $action_id = $obj->{"action_id"};
        my $date = $obj->{"date"};
        my $experiment = $obj->{"experiment_name"};
        next if not defined $experiment or isSkipValue($experiment);
        $experiment_dates{$experiment}{$action_id} = $date;
        say "--- DEBUG ---";
        say Dumper($experiment_dates{$experiment});
    }

    while (my ($sample_id, $sample_obj) = each %store) {
        while (my ($action_id, $action_name) = each %$actions) {
            next unless exists $conversion{$action_id};
            my $experiment = $sample_obj->{"isolation_exp_id"} || NACHAR;
            my $store_key = $conversion{$action_id}{'new_name'};
            $store{$sample_id}{$store_key} = NACHAR;

            if (defined $experiment_dates{$experiment}{$action_id}) {
                if (defined $store{$sample_id}{$store_key}){
                    next;
                }else{
                    $store{$sample_id}{$store_key} = $experiment_dates{$experiment}{$action_id};
                }
            }
        }
    }
    return \%store;
}

sub parseDictFile{
    my ($file, $fileType) = @_;
    say "[INFO]   Parsing input dictionary file $file";
    my %store = ();
    
    open my $dict_fh, "<", $file or die "$!: Unable to open file ($file)\n";
    while ( <$dict_fh> ){
        next if /^#/ms;
        chomp;
        if ( $fileType eq 'center2centername' ){
            my ( $id, $descr, $name ) = split /\t/;
            die "[ERROR] id occurs multiple times ($id) in file ($file)\n" if exists $store{ $id };
            $store{ $id } = $name if ( $id ne EMPTY and $name ne EMPTY );
        }
        else{
            die "[ERROR] File type not set or not recognized ($fileType)\n";
        }
    }
    close $dict_fh;
    
    return \%store;
}

sub addContactInfoToSubmissions{
    my ($submissions, $contact_groups) = @_;
    my @fields = qw( report_contact_name report_contact_email data_contact_name data_contact_email );
    say "[INFO]   Adding contact group information to submissions";
    my %store = %{$submissions};
    foreach my $submission_id (sort keys %store){
        my $submission = $store{$submission_id};
        
        if( exists $submission->{ 'report_contact_email' } ){
            ## skip records from the time when contact info was entered in shipments tab
            next;
        }
        elsif( defined $submission->{ 'group_id' } ){
            my $group_id = $submission->{ 'group_id' };
            next if $group_id eq 'na';
            if ( exists $contact_groups->{ $group_id } ){
                my $group = $contact_groups->{ $group_id };
                foreach my $field ( @fields ){
                    $submission->{$field} = $group->{ $field };
                }
            }
            else{
                warn "[WARN] Submission \"$submission_id\" has group ID \"$group_id\" but not found in contact groups (check Contact tab of FOR-001)\n";
                next;
            }
        }
        else{
            warn "[WARN] No Group ID field in place for submission $submission_id\n";
            next;
        }
    }
    return \%store;
}

sub checkDrupStage3Info{
    my ($submissions, $samples) = @_;
    my $drup_stage3_count = 0;
    my $expected_cohort = 'DRUPstage3';

    say "[INFO]   Checking $expected_cohort submissions";
    while (my ($id,$obj) = each %$submissions){
        my $project_type = $obj->{'project_type'};
        my $project_name = $obj->{'project_name'};

        next unless ($project_type eq 'Cohort' and $project_name =~ /DRUP/);

        $drup_stage3_count++;
        my $patient_id = $project_name;
        $patient_id =~ tr/-//d;
        my $sample_name = $patient_id . "T";
        my $sample_count = 0;
        while (my($barcode,$sample) = each %$samples){
            next unless $sample->{'sample_name'} eq $sample_name;
            next unless $sample->{'original_submission'} eq $id;
            next unless $sample->{'analysis_type'} eq 'Somatic_T';
            my $cohort = $sample->{'cohort'};
            if ( $cohort ne $expected_cohort ){
                say "[WARN]     Found sample for submission $id with incorrect cohort (name:$sample_name id:$barcode cohort:$cohort)";
            }
            $sample_count++;
        }
        say "[WARN]     Found no samples for $expected_cohort submission $id!" if $sample_count < 1;
    }
    say "[INFO]     Summary: $drup_stage3_count submissions encountered and checked for $expected_cohort";
}

sub checkContactInfo{
    my ($contact_groups) = @_;
    my @name_fields = qw(client_contact_name report_contact_name data_contact_name);
    my @mail_fields = qw(report_contact_email data_contact_email);
    say "[INFO]   Checking contact group information for completeness";
    foreach my $id (sort keys %$contact_groups){
        my $info = $contact_groups->{$id};
        my $name = $info->{client_contact_name};
        my $mail = $info->{client_contact_email};
        
        ## These fields should at the very least have content
        foreach my $field (@name_fields, @mail_fields){
            if ( $info->{$field} eq "" ){
                say "[WARN] No content in field \"$field\" for contact group ID \"$id\" (see FOR-001 Contacts tab)";
            }
        }
        ## These fields should contain only (valid) email addresses
        foreach my $field (@mail_fields){
            my @addressess = split( ";", $info->{$field});
            foreach my $address (@addressess){
                if( $address eq NACHAR ){
                    next;
                }
                elsif( not Email::Valid->address($address) ){
                    say "[WARN] No valid email address ($address) in field \"$field\"for contact group ID \"$id\" (see FOR-001 Contacts tab)";
                }
            }
        }
    }
}

sub addAccessSamplesToSamples{
    
    my ($lims, $objects, $submissions, $centers_dict) = @_;
    my %store = %{$lims};
    my %blood_samples_by_name = ();
    my %dna_blood_samples_by_name = ();
    my %tissue_samples_by_coupe = ();

    ## First gather info of certain samples to enrich DNA/RNA samples with later on
    say "[INFO]   Get blood and tissue sample info for later enrichment of DNA/RNA samples";
    foreach my $main_id ( sort keys %$objects ){
        my $row_info = $objects->{$main_id};    
        my $source = $row_info->{ 'sample_source' };        
        my $id     = $row_info->{ 'sample_id' }; # is field "Sample_barcode" in Access LIMS
        my $name   = $row_info->{ 'sample_name' };
        my $coupe  = $row_info->{ 'coupes_barcode' };
        my $status = $row_info->{ 'lab_status' };
        my $failed = $status =~ m/failed/i;

        ## Add to temporary mappings (used later on to enrich DNA/RNA samples)
        if ( $source eq 'BLOOD' and not isSkipValue($name) ){
            if ( exists $blood_samples_by_name{ $name } and $failed ){
                warn "[WARN]     Exclude mapping for BLOOD sample with name \"$name\" and id \"$id\" because status is failed";
            }else{
                $blood_samples_by_name{ $name } = $row_info;
            }
        }
        elsif ( $source eq 'TISSUE' and not isSkipValue($coupe) ){
            if ( exists $tissue_samples_by_coupe{ $coupe } and $failed ){
                warn "[WARN]     Exclude mapping for TISSUE sample with coupe barcode \"$coupe\" and id \"$id\" because status is failed";
            }else{
                $tissue_samples_by_coupe{ $coupe } = $row_info;
            }
        }
        elsif ( $source eq 'DNA-BLOOD' and not isSkipValue($name) ){
            if ( exists $dna_blood_samples_by_name{ $name } and $failed ){
                warn "[WARN]     Exclude mapping for DNA-BLOOD sample with name \"$name\" and id \"$id\" because status is failed";
            }else{
                $dna_blood_samples_by_name{ $name } = $row_info;
            }
        }
    }
    
    ## Now read through all samples again to store all samples 
    ## and enrich DNA and RNA samples with info from their ancestor samples
    while ( my($id, $object) = each %$objects ){
        
        my $name    = $object->{ 'sample_name' };
        my $patient = $object->{ 'patient' };
        my $source  = $object->{ 'sample_source' };
        my $coupe   = $object->{ 'coupes_barcode' } || NACHAR;
        
        (my $ancestor_coupe = $coupe) =~ s/(_DNA|_RNA)$//;
        my $is_dna = $source =~ /^DNA/ ? 1 : 0;
        my $is_rna = $source =~ /^RNA/ ? 1 : 0;
        my $is_blood = $source =~ /BLOOD$/ ? 1 : 0;
        my $is_tissue = $source =~ /TISSUE$/ ? 1 : 0;
        
        ## Only store DNA/RNA samples (all samples should have DNA derivative in LIMS)
        next unless ($is_dna or $is_rna);
        
        ## Source of a sample must be TISSUE or BLOOD or PLASMA (with possible DNA or RNA prefix for derivatives)
        unless ( $source =~ /^((DNA|RNA)\-)?(PLASMA|TISSUE|BLOOD)$/ ){
            warn "[WARN] SKIPPING sample: unable to parse source (\"$source\") for sample (\"$name\"). Need to fix this!\n" and next;
        }
        
        ## Skip particular cases of weird historic sample naming or absent name
        next if not $name or isSkipValue( $name );
        next if $name =~ /^(CPCT\d{8}|DRUP\d{8}|PMC\d{6})A/ms;
        next if $name =~ /^(PMC\d{6})(T|R){1}/ms;
        next unless $name =~ /^((CPCT|DRUP|WIDE|CORE)[0-9A-Z]{2}([0-9A-Z]{2})\d{4})(T|R){1}/ms;

        ## Enrich (DNA|RNA)-TISSUE samples with info from their TISSUE ancestor
        if ( $is_tissue ){
            if ( exists $tissue_samples_by_coupe{ $ancestor_coupe } ){
                my $ancestor = $tissue_samples_by_coupe{ $ancestor_coupe };
                foreach my $field ( FIELDS_TO_COPY_FROM_TISSUE_ANCESTOR ){
                    if ( exists $ancestor->{ $field } ){
                        $object->{ $field } = $ancestor->{ $field };
                    }else{
                        die "[ERROR] Trying to copy field (\"$field\") from ancestor sample but field does not exist.\n";
                    }
                }
            }
        
            ## If ref_sample_id absent: find blood counterpart by name
            my $ref_sample_id = $object->{ 'ref_sample_id' };
            if ( $ref_sample_id eq "" ){
                my $patient_string = $object->{ 'patient' };
                $patient_string =~ s/\-//g; # string this point still with dashes
                my $ref_sample_name = $patient_string . 'R';
                
                if ( exists $dna_blood_samples_by_name{ $ref_sample_name } ){
                    my $ref_sample = $dna_blood_samples_by_name{ $ref_sample_name };
                    $ref_sample_id = $ref_sample->{ 'sample_id' };
                    $object->{ 'ref_sample_id' } = $ref_sample_id;
                }
            }
        }
        elsif ( $is_blood ){
            # DNA-BLOOD samples are not enriched from BLOOD ancestor sample
        }
        else{
            warn "[WARN] Expected either tissue of blood sample but found neither (id:$id name:$name)\n";
        }

        ## All other sample names should follow certain formats
        if ( $name =~ /^((CPCT|DRUP|WIDE|CORE)[0-9A-Z]{2}([0-9A-Z]{2})\d{4})(T|R){1}/ms ){
            my ($patient_id, $study, $center, $tum_or_ref) = ($1, $2, $3, $4);
                
            $object->{ 'label' }   = $study;
            $object->{ 'patient' } = $patient_id;
            
            if ( $is_rna ){
                $object->{ 'analysis_type' } = "RNAanalysis";
            }elsif( $tum_or_ref eq 'T' ){
                $object->{ 'analysis_type' } = "Somatic_T";
            }elsif( $tum_or_ref eq 'R' ){
                $object->{ 'analysis_type' } = "Somatic_R";
            }else{
                $object->{ 'analysis_type' } = 'Unknown';
            }
            
            ## CORE is handled per case/submission
            if ( $study eq 'CORE' ){
                my $submission_id = $object->{ 'submission' };
                ## specifically check for non-ref samples if submission is defined
                if ( $submission_id eq '' ){
                    if ( $object->{ 'analysis_type' } ne 'Somatic_R' ){
                        warn "[WARN] SKIPPING CORE sample because of incorrect submission id \"$submission_id\" (id:$id name:$name)\n";
                    }
                    next;
                }
                $object->{ 'entity' } = $submission_id;
                $object->{ 'project_name' } = $submission_id;
            }
            ## All other samples are clinical study based (CPCT/DRUP/WIDE)
            elsif ( exists $centers_dict->{ $center } ){
                my $centername = $centers_dict->{ $center };
                my $original_submission = $object->{ 'submission' };
                my $register_submission = 'HMFreg' . $study;
                $object->{ 'original_submission' } = $original_submission;
                $object->{ 'submission' } = $register_submission;
                $object->{ 'project_name' } = $register_submission;
                $object->{ 'entity' } = join( "_", $study, $centername );
            }
            else {
                warn "[WARN] SKIPPING sample because not is core and center id unknown \"$center\" (id:$id name:$name)\n";
                next;
            }
        }
        else{
            warn "[WARN] SKIPPING sample from Access lims because of unknown sample name format (id:$id name:$name)\n" and next;
        }
 
        ## Final sanity checks before storing
        my $all_fields_ok = 1;
        ## Some fields need to be set for downstream tools not to crash!
        my @fields_to_check = qw( analysis_type submission );
        foreach my $field ( @fields_to_check ){
            if ( not exists $object->{ $field } or not defined $object->{ $field } or $object->{ $field } eq NACHAR ){
                warn "[WARN] SKIPPING sample because field $field is not present (id:$id name:$name)\n";
                $all_fields_ok = 0;
            }
            elsif ( $object->{ $field } eq NACHAR or $object->{ $field } eq "" ){
                warn "[WARN] SKIPPING sample because $field is not defined (id:$id name:$name)\n";
                $all_fields_ok = 0;
            }
        }
        next if not $all_fields_ok;
        
        ## Store sample unless key not OK
        my $reason_not_to_store = checkKeyToStore( \%store, $id );
        if ( $reason_not_to_store ){
            warn "[WARN] SKIPPING for reason: $reason_not_to_store" and next;
        }
        $store{ $id } = $object;
    }
    
    return \%store;
}

sub addExcelSamplesToSamples{
    
    my ($lims, $objects, $shipments) = @_;
    my %store = %{$lims};
    my %name2id = ();
    
    ## open file and check header before reading data lines
    while ( my($sample_id, $row_info) = each %$objects ){        
        
        my $sample_name = $row_info->{ 'sample_name' } or die "No sample_name in row_info";
        next if isSkipValue( $sample_name );

        my $sample_id = $row_info->{ 'sample_id' } or warn "[WARN] SKIPPING sample ($sample_name): No sample_id found" and next;
        my $submission = $row_info->{ 'submission' } or warn "[WARN] SKIPPING sample ($sample_name): No submission found" and next;
        my $analysis_type = $row_info->{ 'analysis_type' } or warn "[WARN] SKIPPING sample ($sample_name): No analysis_type found" and next;

        $row_info->{ 'label' } = 'RESEARCH';
        $row_info->{ 'patient' } = $row_info->{ 'sample_name' };
        $row_info->{ 'entity' } = $row_info->{ 'submission' };
        
        ## check data analysis type and set accordingly
        if ( $analysis_type =~ /^(Somatic_R|Somatic_T|SingleAnalysis|FASTQ|BCL|LabOnly)$/ ){
            ## Already final status so no further action
        }
        elsif ( $sample_name =~ /^(CORE\d{2}\d{6})(T|R){1}/ms ){
            my ($patient, $tum_or_ref) = ($1, $2);
            $row_info->{ 'label' }         = "CORE";
            $row_info->{ 'patient' }       = $patient;
            $row_info->{ 'entity' }        = $submission;
            $row_info->{ 'analysis_type' } = $tum_or_ref eq 'T' ? 'Somatic_T' : 'Somatic_R';
        }
        elsif ( $analysis_type eq 'SomaticAnalysis' or $analysis_type eq 'SomaticsBFX' ){
            ## SomaticsBFX is the old term, SomaticAnalysis the new
            my $partner = $row_info->{ 'ref_sample_id' };
            if ( $partner ne '' and $partner ne NACHAR ){
                $row_info->{ 'analysis_type' } = 'Somatic_T';
            }
            else{
                $row_info->{ 'analysis_type' } = 'Somatic_R';
            }
            ## TODO remove once shallowseq column is built into lims FOR-001
            ## Harcode Somatic_T samples to not run in shallow mode
            ## Harcode Somatic_T samples to not use existing ref data
            $row_info->{ 'shallowseq' } = 0;
            $row_info->{ 'other_ref' } = 0;
            
        }
        elsif ( $analysis_type eq 'GermlineBFX' or $analysis_type eq 'Germline' ){
            ## GermlineBFX is the old term, SingleAnalysis the new
            $analysis_type = 'SingleAnalysis';
            $row_info->{ 'analysis_type' } = $analysis_type;
        }
        elsif ( $analysis_type eq 'NoBFX' or $analysis_type eq 'NoAnalysis' or $analysis_type eq '' or $analysis_type eq 'NA' ){
            ## NoBFX is the old term, FASTQ the new
            $analysis_type = 'FASTQ';
            $row_info->{ 'analysis_type' } = $analysis_type;
        }
        elsif ( $analysis_type eq 'Labonly' ){
            $analysis_type = 'LabOnly';
            $row_info->{ 'analysis_type' } = $analysis_type;
        }
        elsif ( $analysis_type eq 'SNPgenotyping' or $analysis_type eq 'SNP' ){
            $analysis_type = 'SnpGenotyping';
            $row_info->{ 'analysis_type' } = $analysis_type;
        }
        else {
            warn "[WARN] SKIPPING sample ($sample_name): has unknown analysis type ($analysis_type)\n";
            next;
        }
        
        ## add submission info and parse KG
        if ( exists $shipments->{ $submission } ){
            my $sub = $shipments->{ $submission };
            my $project_name = $sub->{ 'project_name' };
            $row_info->{ 'project_name' } = $project_name;
            
            if ( $sub->{ 'project_type' } eq 'KG production' ){
                
                my @dvo_parts = split( /\-/, $project_name );
                my $center = uc( $dvo_parts[0] );
                $row_info->{ 'entity' } = 'KG_' . $center;
                $row_info->{ 'label' } = 'KG';
            }
            ## assume that all samples of submission need same analysis
            ## so will just overwrite analysis_type of submission
            $sub->{ 'analysis_type' } = $analysis_type;
        }

        my $unique = $row_info->{ 'sample_id' };
        next if isSkipValue( $unique );

        ## checks before storing
        my $regex = '^[0-9a-zA-Z\-]*$';
        warn "[WARN] SKIPPING sample ($sample_name): sample_name contains unacceptable chars\n" and next if $sample_name !~ /$regex/;
        warn "[WARN] SKIPPING sample ($sample_name): sample_id ($sample_id) contains unacceptable chars\n" and next if $sample_id !~ /$regex/;
        warn "[WARN] SKIPPING sample ($sample_name): no submission defined for sample\n" and next unless $row_info->{ 'submission' };
        warn "[WARN] SKIPPING sample ($sample_name): no analysis type defined for sample\n" and next unless $row_info->{ 'analysis_type' };
        warn "[WARN] SKIPPING sample ($sample_name): no project name defined for sample\n" and next unless $row_info->{ 'project_name' };        
        
        ## store at uniqe id
        my $reason_not_to_store = checkKeyToStore( \%store, $unique );
        if ( $reason_not_to_store ){
            warn "[WARN] SKIPPING sample with name \"$sample_name\" for reason: $reason_not_to_store" and next;
        }
        $store{ $unique } = $row_info;
        
    }
    
    return \%store;
}

sub fixIntegerFields{
    my ($obj) = @_;
    foreach my $key ( keys %$obj ){
        if ( $obj->{ $key } =~ /^\d+$/ ){
            $obj->{ $key } = $obj->{ $key } + 0;
        }
    }
}

sub fixBooleanFields{
    my ($obj) = @_;
    foreach my $key ( BOOLEAN_FIELDS ){
        next unless exists $obj->{ $key };
        next unless defined $obj->{ $key };
        my $value = $obj->{ $key };
        if ( $value =~ m/^true$/i ){
            $obj->{ $key } = JSON::XS::true;
        }elsif ( $value =~ m/^false$/i ){
            $obj->{ $key } = JSON::XS::false;
        }else{
            warn "[WARN] Unexpected value ($value) in boolean field ($key)\n"
        }
    }
}

sub fixDateFields{
    my ($obj) = @_;
    my @date_fields = qw( arrival_date sampling_date report_date );
    foreach my $date_field ( @date_fields ){
        
        next unless defined $obj->{ $date_field };
        my $old_date = $obj->{ $date_field };
        my $new_date = $old_date;
        my $identifier = $obj->{ 'sample_name' };
        
        ## date is not always filled in so skip NA fields
        next if isSkipValue( $old_date );
        
        ## Convert all date strings to same format yyyy-mm-dd (eg 2017-01-31)
        if( $old_date =~ /^\w+ (\w{3}) (\d{2}) \d+:\d+:\d+ \w+ (\d{4})$/ ){ 
            ## eg Tue Apr 23 00:00:00 CEST 2019
            my $month_name = $1;
            my $day = $2;
            my $year = $3;
            my $month = getMonthIndexByName( $month_name );
            $new_date = join( "-", $year, $month, $day );
        }
        elsif ( $old_date =~ /^(\d{2})(\d{2})(\d{2})$/ ){
            ## format unclear so need for checks
            warn "[WARN] Date \"$old_date\" in \"$date_field\" has unexpected year ($identifier): please check\n" if ($1 < 8) or ($1 > 20);
            warn "[WARN] Date \"$old_date\" in \"$date_field\" has impossible month ($identifier): please fix\n" if $2 > 12;
            $new_date = join( "-", "20" . $1, $2, $3 );
        }
        elsif ( $old_date =~ /^(\d{2})-(\d{2})-(\d{4})$/ ){
            ## dd-mm-yyyy
            warn "[WARN] Date \"$old_date\" in \"$date_field\" has impossible month ($identifier): please fix\n" if $2 > 12;
            $new_date = join( "-", $3, $2, $1 );
        }
        elsif ( $old_date =~ /^(\d{4})-(\d{2})-(\d{2})$/ ){
            ## case yyyy-mm-dd already ok
            warn "[WARN] Date \"$old_date\" in \"$date_field\" has impossible month ($identifier): please fix\n" if $2 > 12;
        }
        else{
            warn "[WARN] Date string \"$old_date\" in field \"$date_field\" has unknown format for sample ($identifier): kept string as-is but please fix\n";
        }
        
        ## store new format using reference to original location
        $obj->{ $date_field } = $new_date;
    }
}

sub getMonthIndexByName{
    my ($month_name) = @_;
    my %mapping = ( 
        "Jan" => "01", "Feb" => "02", "Mar" => "03", "Apr" => "04",
        "May" => "05", "Jun" => "06", "Jul" => "07", "Aug" => "08",
        "Sep" => "09", "Oct" => "10", "Nov" => "11", "Dec" => "12"
    );
    if ( exists $mapping{ $month_name } ){
        return $mapping{ $month_name };
    }
    else{
        warn "[WARN] Unknown Month name ($month_name): kept as-is but please fix\n";
        return $month_name;
    }
}

sub parseExcelSheet{
    my ($config) = @_;
    
    my $excel = $config->{ 'excel' };
    my $h_val = $config->{ 'h_val' };
    my $h_col = $config->{ 'h_col' };
    my $h_row = $config->{ 'h_row' };
    my $trans = $config->{ 'trans' };
    my $sheet = $config->{ 'sheet' };
    
    say "[INFO] Loading excel file $excel sheet \"$sheet\"";
    my $workbook = Spreadsheet::XLSX->new( $excel ) or die "[ERROR] Unable to load excel file $excel: $!\n";
    my $sheet_obj = $workbook->worksheet( $sheet ) or die "[ERROR] Unable to read sheet \"$sheet\" from file $excel: $!\n";
    
    my @header = ();
    my $max_row = $sheet_obj->{'MaxRow'};
    my $max_col = $sheet_obj->{'MaxCol'};
    
    ## check if header exist where it should be
    my $first_val = EMPTY;
    my $first_cel = $sheet_obj->get_cell( $h_row, $h_col );
    $first_val = $first_cel->unformatted() if defined $first_cel;
    die "[ERROR] Header value ($h_val) cannot be found at set location ($excel)\n" unless $first_val eq $h_val;
    
    ## now read header values for later storage
    foreach my $col ( $h_col .. $max_col ){
        my $cell = $sheet_obj->get_cell( $h_row, $col );
        my $cell_val = NACHAR;
        $cell_val = $cell->unformatted() if defined $cell;
        $cell_val = $trans->{ $cell_val } if defined $trans->{ $cell_val };
        push( @header, $cell_val );
    }
    
    return( \@header, $sheet_obj, $max_row, $max_col );
}

sub isSkipValue{
    my ($value) = @_;
    die "[ERROR] Value to check for skipping is not defined\n" if not defined $value;
    my @to_skip = ( NACHAR, EMPTY, '', 'na', 'naR', 'naT', 'invalid', 'failed', 'nvt', 'no', 'x' );
    foreach my $skip_string ( @to_skip ){
        return 1 if $value =~ /^$skip_string$/i;
    }
    return 0;
}

sub checkKeyToStore{
    my ($store, $key) = @_;    
    my $failReason = 0;
    
    if ( not defined $key ){
        $failReason = "key variable is not defined";
    }
    elsif ( isSkipValue($key) ){
        $failReason = "key is string to skip (key: $key)";
    }
    elsif ( $key =~ /[\n\r]/ ){
        my $woNewlines = $key =~ s/[\n\r\f]/\[ENTER\]/gr;
        $failReason = "key contains newline/enter (key: $woNewlines)";
    }
    elsif ( not $key =~ /^[a-zA-z0-9\-]+$/ ){
        $failReason = "key contains characters that are not allowed (key: $key)";
    }
    elsif ( exists $store->{ $key } ){
        $failReason = "duplicate key that should be unique (key: $key)";
    }
    
    return $failReason;
}

sub checkDefined{
	my ( $key, $hash) = @_;
	if ( not defined $hash->{$key} ){
		warn "[WARN] Value $key is not defined in:\n";
	    print Dumper $hash;
	}
}

sub getFieldNameTranslations{
    ## columns contact sheet in current FOR-001
    my %CONT_DICT = (
      "Group_ID"               => 'group_id',
      "Client_contact_name"    => 'client_contact_name',
      "Client_contact_email"   => 'client_contact_email',
      "On_behalf_of_client_name"  => 'on_behalf_of_client_contact_name',
      "On_behalf_of_client_email" => 'on_behalf_of_client_contact_email',
      "Report_contact_name"    => 'report_contact_name',
      "Report_contact_email"   => 'report_contact_email',
      "Data_contact_name"      => 'data_contact_name',
      "Data_contact_email"     => 'data_contact_email',
      "Lab_contact_name"       => 'lab_contact_name',
      "Lab_contact_email"      => 'lab_contact_email',
    );
    
    ## columns shipments sheet in 2018 rest lims (FOR-001)
    my %SUBM_DICT_2018 = (
      "Arrival_date"      => 'arrival_date',
      "Project_name"      => 'project_name',
      "HMF_reg"           => 'submission',
      "Requested_product" => 'request',
      "Product_category"  => 'project_type',
      "Sample_count"      => 'sample_count',
      "Lab_is_finished"   => 'has_lab_finished',
      "TAT_lab"           => 'turn_around_time',
      "Contact_name"      => 'report_contact_name',
      "Contact_email"     => 'report_contact_email',
      "Remarks"           => 'remarks',
      "Storage_status"    => 'lab_storage_status',
    );

    ## columns shipments sheet in 2019 rest lims (FOR-001)
    my %SUBM_DICT_2019 = (
      "Arrival_date"      => 'arrival_date',
      "Project_name"      => 'project_name',
      "HMF_reg"           => 'submission',
      "Requested_product" => 'request',
      "Product_category"  => 'project_type',
      "Sample_count"      => 'sample_count',
      "Lab_is_finished"   => 'has_lab_finished',
      "Group_ID"          => 'group_id',
      "TAT_lab"           => 'turn_around_time',
      "Contact_name"      => 'report_contact_name',
      "Contact_email"     => 'report_contact_email',
      "Portal_contact_name" => 'data_contact_name',
      "Portal_contact_email" => 'data_contact_email',
      "Remarks"           => 'remarks',
      "Storage_status"    => 'lab_storage_status',
    );
    
    ## columns shipments sheet in rest lims (FOR-001)
    my %SUBM_DICT = (
      "Arrival_date"      => 'arrival_date',
      "Project_name"      => 'project_name',
      "HMF_reg"           => 'submission',
      "Requested_product" => 'request',
      "Product_category"  => 'project_type',
      "Sample_count"      => 'sample_count',
      "Lab_is_finished"   => 'has_lab_finished',
      "Group_ID"          => 'group_id',
      "Remarks"           => 'remarks',
    );

    ## columns samples sheet in 2018 rest lims (FOR-001)
    my %SAMP_DICT_2018 = (
      "Sample_ID"         => 'sample_id',
      "Sample_name"       => 'sample_name',
      "DNA_conc"          => 'conc',
      "Yield"             => 'yield',
      "Q30"               => 'q30',
      "Analysis_type"     => 'analysis_type',
      "Partner_sample"    => 'ref_sample_id',
      "HMF_reg"           => 'submission',
      "SNP_required"      => 'is_snp_required',
      "SNP_exp"           => 'snp_experiment_id',
      "Requested_product" => 'request',
      "State"             => 'lab_status', # lab status
      "Primary_tumor_type"=> 'ptum',
      "Priority"          => 'priority', 
      "Arival_date"       => 'arrival_date',
      "Remarks"           => 'remarks',
    );

    ## columns samples sheet in CURRENT rest lims (FOR-001)
    my %SAMP_DICT = (
      "Sample_ID"         => 'sample_id',
      "Sample_name"       => 'sample_name',
      "DNA_conc"          => 'conc',
      "Yield"             => 'yield',
      "Q30"               => 'q30',
      "Analysis_type"     => 'analysis_type',
      "Partner_sample"    => 'ref_sample_id',
      "HMF_reg"           => 'submission',
      "SNP_required"      => 'is_snp_required',
      "SNP_exp"           => 'snp_experiment_id',
      "Requested_product" => 'request',
      "State"             => 'lab_status', # lab status
      "Priority"          => 'priority', 
      "Arrival_date"      => 'arrival_date',
      "Remarks"           => 'remarks',
    );

    ## columns In Process sheet (HMF-FOR-002)
    my %PROC_DICT = (
      'Sample_ID'         => 'sample_id', # eg FR12345678
      'Sample_name'       => 'sample_name', # eg CPCT1234567R
      'Diluted_library'   => 'library_id', # eg FR12345678 (THIS WAS "barcode_3nm")
      'Sop_tracking_code' => 'lab_sop_versions',
    );

    ## columns MS Access table samples
    my %CPCT_DICT = (
      'Coupes_barc'        => 'coupes_barcode',
      'Sampling_date'      => 'sampling_date',
      'Arrival_HMF'        => 'arrival_date',
      'Patient_name'       => 'patient',
      'Sample_name'        => 'sample_name',
      'Source'             => 'sample_source',
      'Yield'              => 'yield',
      'Sample_barcode'     => 'sample_id', # was Sample_ID_(DNA|RNA|Plasma)
      'Pathology_exp'      => 'pathology_exp',
      'Qiasymphony_exp'    => 'isolation_exp_id',
      'Prep'               => 'preparation_exp_id', # was (DNA|RNA)_prep
      'Purity_shallow'     => 'purity_shallow', # was Purity_shallow_(1|2|3)
      'Primary_tumor_type' => 'ptum',
      'tumor_'             => 'tumor_perc', # % in tumor_% is absent in export
      'Conc'               => 'conc', # was (DNA|RNA)_conc
      'Final_lab_status'   => 'lab_status',
      'Status_prep'        => 'prep_status',
      'Remarks'            => 'lab_remarks',
      'Other_Ref'          => 'other_ref',
      'Sample_ID_DNA_ref'  => 'ref_sample_id',
      'Hospital_patient_ID'=> 'hospital_patient_id',
      'Report_germline'    => 'report_germline',
      'Report_germline_level' => 'report_germline_level',
      'Submission_number'  => 'submission',
      'Date_of_birth'      => 'date_of_birth',
      'Purity'             => 'purity',
      'ShallowSeq'         => 'shallowseq', # Should we expect a ShallowSeq run
      'Database'           => 'add_to_database', # Should biopsy be added to SQL DB
      'Data_request'       => 'add_to_datarequest', # Should biopsy be part of DRs
      'Report_date'        => 'report_date',
      'Cohort'             => 'cohort', # What cohort is sample part of
      'Report_viral'       => 'report_viral', # Should we add viral info to report
      'Report_pgx'         => 'report_pgx', # Should we add pharmacogenomics info to report
      'Hospital_PA_sample_ID' => 'hospital_pa_sample_id',
      'Matching_other_HMF_patient_ID' => 'matching_other_HMF_patient_id',
      'QC_fail_report'     => 'qc_fail_report'
    );

    ## columns MS Access table actions
    my %ACTI_DICT = (
        'Action_ID' => 'action_id',
        'Action_description' => 'action_desc'
    );

    ## columns MS Access table registration
    my %REGI_DICT = (
        'ID' => 'registration_id',
        'Action_ID' => 'action_id',
        'Datum' => 'date',
        'User' => 'user',
        'Sample_name' => 'sample_name',
        'Experiment_name' => 'experiment_id',
    );

    my %translations = (
        'CONT_CURR' => \%CONT_DICT,
        'SUBM_CURR' => \%SUBM_DICT,
        'SUBM_2018' => \%SUBM_DICT_2018,
        'SUBM_2019' => \%SUBM_DICT_2019,
        'SAMP_CURR' => \%SAMP_DICT,
        'SAMP_2018' => \%SAMP_DICT_2018,
        'PROC_CURR' => \%PROC_DICT,
        'CPCT_CURR' => \%CPCT_DICT,
        'ACTI_CURR' => \%ACTI_DICT,
        'REGI_CURR' => \%REGI_DICT,
    );
    
    return \%translations;
}

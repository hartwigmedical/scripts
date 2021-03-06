#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use File::Basename;
use Getopt::Long;
use File::Slurp;
use Time::Piece;
use JSON;
use 5.010.000;

my $DATETIME = localtime;
my $SCRIPT = basename $0;
my $NA_CHAR = "NA";

my $GER_INI = "SingleSample.ini";
my $SOM_INI = "Somatic.ini";
my $SHA_INI = "ShallowSeq.ini";
my $RNA_INI = "Rna.ini";
my $BCL_INI = "BCL.ini";

my $Q30_LIM = 75; # q30 limit is currently fixed for all MS Access LIMS samples
my $YIELD_F = 1e9; # LAB lims contains yield in Gbase which needs to be converted to bases

## could make these more fine grained (lower number is higher prio)
my $NO_PIPELINE_PRIO = 100;
my $YES_PIPELINE_PRIO = 99;

my $LIMS_IN_FILE = '/data/ops/lims/prod/lims.json';
my $OUTPUT_DIR = '/data/ops/api/prod/jsons';
my $USE_EXISTING_REF = 0;
my $USE_EXISTING_TUM = 0;
my $FORCE_OUTPUT = 0;

## -----
## Gather input
## -----
my %opt = ();
GetOptions (
    "samplesheet=s"  => \$opt{ samplesheet },
    "outdir=s"       => \$OUTPUT_DIR,
    "limsjson=s"     => \$LIMS_IN_FILE,
    "useExistingRef" => \$USE_EXISTING_REF,
    "useExistingTum" => \$USE_EXISTING_TUM,
    "forceOutput"    => \$FORCE_OUTPUT,
    "experiment"     => \$opt{ as_experiment },
    "debug"          => \$opt{ debug },
    "help|h"         => \$opt{ help },
    "verbose"        => \$opt{ verbose }
) or die "Error in command line arguments\n";
my @ids = @ARGV;

my $HELP =<<HELP;

  Description
    Creates the JSON(s) to register set/run at HMF API.
    Output directory = $OUTPUT_DIR
    
  Usage
    $SCRIPT -samplesheet \${path-to-samplesheet}
      eg: $SCRIPT -samplesheet /data/run/SampleSheet.csv
    
    $SCRIPT \${sample_id} [\${sample_id2} \${sample_id_n}]
      eg: $SCRIPT FR12345678
      eg: $SCRIPT FR11111111 FR22222222 FR33333333
    
  Options
    -outdir <s>        [$OUTPUT_DIR]
    -limsjson <s>      [$LIMS_IN_FILE]
    -useExistingRef    (add use_existing_sample flag for ref sample in json)
    -useExistingTum    (add use_existing_sample flag for tum sample in json)
    -forceOutput       (write json even if sample already has a run in API)
    -experiment <s>    (resets submission to HMFregVAL and entity to HMF_EXPERIMENT)

HELP

print $HELP and exit(0) if $opt{ help };
print $HELP and exit(0) if scalar(@ids) == 0 and not defined $opt{ samplesheet };
die "[ERROR] Output dir is not writeable ($OUTPUT_DIR)?\n" unless -w $OUTPUT_DIR;

## MAIN
sayInfo("START of script $SCRIPT ($DATETIME)");

if ( defined $opt{ samplesheet } ){
    sayInfo("Reading SampleSheet file ($opt{ samplesheet })");
    my $ids_from_sheet = addSamplesFromSamplesheet( $opt{ samplesheet } );
    push( @ids, @$ids_from_sheet );
}

sayInfo("InputCount: ".scalar(@ids));
sayInfo("Reading LIMS file ($LIMS_IN_FILE)");
my $lims = readJson( $LIMS_IN_FILE );
my $samples = $lims->{ 'samples' };
my %stats = ();

foreach my $sample_id ( @ids ){
    sayInfo("Processing $sample_id");
    my $return = processSample( $sample_id, $samples);
    $stats{ $return }++;
}

sayInfo("STATS -----")  ;

foreach my $reason ( keys %stats ){
    my $count = $stats{ $reason };
    sayInfo("Stats:   $reason = $count");
}

## -----
## /MAIN
## -----
    
sub addSamplesFromSamplesheet{
    my ($file) = @_;
    
    my %head = ();
    my %data = ();
    my $currblock = '';
    
    ## first read file to obtain header fields
    my @header;
    open my $header_fh, '<', $file or die "Unable to open file ($file): $!\n";
    while ( <$header_fh> ){
        chomp;
        if ( $_ =~ /^\[Data\]/ ){
            my $header_line = <$header_fh>;
            die "[ERROR] Header line should contain Sample_ID\n" unless $header_line =~ /Sample_ID/;
            @header = split( ",", $header_line );
        }
    }
    close $header_fh;
    die "[ERROR] No header line was parsed?\n" unless scalar @header;
    
    ## re-read file to get all information
    open my $samplesheet_fh, '<', $file or die "Unable to open file ($file): $!\n";
    while ( <$samplesheet_fh> ){
        chomp;
        next if $_ eq '' or $_ =~ /^[\,\s]+$/;
        if ( $_ =~ /^\[(Header|Reads|Settings|Data)\]/ ){
            $currblock = $1;
        }
        elsif ( $currblock eq 'Header' ){
            my ($field, $content) = split /\,/;
            $head{ $field } = $content;
        }
        elsif ( $currblock eq 'Data' ){
            next if $_ =~ /Sample_ID/; # skip data header line
            my @line_values = split( ',', $_ );
            my %record = ();
            foreach my $field ( @header ){
                $record{ $field } = shift @line_values;
            }
            my $id = $record{ 'Sample_ID' };
            my $name = $record{ 'Sample_Name' };
            my $submission = $record{ 'Sample_Project' };

            ## VAL and GIAB samples are not present in LIMS so need manual work
            if ($submission eq "HMFregVAL"){
                sayWarn("SKIPPING sample ($name, $id) because of unsupported submission in SampleSheet ($submission)");
            }
            elsif ($submission eq "HMFregGIAB"){
                sayWarn("SKIPPING sample ($name, $id) because of unsupported submission in SampleSheet ($submission)");
            }
            else{
                $data{ $id } = 1;
            }
        }
    }
    close $samplesheet_fh;
    
    my $hmfRunName = $head{ 'ExperimentName' } || $NA_CHAR;
    sayInfo("Found run $hmfRunName in SampleSheet");
    my @out = sort keys %data;
    return( \@out );
}

sub processSample{
    my ($sample_id, $lims_samples) = @_;
    my @warn_msg = ();
    if ( not exists $lims_samples->{ $sample_id } ){
        sayWarn("  RESULT: Sample not present in LIMS ($sample_id)");
        return "NoJsonMade_sampleDoesNotExistsInLims";
    }
    my $sample = $lims_samples->{ $sample_id };
    
    my $name       = getValueByKey( $sample, 'sample_name' ); # eg CPCT02010000R
    my $barcode    = getValueByKey( $sample, 'sample_id' ); # eg FR12345678
    my $patient    = getValueByKey( $sample, 'patient' ); # eg CPCT02010000
    my $submission = getValueByKey( $sample, 'submission' ); # eg HMFregCPCT
    my $analysis   = getValueByKey( $sample, 'analysis_type' ); # eg Somatic_T
    my $entity     = getValueByKey( $sample, 'entity' ); # eg HMFreg0001
    my $priority   = getPriorityForSample( $sample );
    my $yield      = getValueByKey( $sample, 'yield' ) * $YIELD_F;
    
    ## reset 0 yield to 1 base in order to avoid samples being ready directly
    ## except for so-called "VirtualSample" samples (these index seqs should be absent)
    if ( $yield == 0 and $name !~ /^VirtualSample\d+/ ){
        $yield = 1;
    }

    ## overwrite submission and entity in case of experiment
    ## this also makes sure output goes to experiment bucket via entity
    if ( $opt{'as_experiment'} ){
        $submission = 'HMFregVAL';
        $entity = 'HMF_EXPERIMENT';
    }
    
    my $use_existing_ref = $USE_EXISTING_REF;
    my $use_existing_tum = $USE_EXISTING_TUM;

    ## not all samples have q30 field because this was added later to lims
    my $q30 = $Q30_LIM;
    if ( defined $sample->{ 'q30' } ){
        $q30 = $sample->{ 'q30' };
    }
    if ( $q30 !~ /^\d+$/ or $q30 < 0 or $q30 > 100 ){
        die "[ERROR] Q30 found for sample ($name) but not an integer percentage ($q30)\n";
    }
    
    ## init the json info
    my %json_data = ();

    sayInfo("  NAME=$name, ENTITY=$entity, ANALYSIS:$analysis");

    my $date = localtime->strftime('%y%m%d');

    ## Setup json content based on analysis type
    if ( $analysis eq 'BCL' ){
        my $set = join("_", $date, $submission, $barcode, $name );
        sayInfo("  SET: $set");
        $json_data{ 'ini' } = "$BCL_INI";
        $json_data{ 'set_name' } = "$set";
        $json_data{ 'entity' } = "$entity";
        addSampleToJsonData( \%json_data, $submission, $barcode, $name, 'ref', $q30, $yield, $use_existing_ref );
    }
    elsif ( $analysis eq 'FASTQ' ){
        my $set = join("_", $date, $submission, $barcode, $name );
        sayInfo("  SET: $set");
        $json_data{ 'set_name' } = "$set";
        $json_data{ 'entity' } = "$entity";

        $json_data{ 'fastq_portal' } = JSON::true;
        addSampleToJsonData( \%json_data, $submission, $barcode, $name, 'ref', $q30, $yield, $use_existing_ref );
    }
    elsif( $analysis eq 'SingleAnalysis' ){
        my $set = join( "_", $date, $submission, $barcode, $name );
        sayInfo("  SET: $set");
        $json_data{ 'ini' } = "$GER_INI";
        $json_data{ 'set_name' } = "$set";
        $json_data{ 'entity' } = "$entity";
        addSampleToJsonData( \%json_data, $submission, $barcode, $name, 'ref', $q30, $yield, $use_existing_ref );
    }
    elsif ( $analysis eq 'RNAanalysis' ){
        my $set = join( "_", $date, "HMFregRNA", $barcode, $name );
        # TODO: remove default 15 GBase yield once yield is correctly set in LIMS for RNA
        my $default_rna_yield = 15e9;
        sayInfo("  SET: $set");
        $yield = $default_rna_yield;
        $json_data{ 'ini' } = "$RNA_INI";
        $json_data{ 'set_name' } = "$set";
        $json_data{ 'entity' } = "$entity";
        addSampleToJsonData( \%json_data, $submission, $barcode, $name, 'tumor-rna', $q30, $yield, $use_existing_ref );
    }
    elsif ( $analysis eq 'Somatic_T' ){
        
        my $ref_obj;
        my $ini = $SOM_INI;
        my $needs_shallow = getValueByKey( $sample, 'shallowseq' ); # 0 | 1
        my $other_ref = getValueByKey( $sample, 'other_ref' ); # Yes | No
        
        ## need to find the ref sample of somatic pair
        if ( exists $sample->{ ref_sample_id } and $sample->{ ref_sample_id } ne "" ){
            ## for somatic samples (biopsy) a ref sample needs to be defined
            my $ref_sample_id = $sample->{ ref_sample_id };
            $ref_obj = getSomaticRSampleByStringForField( $lims_samples, $ref_sample_id, 'sample_id' );
        }
        else{
            ## for clinical studies the partner sample can be found by patientID + R
            my $ref_string = $patient.'R';
            $ref_obj = getSomaticRSampleByStringForField( $lims_samples, $ref_string, 'sample_name' );
        }

        if ( not defined $ref_obj ){
            sayWarn("  RESULT: SKIPPING because somatic R not found for input T (PATIENT=$patient)");
            return "NoJsonMade_RnotFoundForSomaticT";
        }
        
        my $barcode_ref = getValueByKey( $ref_obj, 'sample_id' );
        my $name_ref = getValueByKey( $ref_obj, 'sample_name' );
        my $yield_ref = getValueByKey( $ref_obj, 'yield' );
        my $submission_ref = getValueByKey( $ref_obj, 'submission' );
        $yield_ref = $yield_ref == 0 ? 1 : $yield_ref * $YIELD_F;
        my $set = join( "_", $date, $submission, $barcode_ref, $barcode, $patient );

        ## adjust content in case of ShallowSeq
        if ( $needs_shallow ){
            sayInfo("  ShallowSeq flag set in LIMS");
            my $match_string = join( "_", "ShallowSeq", $barcode_ref, $barcode, $name );

            my $api_shallow_runs = decode_json(`hmf_api_get 'runs?ini=ShallowSeq.ini&barcode=$barcode'`);
            my $run_status = '';
            if ( scalar @$api_shallow_runs > 0 ) {
                my $most_recent_api_run = $api_shallow_runs->[-1];
                my $run_name = $most_recent_api_run->{set}->{name};
                if ($run_name =~ /$match_string/) {
                    $run_status = $most_recent_api_run->{status};
                }
            }

            if ( $run_status eq "" ){
                say "[INFO]   No ShallowSeq run found in api for $match_string: going for ShallowSeq mode";
                $yield = 35 * $YIELD_F;
                $yield_ref = $yield;
                $entity = 'HMF_SHALLOWSEQ';
                $set = join( "_", $date, "ShallowSeq", $barcode_ref, $barcode, $name );
                $ini = $SHA_INI;
                $priority = $YES_PIPELINE_PRIO;
            }
            elsif ( $run_status eq "Waiting" ){
                sayWarn("  RESULT: SKIPPING because run found for $match_string with status $run_status (so assuming extra seq)");
                return "NoJsonMade_ShallowExtraSeq";
            }
            elsif ( $run_status eq "Processing" ){
                sayWarn("  RESULT: SKIPPING because run found for $match_string with status $run_status (so assuming accidental extra seq)");
                return "NoJsonMade_ShallowSeqStillProcessing";
            }
            elsif ( $run_status eq "Finished" ){
                sayInfo("  ShallowSeq run with status $run_status found for $match_string: going for full Somatic mode");
            }
            else{
                die "[ERROR]   ShallowSeq with status $run_status found for $match_string: no idea what to do";
            }
        }

        ## add suffix to ref barcode and use tumor submission in case ref is needed from other existing patientId
        if ( $other_ref eq "Yes" ){
            my $new_name_ref = $patient . 'R';
            my $new_barcode_ref = $barcode_ref . "_" . $new_name_ref;
            push( @warn_msg, "DOUBLE CHECK JSON for $barcode ($name): OtherRef flag is set in LIMS so adding suffix to the REF barcode ($new_barcode_ref)" );
            $name_ref = $new_name_ref;
            $barcode_ref = $new_barcode_ref;
            $submission_ref = $submission;
        }
 
        ## check if barcode already exists in HMF API
        my $tum_exists = `hmf_api_get 'samples?barcode=$barcode' | jq 'select(.[].status != "Unregistered") | length'`;
        my $ref_exists = `hmf_api_get 'samples?barcode=$barcode_ref' | jq 'select(.[].status != "Unregistered") | length'`;
        chomp($tum_exists);
        chomp($ref_exists);
        
        if ( $tum_exists ){
            sayInfo("  TUM barcode ($barcode) already exists in HMF API (so will add use_existing flag)");
            $use_existing_tum = 1;
        }        
        if ( $ref_exists ){
            sayInfo("  REF barcode ($barcode_ref) already exists in HMF API (so will add use_existing flag)");
            $use_existing_ref = 1;
        }

        sayInfo("  Set name constructed to $set");
        $json_data{ 'ini' } = "$ini";
        $json_data{ 'set_name' } = "$set";
        $json_data{ 'entity' } = "$entity";
        $json_data{ 'priority' } = $priority;
        addSampleToJsonData( \%json_data, $submission_ref, $barcode_ref, $name_ref, 'ref', $q30, $yield_ref, $use_existing_ref );
        addSampleToJsonData( \%json_data, $submission, $barcode, $name, 'tumor', $q30, $yield, $use_existing_tum );
    }
    elsif ( $analysis eq 'Somatic_R' ){
        sayInfo("  RESULT for $sample_id: SKIPPING because is somatic ref sample ($name)");
        return "NoJsonMade_isSomaticR";
    }
    else{
        sayWarn("  RESULT for $sample_id: Somehow no (correct) analysis type ($analysis) was defined for input");
        return "NoJsonMade_hasWrongAnalsisType";
    }

    ## output json
    my $json_file = $json_data{ 'set_name' }.'.json';
    my $json_path = $OUTPUT_DIR.'/'.$json_file;
    
    ## check if set was already registered earlier
    my $setname_wo_date = $json_data{ 'set_name' };
    $setname_wo_date =~ s/^\d{6}_//;
    my $existing_runs = decode_json(`hmf_api_get 'runs?set_name_contains=$setname_wo_date'`);
    my $existing_count = scalar @$existing_runs;

    if ( $existing_count > 0 ){
        if ( $FORCE_OUTPUT ){
            push( @warn_msg, "Existing run(s) found in API for $setname_wo_date ($existing_count) but output was enforced" );
        }
        else{
            sayWarn("  RESULT for $sample_id ($name): SKIPPING because a run with '$setname_wo_date' already exists in API");
            return "NoJsonMade_setJsonAlreadyExists";
        }
    }
    
    ## print any stored warnings now
    if ( scalar @warn_msg ){
        foreach my $msg ( @warn_msg ){
            sayWarn("  $msg");
        }
    }
    
    ## all checks were OK: print config file
    printSetJson( \%json_data, $json_path );
    sayInfo("  RESULT for $sample_id: OK");
    return "OK_JSON_MADE";
}

sub printSetJson{
    my ($data, $out_path) = @_;
    my $json_obj = JSON->new->allow_nonref;
    my $json_txt = $json_obj->pretty->encode( $data );
    sayInfo("  Writing json ($out_path)");
    open OUT, '>', $out_path or die "Unable to open output file ($out_path): $!\n";
        print OUT $json_txt;
    close OUT;
}

sub readJson{
    my ($json_file) = @_;
    my $json_txt = read_file( $json_file );
    my $json_obj = decode_json( $json_txt );
    return( $json_obj );
}

sub getSomaticRSampleByStringForField{
    my ($info, $search_string, $search_field) = @_;
    
    foreach my $sample_id ( keys %$info ){
        my $field_value = $info->{ $sample_id }{ $search_field };
        if (( $field_value eq $search_string ) and ( $info->{ $sample_id }{ 'analysis_type' } eq 'Somatic_R')){
            return $info->{ $sample_id };
        }
    }
    sayWarn("$search_string not found in field $search_field of any record");
    return(undef);
}

sub getValueByKey{
    my ($info, $key) = @_;
    if ( not defined $info->{ $key } ){
        say "[ERROR] Cannot find field \"$key\" in data structure:";
        print Dumper( $info );
        die "[ERROR] Unable to get field $key\n"
    }
    else{
        return( $info->{ $key } );
    }
}

sub getPriorityForSample{
    my ($sample_info) = @_;
    ## unfortunately cannot err on key absence 
    ## because not all samples have the prio property
    if ( defined $sample_info->{ 'priority' } and $sample_info->{ 'priority' } =~ /yes/i ){
        return $YES_PIPELINE_PRIO;
    }
    else{
        return $NO_PIPELINE_PRIO;
    }
}

sub addSampleToJsonData{
    my ($store, $submission, $barcode, $name, $type, $q30, $yield, $use_existing) = @_;
    my %tmp = (
        'barcode'    => "$barcode",
        'name'       => "$name",
        'submission' => "$submission",
        'type'       => "$type",
        'q30_req'    => int($q30),
        'yld_req'    => int($yield),
    );
    if ( $use_existing ){
        $tmp{ 'use_existing_sample' } = JSON::true;
    }
    push( @{$store->{ 'samples' }}, \%tmp );
}

sub sayInfo{
    my ($msg) = @_;
    say "[INFO] $msg"
}

sub sayWarn{
    my ($msg) = @_;
    warn "[WARN] $msg\n"
}


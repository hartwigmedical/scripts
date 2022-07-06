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
my $TAR_INI = "Targeted.ini";

my $Q30_LIM = 75; # q30 limit is currently fixed for all MS Access LIMS samples
my $YIELD_F = 1e9; # LAB lims contains yield in Gbase which needs to be converted to bases

my $NO_PIPELINE_PRIO = 100;
my $YES_PIPELINE_PRIO = 99;

my $LIMS_IN_FILE = '/data/ops/lims/prod/lims.json';
my $OUTPUT_DIR = '/data/ops/api/prod/jsons';
my $USE_EXISTING_REF = 0;
my $USE_EXISTING_TUM = 0;
my $SKIP_RECALCULATING_YIELD_REF = 0;
my $SKIP_RECALCULATING_YIELD_TUM = 0;
my $IGNORE_SHALLOWSEQ = 0;
my $FORCE_OUTPUT = 0;

my %opt = ();
GetOptions (
    "samplesheet=s"  => \$opt{ samplesheet },
    "outdir=s"       => \$OUTPUT_DIR,
    "limsjson=s"     => \$LIMS_IN_FILE,
    "useExistingRef" => \$USE_EXISTING_REF,
    "useExistingTum" => \$USE_EXISTING_TUM,
    "noShallow"      => \$IGNORE_SHALLOWSEQ,
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
    -noShallow         (do not create ShallowSeq run even if set in LIMS)
    -experiment        (resets submission to HMFregVAL and entity to HMF_EXPERIMENT)

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
            my $header_line = removeNewlineCharacters(<$header_fh>);
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
        my $line = removeNewlineCharacters($_);
        next if $line eq '' or $line =~ /^[\,\s]+$/;
        if ( $line =~ /^\[(Header|Reads|Settings|Data)\]/ ){
            $currblock = $1;
        }
        elsif ( $currblock eq 'Header' ){
            my ($field, $content) = split /\,/;
            $head{ $field } = $content;
        }
        elsif ( $currblock eq 'Data' ){
            next if $line =~ /Sample_ID/; # skip data header line
            my @line_values = split( ',', $line );
            my %record = ();
            foreach my $field ( @header ){
                $record{ $field } = shift @line_values;
            }
            my $id = $record{ 'Sample_ID' };
            my $name = $record{ 'Sample_Name' };
            my $submission = $record{ 'Sample_Project' };

            if (! defined $submission || $submission eq ""){
                die "[ERROR] No submission found in line: $line\n";
            }

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
    my $lab_status = getValueByKey( $sample, 'lab_status' );
    
    ## reset 0 yield to 1 base in order to avoid samples being ready directly
    ## except for so-called "VirtualSample" samples (these index seqs should be absent)
    if ( $yield == 0 and $name !~ /^VirtualSample\d+/ ){
        $yield = 1;
    }

    if ( $lab_status ne "Finished" and $lab_status ne "finished" and $lab_status ne "In process" and $name !~ /^VirtualSample\d+/){
        sayWarn("  Check JSON for sample $name, since lab status is unexpected: lab_status='$lab_status'")
    }

    ## overwrite submission and entity in case of experiment
    ## this also makes sure output goes to experiment bucket via entity
    if ( $opt{'as_experiment'} ){
        $submission = 'HMFregVAL';
        $entity = 'HMF_EXPERIMENT';
    }
    
    my $use_existing_ref = $USE_EXISTING_REF;
    my $use_existing_tum = $USE_EXISTING_TUM;
    my $skip_recalculating_yield_ref = $SKIP_RECALCULATING_YIELD_REF;
    my $skip_recalculating_yield_tum = $SKIP_RECALCULATING_YIELD_TUM;

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
        addSampleToJsonData(
            \%json_data,
            $submission,
            $barcode,
            $name,
            'ref',
            $q30,
            $yield,
            $use_existing_ref,
            $skip_recalculating_yield_ref,
        );
    }
    elsif ( $analysis eq 'FASTQ' ){
        my $set = join("_", $date, $submission, $barcode, $name );
        sayInfo("  SET: $set");
        $json_data{ 'set_name' } = "$set";
        $json_data{ 'entity' } = "$entity";

        $json_data{ 'fastq_portal' } = JSON::true;
        addSampleToJsonData(
            \%json_data,
            $submission,
            $barcode,
            $name,
            'ref',
            $q30,
            $yield,
            $use_existing_ref,
            $skip_recalculating_yield_ref,
        );
    }
    elsif( $analysis eq 'SingleAnalysis' ){
        my $set = join( "_", $date, $submission, $barcode, $name );
        sayInfo("  SET: $set");
        $json_data{ 'ini' } = "$GER_INI";
        $json_data{ 'set_name' } = "$set";
        $json_data{ 'entity' } = "$entity";
        addSampleToJsonData(
            \%json_data,
            $submission,
            $barcode,
            $name,
            'ref',
            $q30,
            $yield,
            $use_existing_ref,
            $skip_recalculating_yield_ref,
        );
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
        addSampleToJsonData(
            \%json_data,
            $submission,
            $barcode,
            $name,
            'tumor-rna',
            $q30,
            $yield,
            $use_existing_ref,
            $skip_recalculating_yield_ref,
        );
    }
    elsif ( $analysis eq 'Somatic_T' ){
        
        my $ref_obj;
        my $ini = $SOM_INI;
        my $needs_shallow = getValueByKey( $sample, 'shallowseq' ); # 0 | 1
        #my $other_ref = getValueByKey( $sample, 'other_ref' ); # Yes | No
        
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
        my $patient_ref = getValueByKey( $ref_obj, 'patient' );
        my $yield_ref = getValueByKey( $ref_obj, 'yield' );
        my $submission_ref = getValueByKey( $ref_obj, 'submission' );
        my $lab_status_ref = getValueByKey( $ref_obj, 'lab_status' );

        if ( $lab_status_ref ne "Finished" and $lab_status_ref ne "finished"  and $lab_status ne "In process" ){
            sayWarn("  Check JSON for sample $name, since ref lab status is unexpected: lab_status_ref='$lab_status_ref'")
        }

        $yield_ref = $yield_ref == 0 ? 1 : $yield_ref * $YIELD_F;
        my $set = join( "_", $date, $submission, $barcode_ref, $barcode, $patient );

        ## adjust content in case of ShallowSeq
        if ( $needs_shallow and not $IGNORE_SHALLOWSEQ ){
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
                sayWarn("  RESULT: SKIPPING because ShallowSeq runs with status $run_status found for $match_string");
                return "NoJsonMade_DeletedShallowSeqRunFound";
            }
        }

        my $ref_api_status = `hmf_api_get 'samples?barcode=$barcode_ref' | jq '.[0].status' | tr -d '"\n'`;
        if ( $patient ne $patient_ref ){
            ## add suffix to ref barcode and use tumor submission in case ref is needed from other existing patientId
            my $new_name_ref = $patient . 'R';
            my $new_barcode_ref = getCorrectBarcodeWithSuffixForRefSampleName(
                $new_name_ref,
                $barcode_ref,
                $name,
                $barcode,
                \@warn_msg,
                $date,
                "OtherRef flag is set in LIMS.",
            );
            $name_ref = $new_name_ref;
            $barcode_ref = $new_barcode_ref;
            $submission_ref = $submission;
            $skip_recalculating_yield_ref = 1;
        } elsif ( $ref_api_status eq "Deleted" ){
            ## add suffix to ref barcode in case ref fastq has already been deleted
            my $new_barcode_ref = getCorrectBarcodeWithSuffixForRefSampleName(
                $name_ref,
                $barcode_ref,
                $name,
                $barcode,
                \@warn_msg,
                $date,
                "REF sample is 'Deleted'.",
            );
            $barcode_ref = $new_barcode_ref;
            $skip_recalculating_yield_ref = 1;
        }
 
        ## check if barcode already exists in HMF API
        my $tum_exists = `hmf_api_get 'samples?barcode=$barcode' | jq 'map(select(.status != "Unregistered")) | length'`;
        my $ref_exists = `hmf_api_get 'samples?barcode=$barcode_ref' | jq 'map(select(.status != "Unregistered")) | length'`;
        chomp($tum_exists);
        chomp($ref_exists);

        if ( $tum_exists ne "0" ){
            sayInfo("  TUM barcode ($barcode) already exists in HMF API (so will add use_existing flag)");
            $use_existing_tum = 1;
        }
        if ( $ref_exists ne "0" ){
            sayInfo("  REF barcode ($barcode_ref) already exists in HMF API (so will add use_existing flag)");
            $use_existing_ref = 1;
        }

        sayInfo("  Set name constructed to $set");
        $json_data{ 'ini' } = "$ini";
        $json_data{ 'set_name' } = "$set";
        $json_data{ 'entity' } = "$entity";
        $json_data{ 'priority' } = $priority;
        addSampleToJsonData(
            \%json_data,
            $submission_ref,
            $barcode_ref,
            $name_ref,
            'ref',
            $q30,
            $yield_ref,
            $use_existing_ref,
            $skip_recalculating_yield_ref,
        );
        addSampleToJsonData(
            \%json_data,
            $submission,
            $barcode,
            $name,
            'tumor',
            $q30,
            $yield,
            $use_existing_tum,
            $skip_recalculating_yield_tum,
        );
    }
    elsif ( $analysis eq 'Targeted_Tumor_Only' ){
        my $ini = $TAR_INI;
        my $set = join( "_", $date, $submission, $barcode, $patient );
        my $exists_in_api = `hmf_api_get 'samples?barcode=$barcode' | jq 'map(select(.status != "Unregistered")) | length'`;
        chomp($exists_in_api);

        if ( $exists_in_api ne "0" ){
            sayInfo("  Sample barcode ($barcode) already exists in HMF API (so will add use_existing flag)");
            $use_existing_tum = 1;
        }

        sayInfo("  Set name constructed to $set");
        $json_data{ 'ini' } = "$ini";
        $json_data{ 'set_name' } = "$set";
        $json_data{ 'entity' } = "$entity";
        $json_data{ 'priority' } = $priority;
        addSampleToJsonData(
            \%json_data,
            $submission,
            $barcode,
            $name,
            'tumor',
            $q30,
            $yield,
            $use_existing_tum,
            $skip_recalculating_yield_tum,
        );

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
    my $existing_count = `hmf_api_get 'runs?set_name_contains=$setname_wo_date' | jq 'map(select(.status != "Invalidated")) | length' | tr -d '"\n'`;

    if ( $existing_count ne "0" ){
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

sub getCorrectBarcodeWithSuffixForRefSampleName{
    my ($name_ref, $barcode_ref, $name, $barcode, $warn_msg, $date, $change_reason) = @_;
    my $ready_new_ref_sample_count = `hmf_api_get 'samples?name=$name_ref' | jq 'map(select(.status == "Ready")) | length' | tr -d '"\n'`;
    my $new_barcode_ref;
    my $new_warn_msg;
    if ( $ready_new_ref_sample_count eq "1" ){
        $new_barcode_ref = `hmf_api_get 'samples?name=$name_ref' | jq 'map(select(.status == "Ready"))[0].barcode' | tr -d '"\n'`;
        $new_warn_msg = "DOUBLE CHECK JSON for $barcode ($name): " .
            "$change_reason Reusing existing 'Ready' REF barcode ($new_barcode_ref) to replace ($barcode_ref)";
    } elsif ( $ready_new_ref_sample_count eq "0" ) {
        $new_barcode_ref = $barcode_ref . "-c2f" . $date;
        $new_warn_msg = "DOUBLE CHECK JSON for $barcode ($name): " .
            "$change_reason Adding suffix to create new REF barcode ($new_barcode_ref)";
    } else {
        $new_barcode_ref = `hmf_api_get 'samples?name=$name_ref' | jq 'map(select(.status == "Ready"))[0].barcode' | tr -d '"\n'`;
        $new_warn_msg = "DOUBLE CHECK JSON for $barcode ($name): " .
            "$change_reason Out of multiple choices, randomly chose existing 'Ready' REF barcode ($new_barcode_ref) to replace ($barcode_ref)";
    }
    push( @$warn_msg, $new_warn_msg );
    if ( $new_barcode_ref !~ /^$barcode_ref.*$/ ) {
        push( @$warn_msg, "Replacement barcode ($new_barcode_ref) does not start with ($barcode_ref)");
    }
    return $new_barcode_ref;
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
    my ($store, $submission, $barcode, $name, $type, $q30, $yield, $use_existing, $skip_recalculating_yield) = @_;
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
    if ( $skip_recalculating_yield ){
        $tmp{ 'skip_recalculate' } = JSON::true;
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

sub removeNewlineCharacters {
    my $text = shift;
    $text =~ s/[\r\n]+//g;
    return $text;
}

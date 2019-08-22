#!/usr/bin/perl
use strict;
use warnings;
use 5.010.000;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use JSON;
use File::Slurp;

my $SBP_API_SCRIPT = "query_sbp_api";
my $SBP_JSON_ROOT = "/data/lims/sbpfiles";
my $SBP_SETS_JSON = "$SBP_JSON_ROOT/sets.json";
my $SBP_RUNS_JSON = "$SBP_JSON_ROOT/runs.json";
my $SBP_SAMP_JSON = "$SBP_JSON_ROOT/samples.json";

my $NA_CHAR = "";
my $INPUT_DELIM_OVERVIEW = "\t";
my $INPUT_DELIM_CLINICAL = "\t";
my $OUTPUT_DELIM = "\t";

my %REQUIRED_FIELDS_OVERVIEW = (
    'Label' => 'hmfLabel',
    'SampleID' => 'hmfSampleId',
    'ExternalID' => 'hmfExternalId',
    'SubmissionID' => 'hmfSubmissionId',
    'Set Name' => 'hmfSetName',
    'Status' => 'hmfStatus',
);

my %REQUIRED_FIELDS_CLINICAL = (
    "puritySampleId" => "puritySampleId",
    "purityQC" => "purityQC",
    "purityStatus" => "purityStatus",
    "tumorPurity" => "tumorPurity",
    "sampleId" => "sampleId",
    "sampleArrivalDate" => "sampleArrivalDate",
    "patientId" => "patientId",
    "hospital" => "hospital",
    "gender" => "gender",
    "birthYear" => "birthYear",
    "hasSystemicPreTreatment" => "hasSystemicPreTreatment",
    "hasRadiotherapyPreTreatment" => "hasRadiotherapyPreTreatment",
    "preTreatments" => "preTreatments",
    "informedConsentDate" => "informedConsentDate",
    "biopsyDate" => "biopsyDate",
    "treatmentStartDate" => "treatmentStartDate",
    "treatmentEndDate" => "treatmentEndDate",
    "responseDate" => "responseDate",
    "deathDate" => "deathDate",
    "primaryTumorLocation" => "primaryTumorLocation",
    "cancerSubtype" => "cancerSubtype",
    "biopsySite" => "biopsySite",
    "biopsyLocation" => "biopsyLocation",
    "treatmentGiven" => "treatmentGiven",
    "radiotherapyGiven" => "radiotherapyGiven",
    "treatment" => "treatment",
    "treatmentType" => "treatmentType",
    "responseMeasured" => "responseMeasured",
    "firstResponse" => "firstResponse",
    "allPostTreatments" => "allPostTreatments",
    "allPostTreatmentTypes" => "allPostTreatmentTypes",
);

my %REQUIRED_FIELDS_DRUPCPCT = (
    "patientId" => "patientId",
    "CPCT2YN" => "paticipatedInCpct",
    "CPCTCNT" => "cpctCenterNumber",
    "CPCTPN" => "cpctPatientNumber",
);

my @PRINT_FIELDS_OVERVIEW = qw( sbpSetName sbpSetId sbpBucket sbpStatus hmfStatus hmfLabel );

#my @PRINT_FIELDS_CLINICAL = qw(
#  puritySampleId purityQC purityStatus tumorPurity sampleId patientId primaryTumorLocation cancerSubtype gender birthYear deathDate biopsyDate
#  biopsySite biopsyLocation hasSystemicPreTreatment hasRadiotherapyPreTreatment treatmentGiven treatmentStartDate treatmentEndDate
#  treatment treatmentType responseDate responseMeasured firstResponse preTreatments allPostTreatments allPostTreatmentTypes
#);

my @PRINT_FIELDS_CLINICAL = qw(
    puritySampleId purityQC purityStatus tumorPurity sampleId patientId primaryTumorLocation cancerSubtype biopsyDate biopsySite biopsyLocation
    gender birthYear deathDate hasSystemicPreTreatment hasRadiotherapyPreTreatment treatmentGiven treatmentStartDate treatmentEndDate
    treatment treatmentType responseDate responseMeasured firstResponse preTreatments allPostTreatments allPostTreatmentTypes
);

my $SORT_COLUMN = 'hmfExternalId'; # no clinical field possible (is 1 level deeper)

my $samp_in_file;
my $clin_in_file;
my $drup_in_file;
my $output_file;
my $skip_qc;
my $debug;

GetOptions (
    "sample_overview_file|s=s" => \$samp_in_file,
    "db_clinical_file|c=s" => \$clin_in_file,
    "db_drupcpct_file|d=s" => \$drup_in_file,
    "output_file|o=s" => \$output_file,
    "noqc" => \$skip_qc,
    "debug" => \$debug,
) or die("Error in command line arguments\n");

## input checks
die "[EXIT] Run with -s <SampleOverview.tsv> -c <Clinical.tsv> -d <drupcpct.tsv> -o <Output.tsv> \n" unless $samp_in_file and $clin_in_file and $drup_in_file;
die "[EXIT] Input file not found ($samp_in_file): $!\n" unless -f $samp_in_file;
die "[EXIT] Input file not found ($clin_in_file): $!\n" unless -f $clin_in_file;
die "[EXIT] Input file not found ($drup_in_file): $!\n" unless -f $drup_in_file;

## gather inputs
say "[INFO] Parsing SBP json to get runs info";
my $runs_db = parseSBPjson( $SBP_API_SCRIPT, 'runs' );
say "[INFO] Parsing $clin_in_file";
my $clin_db = parseSqlQueryFile( $clin_in_file, $INPUT_DELIM_CLINICAL, \%REQUIRED_FIELDS_CLINICAL, 'puritySampleId' );
#print Dumper $clin_db;
say "[INFO] Parsing $drup_in_file";
my $drup_db = parseSqlQueryFile( $drup_in_file, $INPUT_DELIM_CLINICAL, \%REQUIRED_FIELDS_DRUPCPCT, 'patientId' );
say "[INFO] Parsing $samp_in_file";
my $data_db = parseSampleOverview( $samp_in_file, $INPUT_DELIM_OVERVIEW, \%REQUIRED_FIELDS_OVERVIEW );

## make selection and merge
$data_db = selectAllDoneCPCTandDRUP( $data_db );

#say "[DEBUG] After selection before merge";
#print Dumper $data_db;
#<>;

$data_db = mergeData( $data_db, $clin_db, $runs_db, $drup_db );

## output
if ( $output_file ){
    open my $fh, ">", $output_file or die "[EXIT] Unable to open output file ($output_file):$!\n";
    printSelection( $data_db, $fh );
    close $fh;
    say "[INFO] Output written to: $output_file";
}else{
    printSelection( $data_db, *STDOUT );
}

## potential QC
performQualCheck( $data_db ) unless $skip_qc;

## -----

sub parseSBPjson{
    my ($script, $type) = @_;
    my $cmd = "$script -type $type -json";
    my $txt = `$cmd`;
    my $objects = decode_json( $txt );
    return $objects;
}

sub mergeData{
    my ($data, $clin, $runs, $drup) = @_;

    #print Dumper($data);
    #<>;

    my @merged = @$data;
    foreach my $obj ( @merged ){
        my $biopsy = $obj->{ hmfExternalId };

        my $sbp_obj    = selectSBPrun( $obj, $runs );
        my $sbp_id     = $sbp_obj->{ 'set' }{ 'id' } || $NA_CHAR;
        my $sbp_name   = $sbp_obj->{ 'name' } || $NA_CHAR;
        my $sbp_bucket = $sbp_obj->{ 'bucket' } || $NA_CHAR;
        my $sbp_status = $sbp_obj->{ 'status' } || $NA_CHAR;

        $obj->{ 'sbpSetId' }   = $sbp_id;
        $obj->{ 'sbpSetName' } = $sbp_name;
        $obj->{ 'sbpBucket' }  = $sbp_bucket;
        $obj->{ 'sbpStatus' }  = $sbp_status;

        ## in case of DRUP try to match with CPCT
        #if ( $biopsy =~ /^(DRUP\d{8})T[IVX]*/ ){
        #    my $patient = $1;
        #    $patient = findCpctIdForDrupPatient( $patient, $drup );
        #    $clin->{ $biopsy }{ 'patientId' } = $patient;
        #}

        if ( exists $clin->{ $biopsy } ){
            my %clin_info = %{$clin_db->{ $biopsy }};
            $obj->{ 'clin_info' } = \%clin_info;
        }else{
            warn "[WARN] No clin info found for biopsy ($biopsy)\n";
        }
    }
    return \@merged;
}

sub findCpctIdForDrupPatient{
    my ($drupID, $drupDB) = @_;
    if ( exists $drupDB->{ $drupID } ){

        #my $participated = $drupDB->{ $drupID }{ paticipatedInCpct } eq 'Yes';
        my $cpctCenID = $drupDB->{ $drupID }{ 'cpctCenterNumber' };
        my $cpctPatID = $drupDB->{ $drupID }{ 'cpctPatientNumber' };

        if ( $cpctCenID =~ /^\d{2}$/ and $cpctPatID =~ /^\d{4}$/ ){
            my $cpctID = join( "", "CPCT02", $cpctCenID, $cpctPatID);
            say "[INFO] DRUP/CPCT match found: changing patientID $drupID into $cpctID";
            return $cpctID;
        }
    }
    return $drupID;
}

sub performQualCheck{
    my ($objects) = @_;
    my $idx = 0;

    foreach my $obj ( @$objects ){
        my $biopsy = $obj->{ hmfExternalId };
        my $hmfset = $obj->{ hmfSetName };
        my $sbpset = $obj->{ sbpSetName };
        my $bucket = $obj->{ sbpBucket };
        if ( ($hmfset ne $sbpset) and ($hmfset ne $NA_CHAR) and ($sbpset ne $NA_CHAR) ){
            #warn "[WARN] HMF vs SBP sets mismatch ($biopsy: \"$hmfset\" vs \"$sbpset\")\n";
            say "[WARN] HMF vs SBP sets mismatch ($biopsy: \"$hmfset\" vs \"$sbpset\")";
        }
        if ( $bucket !~ /^hmf/ ){
            warn "[WARN] Found unexpected bucket ($bucket) for set ($sbpset) of biopsy $biopsy\n";
        }
    }
}

sub selectSBPrun{
    my ($hmfobj, $sbp_runs) = @_;
    my $biopsy = $hmfobj->{ hmfExternalId };

    unless ( exists $hmfobj->{ hmfSetName } ){
        print Dumper $hmfobj;
        die "[EXIT] No set defined?\n"
    }

    my @selection = ();
    foreach my $obj ( @$sbp_runs ){

        ## a run should apply to all rules
        next unless exists $obj->{ tumor_sample };
        next unless exists $obj->{ pipeline };
        next unless exists $obj->{ bucket };

        next unless $obj->{ tumor_sample } eq $biopsy;
        #next unless $obj->{ status } eq 'Validated';
        next unless $obj->{ ini } =~ /Somatic/;
        next unless defined $obj->{ bucket };
        next if $obj->{ bucket } eq 'hmf_experiments';
        my $pipeline = $obj->{ pipeline };
        #my $setname = $obj->{ hmfSetName };
        #next unless $obj->{ bucket } eq 'hmf_experiments';

        if ( $pipeline =~ /v(\d+)\.(\d+)/ ){
            my $major = $1;
            my $minor = $2;
            next if $major < 3;
        }
        push( @selection, $obj );
    }

    my $count = scalar(@selection);
    if ( $count > 1 ){
        return( $selection[-1] );
    }
    elsif( $count == 1 ){
        return( $selection[0] );
    }
    else{
        warn "[WARN] Biopsy without any Validated run ($biopsy) at SBP!\n";
        return undef;
    }
}

sub printSelection{
    my ($objects, $fh) = @_;

    ## make order more readable
    @$objects = sort {
        $a->{$SORT_COLUMN} cmp $b->{$SORT_COLUMN}
    } @$objects;

    my $count = scalar @$objects;

    say "[INFO] Final selection count: $count samples";
    say $fh "#".join( $OUTPUT_DELIM, @PRINT_FIELDS_OVERVIEW, @PRINT_FIELDS_CLINICAL );

    foreach my $obj ( @$objects ){
        my @out = map( $obj->{$_} || $NA_CHAR, @PRINT_FIELDS_OVERVIEW );
        if ( exists $obj->{ 'clin_info' } ){
            push( @out, map( $obj->{ 'clin_info' }{$_} || $NA_CHAR, @PRINT_FIELDS_CLINICAL ) );
        }
        say $fh join( $OUTPUT_DELIM, @out );
    }
}

sub selectAllDoneCPCTandDRUP{
    my ($objects) = @_;
    my @out = ();

    foreach my $obj ( @$objects ){
        next unless $obj->{ hmfStatus } eq 'Done';
        next unless $obj->{ hmfLabel } =~ m/^(CPCT|DRUP)$/;
        next unless $obj->{ hmfExternalId } =~ /^(CPCT|DRUP)\d{8}T/;
        push( @out, $obj );
    }

    return \@out;
}

sub parseSqlQueryFile{
    my ($file, $input_delim, $required_fields, $store_field_name) = @_;
    my %out = ();

    say "[INFO] Opening file: $file";
    open my $fh, "<", $file or die "[EXIT] Unable to open file ($file): $!\n";
    my $header = findAndCheckHeader( $fh, $input_delim, $required_fields );

    while ( <$fh> ){
        chomp;
        next if $_ =~ /^#/;
        next unless /\S/;

        my %tmp = ();
        my @values = split( "\t", $_ );

        foreach my $field_name ( @$header ){
            my $field_value = shift @values;
            $field_value = $NA_CHAR if ( not defined $field_value or $field_value eq "NULL");
            if ( exists $required_fields->{ $field_name } ){
                my $store_name = $required_fields->{ $field_name };
                $tmp{ $store_name } = $field_value;
            }
        }

        ## sanity check
        if ( not exists $tmp{ $store_field_name } or $tmp{ $store_field_name } eq '' ){
            say "[WARNING] has no \"$store_field_name\" in input file \"$file\" for following object";
            print Dumper( \%tmp );
        }

        ## all ok: continue by storing at unique id
        my $id = $tmp{ $store_field_name };
        die "[EXIT] Key ($id) already exists in db\n" if exists $out{ $id };
        $out{ $id } = \%tmp;
    }
    close $fh;

    return \%out;
}

sub parseSampleOverview{
    my ($file, $input_delim, $required_fields) = @_;
    my @out = ();

    say "[INFO] Opening file: $file";
    open my $fh, "<", $file or die "[EXIT] Unable to open file ($file): $!\n";
    my $header = findAndCheckHeader( $fh, $input_delim, $required_fields );

    while ( <$fh> ){
        chomp;
        next if $_ =~ /^#/;
        $_ =~ s/\r//g; # remove any windowns ^M carriage return
        my %tmp = ();
        my @values = split( $input_delim, $_ );
        foreach my $field_name ( @$header ){
            my $field_value = shift @values;
            if ( exists $required_fields->{ $field_name } ){
                my $store_name = $required_fields->{ $field_name };
                $tmp{ $store_name } = $field_value;
            }
        }
        push( @out, \%tmp );
    }

    close $fh;
    return \@out;
}

sub findAndCheckHeader{
    my ($fh, $delim, $required_fields_dict) = @_;
    my @requir_fields = keys %$required_fields_dict;
    my $requir_string = join( ", ", @requir_fields );

    while ( <$fh> ){
        chomp;
        next if $_ =~ /^##/;
        $_ =~ s/\r//g; # remove any windowns ^M carriage return
        my $header_line = $_;
        my @header_fields = split( "\t", $header_line );
        my $header_print = substr( join( ',', @header_fields ), 0, 2000);
        my $header_ok = 1;

        say "[INFO] Header line found: checking required fields";

        if ( $debug ){
            say "[INFO] Header fields ($header_print)";
            say "[INFO] Required fields ($requir_string)";
        }

        foreach my $req_field ( @requir_fields ){
            my $req_field_found = 0;
            foreach my $head_field ( @header_fields ){
                $req_field_found = 1 if $req_field eq $head_field;
            }
            unless ( $req_field_found ){
                warn "[WARN] Header field \"$req_field\" not found in header\n";
                $header_ok = 0;
            }
        }

        if ( $header_ok ){
            return \@header_fields;
        }
        else{
            die "[EXIT] Something wrong with header\n";
        }
    }
    die "[EXIT] Reached end of findAndCheckHeader: should not happen\n";
}

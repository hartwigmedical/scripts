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
use POSIX qw(strftime);
use 5.01000;
$| = 1; # Disable stdout buffering

use constant EMPTY => q{ };
use constant NACHAR => 'NA';
use constant DATA_SOURCE_FIELD => "data_source";
use constant ONCOACT_ENTITY => "CORE_01";

# Fields that will actively be set to boolean for json output
use constant BOOLEAN_FIELDS => qw(shallowseq report_germline report_viral report_pgx);
# Fields that will actively be set to integer for json output
use constant INTEGER_FIELDS => qw(yield q30);
# Fields that will be required to exist in input if checked
my %WARN_IF_ABSENT_IN_LAMA_FIELDS = (_id=>1, status=>1);

my $SCRIPT  = basename($0);
my %opt = (
    'center_tsv' => '/data/ops/lims/prod/center2entity.tsv'
);
my $HELP_TEXT = <<"HELP";

  Description
    Converts the sample information from various LIMS sources (Excel forms and LAMA)
    into one coherent JSON output file.

  Usage
    $SCRIPT -lims_dir /data/ops/lims/pilot -out_json /data/ops/lims/pilot/lims.json

  Required:
    -lims_dir <str>  Path to input dir (eg /data/ops/lims/pilot)
    -out_json <str>  Path to output json (eg /data/tmp/lims.json)

  Optional:
    -center_tsv <str>  Path to center dictionary input tsv [$opt{center_tsv}]

HELP

# Get input and setup all paths
GetOptions (
    "lims_dir=s"   => \$opt{lims_dir},
    "out_json=s"   => \$opt{out_json},
    "center_tsv=s" => \$opt{center_tsv},
    "debug"        => \$opt{debug},
    "help|h"       => \$opt{help},
) or die("Error in command line arguments\n");

die $HELP_TEXT if $opt{help};
die $HELP_TEXT unless $opt{lims_dir};
die $HELP_TEXT unless $opt{out_json};

my $CNTR_TSV = $opt{center_tsv};
my $LIMS_DIR = $opt{lims_dir};
my $JSON_OUT = $opt{out_json};
my $LATEST_DIR = $LIMS_DIR . "/lab_files/latest";

# Current LIMS files
my $FOR_001_SUBM_TSV = $LATEST_DIR . '/for001_submissions.tsv';
my $FOR_001_CONT_TSV = $LATEST_DIR . '/for001_contacts.tsv';
my $FOR_001_SAMP_TSV = $LATEST_DIR . '/for001_samples.tsv';
my $FOR_002_PROC_TSV = $LATEST_DIR . '/for002_processing.tsv';

# Current LAMA files
my $LAMA_ISOLATION_JSON = $LATEST_DIR . '/Isolations.json';
my $LAMA_PATIENT_JSON = $LATEST_DIR . '/Patients.json';
my $LAMA_LIBRARYPREP_JSON = $LATEST_DIR . '/LibraryPreparations.json';
my $LAMA_SAMPLESTATUS_JSON = $LATEST_DIR . '/Statuses.json';
my $LAMA_CONTRACTS_JSON = $LATEST_DIR . '/Contracts.json';

# Files from previous years
my $SUBM_TSV_2023 = $LATEST_DIR . '/2023_for001_submissions.tsv';
my $SAMP_TSV_2023 = $LATEST_DIR . '/2023_for001_samples.tsv';
my $PROC_TSV_2023 = $LATEST_DIR . '/2023_for002_processing.tsv';

my $SUBM_TSV_2022 = $LATEST_DIR . '/2022_for001_submissions.tsv';
my $SAMP_TSV_2022 = $LATEST_DIR . '/2022_for001_samples.tsv';
my $PROC_TSV_2022 = $LATEST_DIR . '/2022_for002_processing.tsv';

my $SUBM_TSV_2021 = $LATEST_DIR . '/2021_for001_submissions.tsv';
my $SAMP_TSV_2021 = $LATEST_DIR . '/2021_for001_samples.tsv';
my $PROC_TSV_2021 = $LATEST_DIR . '/2021_for002_processing.tsv';

my $SUBM_TSV_2020 = $LATEST_DIR . '/2020_for001_submissions.tsv';
my $SAMP_TSV_2020 = $LATEST_DIR . '/2020_for001_samples.tsv';
my $PROC_TSV_2020 = $LATEST_DIR . '/2020_for002_processing.tsv';

my $SUBM_TSV_2019 = $LATEST_DIR . '/2019_for001_submissions.tsv';
my $SAMP_TSV_2019 = $LATEST_DIR . '/2019_for001_samples.tsv';
my $PROC_TSV_2019 = $LATEST_DIR . '/2019_for002_processing.tsv';

my $SUBM_TSV_2018 = $LATEST_DIR . '/2018_subm';
my $SAMP_TSV_2018 = $LATEST_DIR . '/2018_samp';
my $PROC_TSV_2018 = $LATEST_DIR . '/2018_proc';

my $PROC_TSV_2017 = $LATEST_DIR . '/2017_proc';
my $LIMS_JSN_2017 = $LATEST_DIR . '/2017_lims.json';

my @ALL_INPUT_FILES = (
    $LAMA_ISOLATION_JSON, $LAMA_PATIENT_JSON, $LAMA_LIBRARYPREP_JSON, $LAMA_SAMPLESTATUS_JSON,
    $FOR_001_SUBM_TSV, $FOR_001_SAMP_TSV, $FOR_002_PROC_TSV, $FOR_001_CONT_TSV,
    $SUBM_TSV_2023, $SAMP_TSV_2023, $PROC_TSV_2023,
    $SUBM_TSV_2022, $SAMP_TSV_2022, $PROC_TSV_2022,
    $SUBM_TSV_2021, $SAMP_TSV_2021, $PROC_TSV_2021,
    $SUBM_TSV_2020, $SAMP_TSV_2020, $PROC_TSV_2020,
    $SUBM_TSV_2019, $SAMP_TSV_2019, $PROC_TSV_2019,
    $SUBM_TSV_2018, $SAMP_TSV_2018, $PROC_TSV_2018,
    $PROC_TSV_2017, $LIMS_JSN_2017,
    $CNTR_TSV
);

foreach ( $LIMS_DIR ){
    die "[ERROR] Input dir does not exist ($_)\n" unless -e $_;
    die "[ERROR] Input dir is not a directory ($_)\n" unless -d $_;
}
foreach ( @ALL_INPUT_FILES ){
    die "[ERROR] Input file does not exist ($_)\n" unless -f $_;
}
foreach ( $JSON_OUT ){
    die "[ERROR] Output file exists and is not writable ($_)\n" if ( -f $_ and not -w $_ );
}

sayInfo(sprintf "Starting with %s", $SCRIPT);
sayInfo(sprintf "Input dir: %s", $LATEST_DIR);

my $name_dict = getFieldNameTranslations();
my $cntr_dict = parseDictFile( $CNTR_TSV, 'center2centername' );

my $proc_objs = {}; # will contain objects from InProcess sheet
my $subm_objs = {}; # will contain objects from Received-Samples shipments sheet
my $cont_objs = {}; # will contain objects from Received-Samples contact sheet
my $samp_objs = {}; # will contain objects from Received-Samples samples sheet
my $lims_objs = {}; # will contain all sample objects

$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2017, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2018, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2019, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2020, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2021, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2022, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $PROC_TSV_2023, "\t" );
$proc_objs = parseTsvCsv( $proc_objs, $name_dict->{'PROC_CURR'}, 'sample_id',  0, $FOR_002_PROC_TSV, "\t" );

$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_2018'}, 'submission', 0, $SUBM_TSV_2018, "\t" );
$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_2019'}, 'submission', 0, $SUBM_TSV_2019, "\t" );
$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_2020'}, 'submission', 0, $SUBM_TSV_2020, "\t" );
$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_2021'}, 'submission', 0, $SUBM_TSV_2021, "\t" );
$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_2022'}, 'submission', 0, $SUBM_TSV_2022, "\t" );
$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_2023'}, 'submission', 0, $SUBM_TSV_2023, "\t" );
$subm_objs = parseTsvCsv( $subm_objs, $name_dict->{'SUBM_CURR'}, 'submission', 0, $FOR_001_SUBM_TSV, "\t" );
$cont_objs = parseTsvCsv( $cont_objs, $name_dict->{'CONT_CURR'}, 'group_id',   1, $FOR_001_CONT_TSV, "\t" );

$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_2018'}, 'sample_id',  1, $SAMP_TSV_2018, "\t" );
$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_2019'}, 'sample_id',  1, $SAMP_TSV_2019, "\t" );
$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_2020'}, 'sample_id',  1, $SAMP_TSV_2020, "\t" );
$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_2021'}, 'sample_id',  1, $SAMP_TSV_2021, "\t" );
$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_2022'}, 'sample_id',  1, $SAMP_TSV_2022, "\t" );
$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_2023'}, 'sample_id',  1, $SAMP_TSV_2023, "\t" );
$samp_objs = parseTsvCsv( $samp_objs, $name_dict->{'SAMP_CURR'}, 'sample_id',  1, $FOR_001_SAMP_TSV, "\t" );

my $lama_status = parseLamaSampleStatus($LAMA_SAMPLESTATUS_JSON);
my $lama_isolation = parseLamaIsolation($LAMA_ISOLATION_JSON);
my $lama_sample = parseLamaPatients($LAMA_PATIENT_JSON);
my $lama_prep = parseLamaLibraryPreps($LAMA_LIBRARYPREP_JSON);
my $lama_contracts = parseLamaContracts($LAMA_CONTRACTS_JSON);

checkContactInfo( $cont_objs );

$subm_objs = addContactInfoToSubmissions( $subm_objs, $cont_objs );
$lims_objs = addExcelSamples( $lims_objs, $samp_objs, $subm_objs );
$lims_objs = addLamaSamples( $lims_objs, $lama_status, $lama_sample, $lama_isolation, $lama_prep, $lama_contracts, $subm_objs, $cntr_dict );
$lims_objs = addLabSopString( $lims_objs, $proc_objs );

fixAddedDateFields( $lims_objs );
checkDrupStage3Info( $subm_objs, $lims_objs );
printLimsToJson( $lims_objs, $subm_objs, $cont_objs, $JSON_OUT );
sayInfo(sprintf "Finished with %s", $SCRIPT);

sub addLamaSamples{
    my ($lims, $statuses, $samples, $isolations, $preps, $contracts, $submissions, $centers_dict) = @_;
    my %store = %{$lims};
    my %dna_reference_samples_by_name = ();
    sayInfo("  Adding LAMA samples");

    SAMPLESTATUS: while (my ($isolate_barcode, $object) = each %$statuses) {
        my %sample_to_store = %{$object};

        if (not exists $sample_to_store{sample_barcode}) {
            sayWarn(sprintf "NOTIFY: SKIPPING LAMA sample because sample_barcode not present in status object [%s]", $isolate_barcode);
            next;
        }

        my $sample_barcode = $sample_to_store{sample_barcode};

        # add fields for backwards compatibility
        $sample_to_store{received_sample_id} = $sample_barcode;
        $sample_to_store{report_germline_level} = NACHAR;
        $sample_to_store{report_viral} = JSON::XS::true;
        $sample_to_store{report_pgx} = JSON::XS::true;
        $sample_to_store{shallowseq} = JSON::XS::false if not exists $sample_to_store{shallowseq};

        # adding sample info
        if (exists $samples->{$sample_barcode}) {
            addRecordFieldsToTargetRecord($samples->{$sample_barcode}, \%sample_to_store, "merge of sample info for $isolate_barcode");
            # retaining the roman naming for older samples for the time being (can be removed once anonymization project is finished)
            if (exists $samples->{$sample_barcode}{legacy_sample_name}){
                $sample_to_store{sample_name} = $samples->{$sample_barcode}{legacy_sample_name};
            }
        }

        # adding isolate info
        if (exists $isolations->{$isolate_barcode}) {
            addRecordFieldsToTargetRecord($isolations->{$isolate_barcode}, \%sample_to_store, "merge of isolate info for $isolate_barcode");
        }

        # adding prep info
        if (exists $preps->{$isolate_barcode}) {
            addRecordFieldsToTargetRecord($preps->{$isolate_barcode}, \%sample_to_store, "merge of prep info for $isolate_barcode");
        }

        # by now we require certain fields to be present
        my @required_fields = qw(sample_name sample_id contract_code);
        foreach my $field (@required_fields){
            unless (exists $sample_to_store{$field}){
                sayWarn(sprintf "NOTIFY: SKIPPING LAMA status because field [%s] not present [%s,%s]", $field, $sample_barcode, $isolate_barcode);
                next SAMPLESTATUS;
            }
        }

        my $sample_name = $sample_to_store{sample_name};
        my $sample_id = $sample_to_store{sample_id};
        my $contract_code = $sample_to_store{contract_code};
        my $sample_print_info = "$sample_id,$sample_name,$contract_code";

        # adding contract info
        if (exists $contracts->{$contract_code}) {
            addRecordFieldsToTargetRecord($contracts->{$contract_code}, \%sample_to_store, "merge of prep info for $isolate_barcode");
        }else{
            sayWarn("Unable to locate contract [$contract_code] for sample [$sample_print_info]");
        }

        # adding cohort info for backwards compatibility (the concept "cohort" does no longer exist since LAMA v2)
        if (exists $sample_to_store{contract_display_name} and exists $sample_to_store{contract_sample_id_start}){
            $sample_to_store{'contract_sample_id_start'} = $sample_to_store{'contract_sample_id_start'};
            $sample_to_store{'cohort'} = contractToCohort($sample_to_store{'contract_display_name'});
            $sample_to_store{'cohort_code'} = $sample_to_store{'contract_sample_id_start'};
        }else{
            $sample_to_store{'cohort'} = NACHAR;
            $sample_to_store{'cohort_code'} = NACHAR;
        }

        my ($patient_id, $study, $center, $tum_or_ref);
        my $name_regex = '^((CPCT|DRUP|WIDE|ACTN|CORE|SHRP|GAYA|TARG|OMIC|OPTC|GLOW|QUAL)[0-9A-Z]{2}([0-9A-Z]{2})\d{4})(T|R){1}';
        if ($sample_name =~ /$name_regex/ms) {
            ($patient_id, $study, $center, $tum_or_ref) = ($1, $2, $3, $4);
            $sample_to_store{label} = $study;
        }
        elsif ($sample_name =~ /^XXXXXX[0-9A-Z]{2}\d{4}R/ms) {
            # These are reference samples with a temporary name until tumor sample arrives
            next;
        }
        else {
            sayWarn("NOTIFY: Unable to use LAMA sample because name ($sample_name) does not fit regex $name_regex [$sample_print_info]");
            next;
        }

        my $original_submission = $sample_to_store{submission};
        my $isolation_type = $sample_to_store{isolation_type};
        my $prep_type = $sample_to_store{prep_type} || NACHAR;
        my $analysis_type = $sample_to_store{isolation_type};
        my $final_target_yield = NACHAR;

        if (not defined $isolation_type) {
            sayWarn("NOTIFY: LAMA sample has no isolation type defined for $isolate_barcode [$sample_print_info]");
            next;
        }
        elsif ($isolation_type eq 'TUMOR_FFPE_DNA_ISOLATE' or $isolation_type eq 'Tumor FFPE' or $prep_type eq "PANEL") {
            $analysis_type = 'Targeted_Tumor_Only'; # DNA from tumor tissue
            # sanity check that we are indeed dealing with TO (tumor only)
            my $expected_cohort = 'TARGTO';
            if (not exists $sample_to_store{cohort}){
                sayWarn(sprintf "NOTIFY: Found FFPE sample but NO COHORT present for $isolate_barcode (pls fix in LAMA) [$sample_print_info]");
                next;
            }elsif ($sample_to_store{cohort} ne $expected_cohort){
                sayWarn(sprintf "NOTIFY: Found FFPE sample but cohort [%s] is not [%s] pls check LAMA [%s]",
                    $sample_to_store{cohort}, $expected_cohort, $sample_print_info);
                next;
            }
            $final_target_yield = 50;
        }
        elsif ($isolation_type eq 'TUMOR_DNA_ISOLATE' or $isolation_type eq 'Tissue') {
            $analysis_type = 'Somatic_T'; # DNA from tumor tissue
            $final_target_yield = 300;
        }
        elsif ($isolation_type eq 'TUMOR_RNA_ISOLATE' or $isolation_type eq 'RNA') {
            $analysis_type = 'RNAanalysis'; # RNA from tumor tissue
            $final_target_yield = 15;
        }
        elsif ($isolation_type eq 'REFERENCE_DNA_ISOLATE' or $isolation_type eq 'Blood') {
            $analysis_type = 'Somatic_R'; # DNA from blood
            $dna_reference_samples_by_name{ $sample_name } = \%sample_to_store;
            $final_target_yield = 100;
        }
        elsif ($isolation_type eq 'Plasma') {
            $analysis_type = 'PlasmaAnalysis'; # Plasma from blood
        }
        else {
            sayWarn("NOTIFY: unknown isolation type '$isolation_type' for $isolate_barcode (pls fix in LAMA) [$sample_print_info]");
            next
        }

        # The submission is no longer set in LAMAv2 but is required for registration so added
        if (not defined $original_submission or $original_submission eq '') {
            $original_submission = "HMFreg" . $sample_to_store{label};
        }

        $sample_to_store{original_submission} = $original_submission;
        $sample_to_store{submission} = $original_submission;
        $sample_to_store{project_name} = $original_submission;

        if (exists $centers_dict->{ $center } and $sample_name !~ /^CORE01/) {
            # All meant-for-database should by from a known center
            my $centername = $centers_dict->{ $center };
            $sample_to_store{entity} = join("_", $study, $centername);
        }
        else {
#            sayInfo(sprintf "Unknown center [%s] or CORE01 so configured [%s] as entity [%s]", $center, ONCOACT_ENTITY, $sample_print_info);
            $sample_to_store{entity} = ONCOACT_ENTITY;
        }

        # Add the missing fields and store final
        $sample_to_store{analysis_type} = $analysis_type;
        $sample_to_store{original_submission} = $original_submission;
        $sample_to_store{yield} = $final_target_yield;

        # Fix various formats of date fields
        fixDateFields(\%sample_to_store);

        # Fix/translate various field contents
        fixFieldContents(\%sample_to_store, $name_dict->{lama_content_translations_by_field_name});
        $sample_to_store{patient} =~ s/\-//g;
        $sample_to_store{cohort} =~ s/\-//g;

        # Add non-existing fields that might be required for downstream tools to work
        my @fields_that_must_be_present = qw(arrival_date biopsy_site lab_sop_versions ptum report_germline_level);
        foreach my $field (@fields_that_must_be_present){
            if ( not exists $sample_to_store{$field} ){
                $sample_to_store{$field} = "";
            }
        }

        # Define data source and store the final result
        $sample_to_store{"".DATA_SOURCE_FIELD} = "LAMA";
        storeRecordByKey(\%sample_to_store, $isolate_barcode, \%store, "final storing of $isolate_barcode", 1);
    }

    # Need another loop over all tumor samples to complete ref_sample_id for older samples
    while ( my($barcode, $sample) = each %store ){

        # We can skip non tumor samples and tumor samples where ref_sample_id info is already present
        next unless $sample->{analysis_type} eq 'Somatic_T';
        next if (defined $sample->{'ref_sample_id'} and $sample->{'ref_sample_id'} ne "");

        # Otherwise try to complete info by searching for R sample
        my $patient_string = $sample->{ 'patient' };
        $patient_string =~ s/\-//g; # string this point still with dashes
        my $ref_sample_name = $patient_string . 'R';
        if ( exists $dna_reference_samples_by_name{ $ref_sample_name } ){
            my $ref_sample_id = $dna_reference_samples_by_name{ $ref_sample_name }{ 'sample_id' };
            $sample->{'ref_sample_id'} = $ref_sample_id;
        }
    }
    return \%store;
}

sub contractToCohort{
    my ($contract_contract_display_name) = @_;
    # Note: Granularity within CORE does not exist in LAMAv2
    # So CORE COREDB COREDB08 COREDB11 CORELR02 CORELR11 CORERI02 CORESC11 are all configured to COREDB
    my %translation = (
        'CORE' => 'COREDB',
        'Last Resort' => 'COREDB',
        'Radium Insight' => 'COREDB',
        'CPCT' => 'CPCT',
        'CPCT Blinc' => 'CPCTBLINC',
        'CPCT Pegasus' => 'CPCTpancreas',
        'DRUP 3rd stage' => 'DRUPstage3',
        'Panel' => 'TARGTO' # Targeted Tumor Only
    );
    if (exists $translation{$contract_contract_display_name}){
        return $translation{$contract_contract_display_name};
    }
    else {
        return $contract_contract_display_name;
    }
}

sub fixFieldContents{
    my ($record, $translation_dict) = @_;
    while ( my($key, $translations) = each %$translation_dict ){
        if ( exists $record->{$key} ){
            my $val = $record->{$key};
            if ( exists $translations->{$val} ){
                $record->{$key} = $translations->{$val};
            }
        }
    }
}

sub storeRecordByKey{
    my ($record, $key, $store, $info_tag, $do_not_warn_if_exists) = @_;
    if (exists $store->{$key} and not $do_not_warn_if_exists){
        sayWarn("Store key already exists in store and will be overwritten ($key for $info_tag)");
    }
    my %copy_of_record = %$record;
    $store->{$key} = \%copy_of_record;
}

sub addRecordFieldsToTargetRecord{
    my ($record, $target_record) = @_;
    while (my ($key, $val) = each %$record){
        $target_record->{$key} = $val;
    }
}

sub copyFieldsFromObject{
    my ($object, $info_tag, $fieldsTranslationTable, $store) = @_;

    if (scalar keys %$fieldsTranslationTable == 0){
        die "[ERROR] No fields to translate ($info_tag)\n";
    }

    while (my ($src_key, $tgt_key) = each %$fieldsTranslationTable){
        if (exists $object->{$src_key}){
            if (ref($object->{$src_key}) eq 'ARRAY'){
                $store->{$tgt_key} = join(",", @{$object->{$src_key}});
            }else{
                $store->{$tgt_key} = $object->{$src_key};
            }
        }
        elsif(exists $WARN_IF_ABSENT_IN_LAMA_FIELDS{$src_key}){
            sayWarn("No '$src_key' field in object ($info_tag)");
        }
    }
}

sub epochToDate{
    my ($epoch) = @_; # epoch time in milliseconds
    my $registrationDate = strftime "%Y-%m-%d", localtime $epoch/1000;
    return $registrationDate;
}

sub parseLamaPatients {
    my ($inputJsonFile) = @_;
    my %store = ();
    my $objects = parseJson($inputJsonFile);

    foreach my $patient (@$objects) {
        foreach my $sample (@{$patient->{tumorSamples}}) {
            processSampleOfLamaPatient($patient, $sample, 'tumor', \%store);
        }
        foreach my $sample (@{$patient->{referenceSamples}}) {
            processSampleOfLamaPatient($patient, $sample, 'reference', \%store);
        }
    }
    return \%store;
}

sub processSampleOfLamaPatient {
    my ($patient, $sample, $sample_origin, $store) = @_;
    my $sample_field_translations;
    my $patient_id = $patient->{patientId};
    my @sampleBarcodes;

    my %info = ();
    my $info_tag = "patients->" . $patient_id;

    if( $sample_origin eq 'tumor' ){
        $sample_field_translations = $name_dict->{lama_patient_tumor_sample_dict};
        @sampleBarcodes = @{$sample->{sampleBarcodes}};
        copyFieldsFromObject($sample->{'tumorType'}, $info_tag, $name_dict->{lama_patient_tumor_sample_tumor_type_dict}, \%info);
        copyFieldsFromObject($sample->{'biopsy'}, $info_tag, $name_dict->{lama_patient_tumor_sample_biopsy_dict}, \%info);
        $info{'ptum'} = constructPrimaryTumorTypeField($sample);
        $info{'biopsy_site'} = constructBiopsyField($sample);
    }
    elsif( $sample_origin eq 'reference' ){
        $sample_field_translations = $name_dict->{lama_patient_reference_sample_dict};
        @sampleBarcodes = ($sample->{sampleBarcode});
    }
    else{
        my $sample_barcode = $sample->{sampleBarcode};
        die "[ERROR] Unknown sample origin provided to processSampleOfPatient ($sample_origin for $sample_barcode)\n";
    }

    copyFieldsFromObject($sample, $info_tag, $sample_field_translations, \%info);
    copyFieldsFromObject($patient, $info_tag, $name_dict->{lama_patient_dict}, \%info);

    foreach my $sampleBarcode (@sampleBarcodes) {
        storeRecordByKey(\%info, $sampleBarcode, $store, "patient_samples");
    }
}

sub constructPrimaryTumorTypeField{
    my ($sample) = @_;
    my $top_level_key = "tumorType";
    my @sub_level_keys = qw(location type extra);
    my @values = ();
    foreach my $key (@sub_level_keys){
        if (exists $sample->{$top_level_key}{$key} and $sample->{$top_level_key}{$key} ne ""){
            push @values, $sample->{$top_level_key}{$key};
        }
    }
    return join(" | ", @values);
}

sub constructBiopsyField{
    my ($sample) = @_;
    my $top_level_key = "biopsy";
    my @sub_level_keys = qw(location subLocation lateralisation);
    my @values = ();
    foreach my $key (@sub_level_keys){
        if (exists $sample->{$top_level_key}{$key}){
            if ($sample->{$top_level_key}{$key} ne "" and $sample->{$top_level_key}{$key} ne "Other/unknown"){
                push @values, $sample->{$top_level_key}{$key};
            }
        }
    }
    return join(" | ", @values);
}

sub parseLamaLibraryPreps{
    my ($inputJsonFile) = @_;
    my %store = ();
    my $objects = parseJson($inputJsonFile);

    foreach my $experiment (@$objects) {

        foreach my $object (@{$experiment->{libraries}}) {
            my $store_key = $object->{'isolationBarcode'};
            my $status = $object->{'status'};

            # Only store prep info when OK
            next if $status =~ m/failed/i;

            my %info = ();
            my $info_tag = "libraryprep->$store_key";
            copyFieldsFromObject($object, $info_tag, $name_dict->{lama_libraryprep_library_dict}, \%info);
            storeRecordByKey(\%info, $store_key, \%store, "libraryprep", 1);
        }
    }
    return \%store;
}

sub parseLamaIsolation{
    my ($inputJsonFile) = @_;
    my %store = ();
    my $objects = parseJson($inputJsonFile);

    foreach my $experiment (@$objects) {

        # When isolation experiment is Processing there are no frBarcodes yet so skip all isolates
        next if $experiment->{status} eq "PROCESSING";

        foreach my $isolate (@{$experiment->{isolates}}) {

            # Once a isolation experiment is no longer processing there should be a an isolation barcode
            unless ( exists $isolate->{isolationBarcode} ) {
                print Dumper $isolate;
                die "[ERROR] No frBarcode present for above isolate object";
            }

            my $barcode = $isolate->{isolationBarcode};
            my $new_status = $isolate->{status};

            if ( exists $store{$barcode} ) {
                my $old_status = $store{$barcode}{'isolation_status'};
                my $old_is_finished = $old_status eq 'FINISHED';
                my $new_is_finished = $new_status eq 'FINISHED';
                if ( $old_is_finished and $new_is_finished ){
                    sayWarn("SKIPPING isolate: encountered duplicate Finished isolate for $barcode (pls fix in LAMA)");
                    next;
                }
                elsif ( $old_is_finished and not $new_is_finished ){
                    next;
                }
                else{
                    # In this case we simply overwrite a non-Finished status with a Finished one
                }
            }
            my %info = ();
            my $info_tag = "isolation->$barcode";
            copyFieldsFromObject($isolate, $info_tag, $name_dict->{lama_isolation_isolate_dict}, \%info);
            storeRecordByKey(\%info, $barcode, \%store, "isolation", 1);
        }
    }
    return \%store;
}

sub parseLamaSampleStatus{
    my ($inputJsonFile) = @_;
    my %store = ();
    my $objects = parseJson($inputJsonFile);

    foreach my $object (@$objects){

        if (! (exists $object->{sampleBarcode} && exists $object->{sampleId}) ){
            print Dumper $object;
            die "[ERROR] No sampleBarcode/sampleId in object!";
        }

        my $sampleBarcode = $object->{sampleBarcode};
        my $sampleId = $object->{sampleId};
        my $infoTag = "Status:$sampleBarcode/$sampleId";

        # No DNA frBarcode means sample has not been isolated so skip
        if (not defined $sampleBarcode or $sampleBarcode eq ""){
            sayWarn("SKIPPING: No barcode for [$infoTag]");
            next;
        }

        # Collect all info into one object
        my %status = ();
        copyFieldsFromObject($object, $infoTag, $name_dict->{lama_status_dict}, \%status);
        copyFieldsFromObject($object->{labWorkflow}, $infoTag, $name_dict->{lama_status_labworkflow_dict}, \%status);

        # Store
        foreach my $isolationBarcode (@{$object->{isolationBarcodes}}) {
            if (not defined $isolationBarcode or $isolationBarcode eq ""){
                sayWarn("SKIPPING isolation because undefined barcode [$infoTag]");
                next;
            }
            $status{'sample_id'} = $isolationBarcode;
            storeRecordByKey(\%status, $isolationBarcode, \%store, $infoTag);
        }
    }
    return \%store;
}

sub parseLamaContracts{
    my ($inputJsonFile) = @_;
    my %store = ();
    my $objects = parseJson($inputJsonFile);

    foreach my $object (@$objects) {
        my $store_key = $object->{'code'};

        my %info = ();
        my $info_tag = "contracts->$store_key";
        copyFieldsFromObject($object, $info_tag, $name_dict->{lama_contracts_dict}, \%info);
        copyFieldsFromObject($object->{reportSettings}, $info_tag, $name_dict->{lama_contracts_report_dict}, \%info);
        storeRecordByKey(\%info, $store_key, \%store, "contracts", 0);
    }
    return \%store;
}

sub parseTsvCsv{
    my ($objects, $fields, $store_field_name, $should_be_unique, $file, $sep) = @_;
    my $csv = Text::CSV->new({ binary => 1, auto_diag => 1, sep_char => $sep });
    my %store = %$objects;
    my $line_count;
    open(COUNT_LINES, "< $file") or die "can't open $file: $!";
    while ( <COUNT_LINES> ){
        next if /^\s*$/;
        $line_count++
    }
    close(COUNT_LINES);

    sayInfo(sprintf "  Parsing CSV/TSV [%s,nonEmptyLineCount=%s]", basename($file), $line_count);
    open IN, "<", $file or die "[ERROR] Unable to open file ($file): $!\n";
    my $header_line = <IN>; chomp($header_line);
    die "[ERROR] Cannot parse line ($header_line)\n" unless $csv->parse($header_line);
    my @header_fields = $csv->fields();
    my %fields_map = map { $_ => 1 } @header_fields;

    # Checking header content
    my $header_misses_field = 0;
    foreach my $field (keys %$fields) {
        if ( not exists $fields_map{ $field } ){
            sayWarn("Missing header field ($field) in file ($file)");
            $header_misses_field = 1;
        }
    }
    if ( $header_misses_field ){
        print Dumper \%fields_map and die "[ERROR] Header incomplete ($file)\n";
    }

    # Header OK: continue reading in all data lines
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
            sayWarn("SKIPPING object (name: $name) from $file for reason: $reason_not_to_store") and next;
        }

        # Checks OK: fix some fields and store object
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

sub parseJson{
    my ($file) = @_;
    my $size = `du -shL $file | cut -f1`;
    chomp($size);
    sayInfo(sprintf "  Parsing JSON [%s,size=%s]", basename($file), $size);
    return(decode_json(read_file($file)));
}

sub printLimsToJson{
    my ($samples, $submissions, $contact_groups, $lims_file) = @_;
    my $samp_count = scalar keys %$samples;
    my $subm_count = scalar keys %$submissions;
    my $cont_count = scalar keys %$contact_groups;

    my %lims = ( 'samples' => $samples, 'submissions' => $submissions, 'contact_groups' => $contact_groups );
    my $coder = JSON::XS->new->utf8->canonical;
    my $lims_txt = $coder->encode(\%lims);

    sayInfo("  Writing output to $lims_file ($cont_count contact groups, $subm_count submissions and $samp_count samples)");
    open my $lims_json_fh, '>', $lims_file or die "[ERROR] Unable to open output file ($lims_file): $!\n";
    print $lims_json_fh $lims_txt;
    close $lims_json_fh;
}

sub addLabSopString{
    my ($samples, $inprocess) = @_;
    sayInfo("  Adding SOP string information to samples");
    my %store = %$samples;
    my $sop_field_name = 'lab_sop_versions';
    foreach my $id ( keys %store ){
        if ( exists $inprocess->{ $id } ){
            # format: PREP(\d+)V(\d+)-QC(\d+)V(\d+)-SEQ(\d+)V(\d+)
            $store{ $id }{ $sop_field_name } = $inprocess->{ $id }{ $sop_field_name };
        }
        elsif ( defined $samples->{ $id }{ $sop_field_name } ){
            # keep whatever is present
        }
        else{
            # fallback to NA default
            $store{ $id }{ $sop_field_name } = NACHAR;
        }
    }
    return \%store;
}

sub parseDictFile{
    my ($file, $fileType) = @_;
    sayInfo(sprintf "  Parsing DICTIONARY [%s]", basename($file));
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
    my @required_fields = qw(
        requester_report_contact_name
        requester_report_contact_email
        report_contact_name
        report_contact_email
        data_contact_name
        data_contact_email
    );
    sayInfo("  Adding contact group information to submissions");
    my %store = %{$submissions};
    foreach my $submission_id (sort keys %store){
        my $submission = $store{$submission_id};

        if( exists $submission->{ 'report_contact_email' } ){
            # Submission is from the time when contact info was entered in shipments tab
            foreach my $field ( @required_fields ){
                unless ( exists $submission->{$field} ){
                    $submission->{$field} = NACHAR;
                }
            }
            next;
        }
        elsif( defined $submission->{ 'group_id' } ){
            # Retrieve contact info from the contact group table
            my $group_id = $submission->{ 'group_id' };
            next if $group_id eq 'na';
            if ( exists $contact_groups->{ $group_id } ){
                my $group = $contact_groups->{ $group_id };
                foreach my $field ( @required_fields ){
                    $submission->{$field} = $group->{ $field };
                }
            }
            else{
                sayWarn("Submission \"$submission_id\" has group ID \"$group_id\" but not found in contact groups (check Contact tab of FOR-001)");
                next;
            }
        }
        else{
            sayWarn("No Group ID defined for submission $submission_id");
            next;
        }
    }
    return \%store;
}

sub checkDrupStage3Info{
    my ($submissions, $samples) = @_;
    my $counter = 0;
    my $cohort = 'DRUPstage3';

    sayInfo(sprintf "  Checking submissions of %s cohort", $cohort);
    while (my ($id,$obj) = each %$submissions){
        my $project_type = $obj->{'project_type'};
        my $project_name = $obj->{'project_name'};

        next unless ($project_type eq 'Cohort' and $project_name =~ /^DRUP/);

        $counter++;
        my $patient_id = $project_name;
        $patient_id =~ tr/-//d;
        my $sample_name = $patient_id . "T";
        my $sample_count = 0;
        while (my($barcode,$sample) = each %$samples){
            next unless $sample->{'sample_name'} eq $sample_name;
            next unless $sample->{'original_submission'} =~ m/^$id$/i;
            next unless $sample->{'analysis_type'} eq 'Somatic_T';
            my $sample_cohort = $sample->{'cohort'};
            if ( uc($sample_cohort) ne uc($cohort) ){
                sayWarn(sprintf "Found sample [%s,%s] of submission [%s] with incorrect cohort [%s] instead of expected cohort [%s]",
                  $sample_name, $barcode, $id, $sample_cohort, $cohort);
            }
            $sample_count++;
        }
        sayWarn(sprintf "Found no samples for %s submission %s!", $cohort, $id) if $sample_count < 1;
    }
    sayInfo(sprintf "Summary: %d submissions encountered and checked for %s", $counter, $cohort);
}

sub checkContactInfo{
    my ($contact_groups) = @_;
    my @name_fields = qw(client_contact_name requester_report_contact_name report_contact_name data_contact_name);
    my @mail_fields = qw(requester_report_contact_email report_contact_email data_contact_email);
    sayInfo("  Checking contact group information for completeness");
    foreach my $id (sort keys %$contact_groups){
        my $info = $contact_groups->{$id};

        # These fields should at the very least have content
        foreach my $field (@name_fields, @mail_fields){
            if ( $info->{$field} eq "" ){
                sayWarn("No content in field [$field] for contact group ID [$id] (see FOR-001 Contacts tab)");
            }
        }
        # These fields should contain only (valid) email addresses
        foreach my $field (@mail_fields){
            my @addresses = split( ";", $info->{$field});
            foreach my $address (@addresses){
                if( $address eq NACHAR ){
                    next;
                }
                elsif( not Email::Valid->address($address) ){
                    sayWarn("Invalid email address [$address] in field [$field] for contact group ID [$id] (see FOR-001 Contacts tab)");
                }
            }
        }
    }
}

sub addExcelSamples{

    my ($lims, $objects, $shipments) = @_;
    my %store = %{$lims};

    # Open file and check header before reading data lines
    while ( my($row_key, $row_info) = each %$objects ){

        my $sample_name = $row_info->{ 'sample_name' } or die "[ERROR] No sample_name in row_info";
        next if isSkipValue( $sample_name );

        my $sample_id = $row_info->{ 'sample_id' } or sayWarn("SKIPPING sample ($sample_name): No sample_id found") and next;
        my $submission = $row_info->{ 'submission' } or sayWarn("SKIPPING sample ($sample_name, $sample_id): No submission found") and next;
        my $analysis_type = $row_info->{ 'analysis_type' } or sayWarn("SKIPPING sample ($sample_name, $sample_id): No analysis_type found") and next;

        $row_info->{ 'label' } = 'RESEARCH';
        $row_info->{ 'patient' } = $row_info->{ 'sample_name' };
        $row_info->{ 'reporting_id' } = $row_info->{ 'sample_name' };
        $row_info->{ 'entity' } = $row_info->{ 'submission' };

        # Check data analysis type and set accordingly
        if ( $analysis_type =~ /^(Somatic_R|Somatic_T|SingleAnalysis|FASTQ|BCL|LabOnly)$/ ){
            # Already final status so no further action
        }
        elsif ( $sample_name =~ /^(CORE\d{2}\d{6})(T|R){1}/ms ){
            my ($patient, $tum_or_ref) = ($1, $2);
            $row_info->{ 'label' }         = "CORE";
            $row_info->{ 'patient' }       = $patient;
            $row_info->{ 'entity' }        = $submission;
            $row_info->{ 'analysis_type' } = $tum_or_ref eq 'T' ? 'Somatic_T' : 'Somatic_R';
        }
        elsif ( $analysis_type eq 'SomaticAnalysis' or $analysis_type eq 'SomaticsBFX' ){
            # SomaticsBFX is the old term, SomaticAnalysis the new
            my $partner = $row_info->{ 'ref_sample_id' };
            if ( $partner ne '' and $partner ne NACHAR ){
                $row_info->{ 'analysis_type' } = 'Somatic_T';
            }
            else{
                $row_info->{ 'analysis_type' } = 'Somatic_R';
            }

            # Hardcode Somatic_T samples to not use existing ref data for FOR-001 samples
            $row_info->{ 'other_ref' } = "";

            # Hardcode old Somatic_T samples to not run in shallow mode (config was added only in FOR-001 v5.10)
            if ( not exists $row_info->{ 'shallowseq' }) {
                $row_info->{ 'shallowseq' } = JSON::XS::false;
            }

        }
        elsif ( $analysis_type eq 'GermlineBFX' or $analysis_type eq 'Germline' ){
            # GermlineBFX is the old term, SingleAnalysis the new
            $analysis_type = 'SingleAnalysis';
            $row_info->{ 'analysis_type' } = $analysis_type;
        }
        elsif ( $analysis_type eq 'NoBFX' or $analysis_type eq 'NoAnalysis' or $analysis_type eq '' or $analysis_type eq 'NA' ){
            # NoBFX is the old term, FASTQ the new
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
        elsif ( $analysis_type eq 'SOMATIC PANEL' ){
            $analysis_type = 'SomaticPanel';
            $row_info->{ 'analysis_type' } = $analysis_type;
        }
        else {
            sayWarn("SKIPPING sample ($sample_name): has unknown analysis type ($analysis_type)");
            next;
        }

        # Add submission info and parse KG
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
            elsif( $sub->{ 'group_id' } and $sub->{ 'group_id' } eq 'HMF-INNOVATION' ){
                $row_info->{ 'entity' } = 'INNOVATION';
                $row_info->{ 'label' } = 'INNOVATION';
            }
            # Assumes that all samples of submission need same analysis
            $sub->{ 'analysis_type' } = $analysis_type;
        }

        my $unique = $row_info->{ 'sample_id' };
        next if isSkipValue( $unique );

        # Checks before storing
        my $regex = '^[0-9a-zA-Z\-]*$';
        sayWarn("SKIPPING sample ($sample_name): sample_name contains unacceptable chars") and next if $sample_name !~ /$regex/;
        sayWarn("SKIPPING sample ($sample_name): sample_id ($sample_id) contains unacceptable chars") and next if $sample_id !~ /$regex/;
        sayWarn("SKIPPING sample ($sample_name): no submission defined for sample") and next unless $row_info->{ 'submission' };
        sayWarn("SKIPPING sample ($sample_name): no analysis type defined for sample") and next unless $row_info->{ 'analysis_type' };
        sayWarn("SKIPPING sample ($sample_name): no project name defined for sample") and next unless $row_info->{ 'project_name' };

        # Store at unique id
        my $reason_not_to_store = checkKeyToStore( \%store, $unique );
        if ( $reason_not_to_store ){
            sayWarn("SKIPPING sample with name \"$sample_name\" for reason: $reason_not_to_store") and next;
        }
        $row_info->{ "".DATA_SOURCE_FIELD } = "EXCEL";
        $store{ $unique } = $row_info;

    }

    return \%store;
}

sub fixIntegerFields{
    my ($obj) = @_;
    foreach my $key ( INTEGER_FIELDS ){
        # Make sure all integer values are stored as such for json export
        if ( exists $obj->{$key} and $obj->{$key} =~ /^\d+$/ ){
            $obj->{$key} = $obj->{$key} + 0;
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
            my $name = $obj->{ 'sample_name' } || "SampleNameUnknown";
            sayWarn("Unexpected value ($value) in boolean field ($key) for sample ($name)");
        }
    }
}

sub fixAddedDateFields{
    my ($sample_objects) = @_;
    while( my($key, $obj) = each %$sample_objects){
        fixDateFields( $obj );
    }
}

sub fixDateFields{
    my ($obj) = @_;
    my @date_fields = qw( arrival_date sampling_date report_date isolation_date libraryprep_date snpgenotype_date );

    foreach my $date_field ( @date_fields ){

        next unless defined $obj->{ $date_field };
        my $old_date = $obj->{ $date_field };
        my $new_date = $old_date;
        my $identifier = $obj->{ 'sample_name' };

        # Date is not always filled in so skip NA fields
        next if isSkipValue( $old_date );

        # Convert all date strings to same format yyyy-mm-dd (eg 2017-01-31)
        if( $old_date eq '1' ) {
            $new_date = NACHAR;
        }
        elsif( $old_date =~ /^\d{13}$/ ) {
            # eg 1516575600000
            $new_date = epochToDate($old_date)
        }
        elsif( $old_date =~ /^\w+ (\w{3}) (\d{2}) \d+:\d+:\d+ \w+ (\d{4})$/ ){
            # eg Tue Apr 23 00:00:00 CEST 2019
            my $month_name = $1;
            my $day = $2;
            my $year = $3;
            my $month = getMonthIndexByName( $month_name );
            $new_date = join( "-", $year, $month, $day );
        }
        elsif ( $old_date =~ /^(\d{2})(\d{2})(\d{2})$/ ){
            # Format unclear so need for checks
            sayWarn("Date \"$old_date\" in \"$date_field\" has unexpected year ($identifier): please check") if ($1 < 8) or ($1 > 20);
            sayWarn("Date \"$old_date\" in \"$date_field\" has impossible month ($identifier): please fix") if $2 > 12;
            $new_date = join( "-", "20" . $1, $2, $3 );
        }
        elsif ( $old_date =~ /^(\d{2})-(\d{2})-(\d{4})$/ ){
            # case dd-mm-yyyy
            sayWarn("Date \"$old_date\" in \"$date_field\" has impossible month ($identifier): please fix") if $2 > 12;
            $new_date = join( "-", $3, $2, $1 );
        }
        elsif ( $old_date =~ /^(\d{4})-(\d{2})-(\d{2})$/ ){
            # case yyyy-mm-dd already ok
            sayWarn("Date \"$old_date\" in \"$date_field\" has impossible month ($identifier): please fix") if $2 > 12;
        }
        elsif ( exists $old_date->{'$numberLong'}){
            # Older versions of mongo-export use canonical mode and return a hash with numberLong key
            $new_date = epochToDate($old_date->{'$numberLong'})
        }
        else{
            sayWarn("Date string \"$old_date\" in field \"$date_field\" has unknown format for sample ($identifier): kept string as-is but please fix");
        }

        # Store new format using reference to original location
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
        sayWarn("Unknown Month name ($month_name): kept as-is but please fix");
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

    sayInfo("Loading excel file $excel sheet '$sheet'");
    my $workbook = Spreadsheet::XLSX->new( $excel ) or die "[ERROR] Unable to load excel file $excel: $!\n";
    my $sheet_obj = $workbook->worksheet( $sheet ) or die "[ERROR] Unable to read sheet \"$sheet\" from file $excel: $!\n";

    my @header = ();
    my $max_row = $sheet_obj->{'MaxRow'};
    my $max_col = $sheet_obj->{'MaxCol'};

    # Check if header exist where it should be
    my $first_val = EMPTY;
    my $first_cel = $sheet_obj->get_cell( $h_row, $h_col );
    $first_val = $first_cel->unformatted() if defined $first_cel;
    die "[ERROR] Header value ($h_val) cannot be found at set location ($excel)\n" unless $first_val eq $h_val;

    # Now read header values for later storage
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
    my @to_skip = ( NACHAR, EMPTY, '', 'na', 'naR', 'naT', 'invalid', 'failed', 'nvt', 'no', 'x', '#N/A' );
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
        sayWarn("Value $key is not defined in:");
        print Dumper $hash;
    }
}

sub sayInfo{
    my ($msg) = @_;
    say "[INFO] " . (strftime "%y%m%d %H:%M:%S", localtime) . " - " . $msg;
}

sub sayWarn{
    my ($msg) = @_;
    warn "[WARN] " . (strftime "%y%m%d %H:%M:%S", localtime) . " - " . $msg . "\n";
}

sub getFieldNameTranslations{
    # Columns contact sheet in current FOR-001
    my %CONT_DICT = (
        "Group_ID"                      => 'group_id',
        "Client_contact_name"           => 'client_contact_name',
        "Client_contact_email"          => 'client_contact_email',
        "On_behalf_of_client_name"      => 'on_behalf_of_client_contact_name',
        "On_behalf_of_client_email"     => 'on_behalf_of_client_contact_email',
        "Report_contact_name"           => 'report_contact_name',
        "Report_contact_email"          => 'report_contact_email',
        "Requester_report_contact_name" => 'requester_report_contact_name',
        "Requester_report_contact_email" => 'requester_report_contact_email',
        "Data_contact_name"      => 'data_contact_name',
        "Data_contact_email"     => 'data_contact_email',
        "Lab_contact_name"       => 'lab_contact_name',
        "Lab_contact_email"      => 'lab_contact_email',
    );

    # Columns shipments sheet in 2018 rest lims (FOR-001)
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

    # Columns shipments sheet in 2019 rest lims (FOR-001)
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

    # Columns shipments sheet in rest lims (FOR-001)
    my %SUBM_DICT_2020 = (
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

    # Columns shipments sheet in rest lims (FOR-001)
    my %SUBM_DICT_2021 = (
        "Arrival_date"      => 'arrival_date',
        "Project_name"      => 'project_name',
        "HMF_reg"           => 'submission',
        "Requested_product" => 'request',
        "Product_category"  => 'project_type',
        "Sample_count"      => 'sample_count',
        "Lab_is_finished"   => 'has_lab_finished',
        "Group_ID"          => 'group_id',
        "Total_yield_required" => 'total_yield_required',
        "Remarks"           => 'remarks',
    );

    my %SUBM_DICT_2022 = %SUBM_DICT_2021;
    my %SUBM_DICT_2023 = %SUBM_DICT_2022;
    my %SUBM_DICT = %SUBM_DICT_2023;

    # Columns samples sheet in 2018 FOR-001
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

    # Columns samples sheet in 2019 FOR-001
    my %SAMP_DICT_2019 = (
        "Sample_ID"           => 'sample_id',
        "Sample_name"         => 'sample_name',
        "DNA_conc"            => 'conc',
        "Yield"               => 'yield',
        "Q30"                 => 'q30',
        "Analysis_type"       => 'analysis_type',
        "Partner_sample"      => 'ref_sample_id',
        "HMF_reg"             => 'submission',
        "SNP_required"        => 'is_snp_required',
        "SNP_exp"             => 'snp_experiment_id',
        "Requested_product"   => 'request',
        "State"               => 'lab_status', # lab status
        "Priority"            => 'priority',
        "Arrival_date"        => 'arrival_date',
        "Remarks"             => 'remarks',
    );

    my %SAMP_DICT_2020 = %SAMP_DICT_2019;
    my %SAMP_DICT_2021 = %SAMP_DICT_2020;
    $SAMP_DICT_2021{"ShallowSeq_required"} = 'shallowseq';
    my %SAMP_DICT_2022 = %SAMP_DICT_2021;
    my %SAMP_DICT_2023 = %SAMP_DICT_2022;
    my %SAMP_DICT = %SAMP_DICT_2023;

    # Columns In Process sheet (HMF-FOR-002)
    my %PROC_DICT_2022 = (
        'Sample_ID'         => 'sample_id', # eg FR12345678
        'Sample_name'       => 'sample_name', # eg CPCT1234567R
        'Diluted_library'   => 'library_id', # eg FR12345678 (THIS WAS "barcode_3nm")
        'Sop_tracking_code' => 'lab_sop_versions',
    );
    my %PROC_DICT_2023 = %PROC_DICT_2022;
    my %PROC_DICT = %PROC_DICT_2023;

    my %lama_patient_dict = (
        '_id' => 'lama_patient_object_id',
        'patientId' => 'patient',
        'hospitalPatientId' => 'hospital_patient_id',
        'reportingId' => 'reporting_id',
    );

    my %lama_patient_tumor_sample_dict = (
        'legacySampleId'        => 'legacy_sample_name',
        'refIsolationBarcode'   => 'ref_sample_id', # was refFrBarcode pre-lama-v2
        'pathologyNumber'       => 'hospital_pa_sample_id',
        'hospitalSampleLabel'   => 'hospital_sample_label',
        'patientGermlineChoice' => 'report_germline_level',
        'primaryTumorType'      => 'ptum',
        'biopsySite'            => 'lama_tumor_sample_biopsy_site',
        'sopVersion'            => 'sop_version',
        'samplingDate'          => 'sampling_date', # was collectionDate pre-lama-v2
        'isCUP'                 => 'is_cup',
        'arrivalHmf'            => 'arrival_date',
        'submissionNr'          => 'submission',
        'contractCode'          => 'contract_code'
    );

    my %lama_patient_tumor_sample_tumor_type_dict = (
        'location' => 'tumor_location',
        'type'  => 'tumor_type',
        'extra' => 'tumor_extra',
        'doids' => 'tumor_doids',
    );

    my %lama_patient_tumor_sample_biopsy_dict = (
        'biopsyLocation' => 'biopsy_location',
        'biopsySubLocation'  => 'biopsy_sub_location',
        'lateralisation' => 'biopsy_lateralisation',
    );

    my %lama_patient_reference_sample_dict = (
        'legacySampleId'  => 'legacy_sample_name',
        'sopVersion'      => 'sop_version',
        'samplingDate'    => 'sampling_date', # was collectionDate pre-lama-v2
        'arrivalHmf'      => 'arrival_date',
        'sampleBarcode'   => 'sample_barcode',
        'submissionNr'    => 'submission',
        'originalBarcode' => 'original_barcode',
        'contractCode'    => 'contract_code'
    );

    my %lama_isolation_isolate_dict = (
        'sampleBarcode' => 'sample_barcode',
        'isolationBarcode' => 'isolation_barcode',
        'experimentNr'  => 'isolation_id', # was isolationNr pre-lama-v2
        'coupeBarcode'  => 'coupes_barcode', # confirmed the same in lama-v2
        'status'        => 'isolation_status', # confirmed the same in lama-v2
        'resultType'    => 'isolation_type', # content change from Blood|RNA|Tissue to
        'concentration' => 'conc'
    );

    my %lama_libraryprep_library_dict = (
        'prepType'     => 'prep_type',
        'experimentNr' => 'prep_id', # was prepNr pre-lama-v2
        'status'       => 'prep_status',
        'isolationBarcode' => 'isolation_barcode',
        'undilutedBarcode' => 'undiluted_barcode',
        'dilutedBarcode' => 'diluted_barcode'
    );

    my %lama_status_dict = (
        '_id'                  => 'lama_status_object_id',
        'libraryPrepStatus'    => 'lab_status',
        'registrationDate'     => 'registration_date',
        'sampleId'             => 'sample_name',
        'sampleBarcode'        => 'sample_barcode',
        'type'                 => 'status_type'
    );

    my %lama_status_labworkflow_dict = (
        'doShallow'    => 'shallowseq',
        'isTumorPanel' => 'is_tumor_panel',
        'isTumorOnly'  => 'is_tumor_only'
    );

    my %lama_contracts_dict = (
        '_id' => 'lama_contract_object_id',
        'type' => 'contract_type',
        'code' => 'contract_code',
        'displayName' => 'contract_display_name',
        'sampleIdStart' => 'contract_sample_id_start',
        'startDate' => 'contract_start_date',
        'reportingDataRetention' => 'reporting_data_retention',
        'reportingDataTimeUnit' => 'reporting_data_time_unit',
        'hasDatabaseAgreement' => 'has_database_agreement'
    );

    my %lama_contracts_report_dict = (
        'reportGermline' => 'report_germline',
        'flagGermlineOnReport' => 'flagGermlineOnReport'
    );

    my %lama_content_translations_by_field_name = (
        'report_germline_level' => {
            'YES_ONLY_TREATABLE_FINDINGS' => '1: Behandelbare toevalsbevindingen',
            'YES_ALL_FINDINGS' => '2: Alle toevalsbevindingen',
            'NO' => '2: No',
            'YES' => '1: Yes'
        },
        'lab_status' => {
            'Processing' => 'In process',
        }
    );

    my %translations = (
        'CONT_CURR' => \%CONT_DICT,
        'SUBM_CURR' => \%SUBM_DICT,
        'SUBM_2023' => \%SUBM_DICT_2023,
        'SUBM_2022' => \%SUBM_DICT_2022,
        'SUBM_2021' => \%SUBM_DICT_2021,
        'SUBM_2020' => \%SUBM_DICT_2020,
        'SUBM_2019' => \%SUBM_DICT_2019,
        'SUBM_2018' => \%SUBM_DICT_2018,
        'SAMP_CURR' => \%SAMP_DICT,
        'SAMP_2023' => \%SAMP_DICT_2023,
        'SAMP_2022' => \%SAMP_DICT_2022,
        'SAMP_2021' => \%SAMP_DICT_2021,
        'SAMP_2020' => \%SAMP_DICT_2020,
        'SAMP_2019' => \%SAMP_DICT_2019,
        'SAMP_2018' => \%SAMP_DICT_2018,
        'PROC_CURR' => \%PROC_DICT,
        'PROC_2023' => \%PROC_DICT_2023,
        'PROC_2022' => \%PROC_DICT_2022,
        'lama_content_translations_by_field_name' => \%lama_content_translations_by_field_name,
        'lama_patient_dict' => \%lama_patient_dict,
        'lama_patient_tumor_sample_dict' => \%lama_patient_tumor_sample_dict,
        'lama_patient_tumor_sample_tumor_type_dict' => \%lama_patient_tumor_sample_tumor_type_dict,
        'lama_patient_tumor_sample_biopsy_dict' => \%lama_patient_tumor_sample_biopsy_dict,
        'lama_patient_reference_sample_dict' => \%lama_patient_reference_sample_dict,
        'lama_isolation_isolate_dict' => \%lama_isolation_isolate_dict,
        'lama_libraryprep_library_dict' => \%lama_libraryprep_library_dict,
        'lama_status_dict' => \%lama_status_dict,
        'lama_status_labworkflow_dict' => \%lama_status_labworkflow_dict,
        'lama_contracts_dict' => \%lama_contracts_dict,
        'lama_contracts_report_dict' => \%lama_contracts_report_dict,
    );

    return \%translations;
}

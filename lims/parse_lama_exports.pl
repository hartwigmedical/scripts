#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

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

use constant EMPTY => q{ };
use constant NACHAR => 'NA';

my %WARN_IF_ABSENT_FIELDS = (
    _id=>1,
    status=>1
);

my %lama_status_cohort_dict = (
    'cohortCode'           => 'cohort',
    'reportViral'          => 'report_viral',
    'reportGermline'       => 'report_germline',
    'flagGermlineOnReport' => 'flag_germline_on_report',
    'reportConclusion'     => 'report_conclusion',
    'reportPGX'            => 'report_pgx',
    'isShallowStandard'    => 'shallowseq',
    'addToDatabase'        => 'add_to_database',
    'addToDatarequests'    => 'add_to_datarequests',
    'sendPatientReport'    => 'send_patient_report',
    'hmfRegEnding'         => 'hmf_reg_ending'
);

my %lama_patient_dict = (
    '_id' => 'patient',
    'hospitalPatientId' => 'hospital_patient_id'
);

my %lama_patient_tumor_sample_dict = (
    'refFrBarcode'          => 'ref_sample_id',
    'sampleId'              => 'sample_name',
    'hospitalPaSampleId'    => 'hospital_pa_sample_id',
    'patientGermlineChoice' => 'report_germline_level',
    'primaryTumorType'      => 'ptum',
    'biopsySite'            => 'biopsy_site',
    'sopVersion'            => 'blood_registration_sop',
    'collectionDate'        => 'sampling_date',
    'isCUP'                 => 'is_cup',
    'arrivalHmf'            => 'arrival_date',
    'submissionNr'          => 'submission',
);

my %lama_patient_blood_sample_dict = (
    'sopVersion'      => 'blood_registration_sop',
    'collectionDate'  => 'sampling_date',
    'arrivalHmf'      => 'arrival_date',
    'sampleBarcode'   => 'sample_barcode',
    'sampleId'        => 'sample_name',
    'submissionNr'    => 'submission',
    'originalBarcode' => 'original_barcode'
);

my %lama_isolation_isolate_dict = (
    '_id'           => 'original_container_id',
    'isolationNr'   => 'isolation_id',
    'coupeBarcode'  => 'coupes_barcode',
    'sampleId'      => 'sample_name',
    'frBarcode'     => 'sample_id',
    'status'        => 'isolation_status',
    'type'          => 'isolation_type', # currently Blood|RNA|Tissue
    'concentration' => 'conc'
);

my %lama_libraryprep_library_dict = (
    '_id'          => 'sample_id',
    'isShallowSeq' => 'shallowseq',
    'yield'        => 'yield',
    'prepType'     => 'prep_type',
    'prepNr'       => 'prep_id',
    'sampleId'     => 'sample_name',
    'status'       => 'prep_status'
);

my %lama_status_dict = (
    '_id'                  => 'received_sample_id',
    'prepStatus'           => 'lab_status',
    'registrationDateTime' => 'registration_date_epoch',
    'sampleId'             => 'sample_name',
    'frBarcodeDNA'         => 'sample_id',
    'isTissue'             => 'is_tissue',
    'shallowPurity'        => 'shallow_purity',
    'finalPurity'          => 'purity'
);

my $input_dir = "./jsons";

my $lamaIsolationJson = $input_dir . '/Isolations.json';
my $lamaPatientJson = $input_dir . '/Patients.json';
my $lamaLibraryPrepJson = $input_dir . '/LibraryPreps.json';
my $lamaSampleStatusJson = $input_dir . '/SampleStatus.json';

sayInfo("Reading translation tables");
my $CNTR_TSV = './center2entity.tsv';
my $cntr_dict = parseDictFile( $CNTR_TSV, 'center2centername' );

sayInfo("Reading inputs");
my $lamaIsolation = readJson($lamaIsolationJson);
my $lamaPatient = readJson($lamaPatientJson);
my $lamaSampleStatus = readJson($lamaSampleStatusJson);
my $lamaLibraryPrep = readJson($lamaLibraryPrepJson);

my $statusObjects = parseLamaSampleStatus($lamaSampleStatus);
my $isolationObjects = parseLamaIsolation($lamaIsolation);
my $sampleObjects = parseLamaPatients($lamaPatient);
my $prepObjects = parseLamaLibraryPreps($lamaLibraryPrep);
my $submissionObjects = {};

my $samples = combineObjects($statusObjects, $sampleObjects, $isolationObjects, $prepObjects);
addMissingFieldsToSamples($samples, $submissionObjects, $cntr_dict);

printJson($samples, "./lama.json");

printJsonDebug($statusObjects, $sampleObjects, $isolationObjects, $prepObjects, "./lama_raw.json");

sub addMissingFieldsToSamples{
    my ($samples, $submissions, $centers_dict) = @_;

    while (my ($barcode, $sample) = each %$samples){
        my $name = $sample->{sample_name};

        my $analysis_type = NACHAR;
        my $submission = NACHAR;
        my $project_name = NACHAR;
        my $label = NACHAR;

        my ($patient_id, $study, $center, $tum_or_ref);
        my $name_regex = '^((CPCT|DRUP|WIDE|ACTN|CORE)[0-9A-Z]{2}([0-9A-Z]{2})\d{4})(T|R){1}';
        if ( $name =~ /$name_regex/ms ){
            ($patient_id, $study, $center, $tum_or_ref) = ($1, $2, $3, $4);
            $label = $study;
        }
        else{
            sayWarn("LAMA sample does ($name) because name does not fit $name_regex");
            $sample->{analysis_type} = 'GaveWarningUponImportingLama';
            next;
        }

        my $original_submission = $sample->{submission};
        my $isolation_type = $sample->{isolation_type};

        if ($isolation_type eq 'Tissue') {
            # this is DNA from tumor tissue
            $analysis_type = 'Somatic_T';
        }
        elsif ($isolation_type eq 'RNA') {
            # this is RNA from tumor tissue
            $analysis_type = 'RNAanalysis';
        }
        elsif ($isolation_type eq 'Blood') {
            # this is DNA from blood
            $analysis_type = 'Somatic_R';
        }
        elsif ($isolation_type eq 'Plasma') {
            # this is Plasma from blood
            $analysis_type = 'PlasmaAnalysis';
        }
        else{
            die "Should not happen: encountered unknown isolation type '$isolation_type' ($barcode)"
        }

        if ( $study eq 'CORE' and $name !~ /^COREDB/ ){
            ## specifically check for non-ref samples if submission is defined
            if ( $original_submission eq '' ){
                if ( $analysis_type ne 'Somatic_R' ){
                    sayWarn("SKIPPING CORE sample because of incorrect submission id \"$original_submission\" (id:$barcode name:$name)");
                }
                ## we only print warning for non R samples but skip them altogether
                next;
            }
            $sample->{ 'entity' } = $original_submission;
            $sample->{ 'project_name' } = $original_submission;

            ## Set the analysis type for CORE submissions to align with Excel LIMS samples
            if ( exists $submissions->{ $original_submission } ) {
                my $submission_object = $submissions->{ $original_submission };
                $project_name = $submission_object->{ 'project_name' };
                ## Reset project name for sample (from submission)
                $sample->{ 'project_name' } = $project_name;
                ## Add an analysis type to submission
                $submission->{ 'analysis_type' } = "OncoAct";
            }
            else{
                sayWarn("Unable to update submission \"$original_submission\" because not found in submissions (id:$barcode name:$name)");
            }
        }
        ## All other samples are clinical study based (CPCT/DRUP/WIDE/ACTN/COREDB)
        elsif ( exists $centers_dict->{ $center } ){
            my $centername = $centers_dict->{ $center };
            $sample->{ 'entity' } = join( "_", $study, $centername );
        }
        else {
            sayWarn("SKIPPING sample because is not CORE but center ID is unknown \"$center\" (id:$barcode name:$name)");
            next;
        }

        # Add the missing fields
        $sample->{analysis_type} = $analysis_type;
        $sample->{submission} = $submission;
        $sample->{original_submission} = $original_submission;
        $sample->{project_name} = $project_name;
        $sample->{label} = $label;
    }
}

sub combineObjects{
    my ($statuses, $samples, $isolations, $preps) = @_;
    my %lama_samples = ();

    my $missing_sample_count = 0;
    my $missing_isolate_count = 0;
    my $missing_prep_count = 0;

    while (my ($isolate_barcode, $object) = each %$statuses){
        my %sample_to_store = %{$object};
        my $sample_barcode = $sample_to_store{received_sample_id};

        # adding sample info to statuses
        if ( exists $samples->{$sample_barcode} ){
            #storeRecordByKey($samples->{$sample_barcode}, 'lama_sample', \%sample_to_store, "adding sample info for $isolate_barcode", 1);
            addRecordFieldsToTargetRecord($samples->{$sample_barcode}, \%sample_to_store, "merging in sample info for $isolate_barcode");
        }else{
            $missing_sample_count++;
        }

        # adding isolate info to statuses
        if ( exists $isolations->{$isolate_barcode} ){
            #storeRecordByKey($isolations->{$isolate_barcode}, 'lama_isolate', \%sample_to_store, "adding isolate info for $isolate_barcode", 1);
            addRecordFieldsToTargetRecord($isolations->{$isolate_barcode}, \%sample_to_store, "merging in isolate info for $isolate_barcode");
        }
        else{
            $missing_isolate_count++;
        }

        # adding prep info to statuses
        if ( exists $preps->{$isolate_barcode} ){
            #storeRecordByKey($preps->{$isolate_barcode}, 'lama_prep', \%sample_to_store, "adding prep info for $isolate_barcode", 1);
            addRecordFieldsToTargetRecord($preps->{$isolate_barcode}, \%sample_to_store, "merging in prep info for $isolate_barcode");
        }
        else{
            $missing_prep_count++;
        }

        storeRecordByKey(\%sample_to_store, $isolate_barcode, \%lama_samples, "final store $isolate_barcode", 1);
    }

    my $count = scalar keys %lama_samples;

    sayInfo("Total of $count samples stored from LAMA");
    sayInfo("  $missing_sample_count without sample info");
    sayInfo("  $missing_isolate_count without isolate info");
    sayInfo("  $missing_prep_count without prep info");

    return \%lama_samples;
}

sub printJson{
    my ($samples, $output_file) = @_;

    my %lama = ('samples' => $samples);
    my $coder = JSON::XS->new->utf8->canonical;
    my $lama_txt = $coder->encode(\%lama);

    open my $fh, '>', $output_file or die "Unable to open output file ($output_file): $!\n";
        print $fh $lama_txt;
    close $fh;
}

sub printJsonDebug{
    # TODO: can be removed after testing phase
    my ($statuses, $samples, $isolations, $preps, $lims_file) = @_;

    my $status_count = scalar keys %$statuses;
    my $sample_count = scalar keys %$samples;
    my $isolation_count = scalar keys %$isolations;
    my $prep_count = scalar keys %$preps;

    my %lims = ('sample' => $samples, 'status' => $statuses, 'isolation' => $isolations, 'prep' => $preps);
    my $coder = JSON::XS->new->utf8->canonical;
    my $lims_txt = $coder->encode(\%lims);

    sayInfo("Writing output to $lims_file:");
    sayInfo("  $sample_count samples");
    sayInfo("  $status_count statuses");
    sayInfo("  $isolation_count isolations");
    sayInfo("  $prep_count preps");

    sayInfo("Investigate with:");
    print " cat $lims_file | jq '.sample' | less\n";
    print " cat $lims_file | jq '.status' | less\n";
    print " cat $lims_file | jq '.isolation' | less\n";
    print " cat $lims_file | jq '.prep' | less\n";
    print " cat $lims_file | jq '.sample.FB05298354'\n";
    print " cat $lims_file | jq '.status.FR30543060'\n";
    print " cat $lims_file | jq '.isolation.FR30543060'\n";
    print " cat $lims_file | jq '.prep.FR30543060'\n";
    print " cat $lims_file | jq | less\n";

    open my $lims_json_fh, '>', $lims_file or die "Unable to open output file ($lims_file): $!\n";
        print $lims_json_fh $lims_txt;
    close $lims_json_fh;
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

    while (my ($src_key, $tgt_key) = each %$fieldsTranslationTable){
        if (exists $object->{$src_key}){
            $store->{$tgt_key} = $object->{$src_key};
        }
        elsif(exists $WARN_IF_ABSENT_FIELDS{$src_key}){
            sayWarn("No '$src_key' field in object ($info_tag)");
        }
        else{
            # do not warn about missing fields by default
        }
    }
}

sub epochToDate{
    # epoch time in milliseconds
    my ($epoch) = @_;
    my $registrationDate = strftime "%Y-%m-%d", localtime $epoch/1000;
    return $registrationDate;
}

sub parseLamaPatients {
    my ($objects) = @_;
    my %store = ();

    foreach my $patient (@$objects) {
        foreach my $sample (@{$patient->{tumorSamples}}) {
            processSampleOfPatient($patient, $sample, 'tumor', \%store);
        }
        foreach my $sample (@{$patient->{bloodSamples}}) {
            processSampleOfPatient($patient, $sample, 'blood', \%store);
        }
        # At time of writing no plasma samples are being used
        foreach my $sample (@{$patient->{plasmaSamples}}) {
             # do stuff
        }
    }
    return \%store;
}

sub processSampleOfPatient {
    my ($patient, $sample, $sample_origin, $store) = @_;
    my $sample_field_translations;
    my @sampleBarcodes;

    if( $sample_origin eq 'tumor' ){
        $sample_field_translations = \%lama_patient_tumor_sample_dict;
        @sampleBarcodes = @{$sample->{sampleBarcodes}};
    }
    elsif( $sample_origin eq 'blood' ){
        $sample_field_translations = \%lama_patient_blood_sample_dict;
        @sampleBarcodes = ($sample->{sampleBarcode});
    }
    else{
        my $sample_barcode = $sample->{sampleBarcode};
        die "Unknown sample origin provided to processSampleOfPatient ($sample_origin for $sample_barcode)\n";
    }

    my $info_tag = "patients->" . join("|", @sampleBarcodes);

    my %info = ();
    copyFieldsFromObject($sample, $info_tag, $sample_field_translations, \%info);
    copyFieldsFromObject($patient, $info_tag, \%lama_patient_dict, \%info);

    foreach my $sampleBarcode (@sampleBarcodes) {
        storeRecordByKey(\%info, $sampleBarcode, $store, "patient_samples");
    }
}

sub parseLamaLibraryPreps{
    my ($objects) = @_;
    my %store = ();

    foreach my $experiment (@$objects) {
        my $prep_sop = $experiment->{sopVersion};
        my $prep_id = $experiment->{_id};

        foreach my $object (@{$experiment->{libraries}}) {
            my $store_key = $object->{_id};
            my $status = $object->{status};

            # Only store prep info when OK
            next unless $status eq "Finished";

            my %info = ();
            my $info_tag = "libraryprep->$store_key";
            copyFieldsFromObject($object, $info_tag, \%lama_libraryprep_library_dict, \%info);
            storeRecordByKey(\%info, $store_key, \%store, "libraryprep", 1);
        }
    }
    return \%store;
}

sub parseLamaIsolation{
    my ($objects) = @_;
    my %store = ();

    foreach my $experiment (@$objects) {

        # When isolation experiment is Processing there are no frBarcodes yet so skip all isolates
        next if $experiment->{status} eq "Processing";

        foreach my $isolate (@{$experiment->{isolates}}) {

            # Once a isolation experiment is no longer Processing there should be a frBarcode
            if ( not exists $isolate->{frBarcode} ) {
                print Dumper $isolate;
                die "Isolate store key not defined for above object";
            }

            my $barcode = $isolate->{frBarcode};
            my $new_status = $isolate->{status};

            if ( exists $store{$barcode} ) {
                my $existing_status = $store{$barcode}{'isolation_status'};
                if ( $existing_status eq 'Finished' ){
                    if ($new_status eq 'Finished' ){
                        die "Should not happen: encountered duplicate Finished isolate $barcode ($new_status)";
                    }
                    else {
                        next;
                    }
                }
            }
            my %info = ();
            my $info_tag = "isolation->$barcode";
            copyFieldsFromObject($isolate, $info_tag, \%lama_isolation_isolate_dict, \%info);
            storeRecordByKey(\%info, $barcode, \%store, "isolation", 1);
        }
    }
    return \%store;
}

sub parseLamaSampleStatus{
    my ($objects) = @_;
    my %store = ();

    foreach my $object (@$objects){

        # there is no barcode in case a sample has not been prepped yet
        my $sampleBarcode = $object->{frBarcodeDNA};
        next unless defined $sampleBarcode;

        # debug
        my $debug_barcode = "";
        #my $debug_barcode = "FR30543342";
        if ( $debug_barcode ne "" and $sampleBarcode eq $debug_barcode ){
            print Dumper $object;
            die;
        }

        # Collect all info into one object
        my %status = ();
        my $info_tag = "samplestatus->$sampleBarcode";
        copyFieldsFromObject($object, $info_tag, \%lama_status_dict, \%status);
        copyFieldsFromObject($object->{cohort}, $info_tag, \%lama_status_cohort_dict, \%status);

        # Fix date fields
        $status{registration_date} = epochToDate($status{registration_date_epoch});

        # Store
        storeRecordByKey(\%status, $sampleBarcode, \%store, "patient->samples");

        # DEBUG
        #my @fields = sort (values %lama_status_dict, values %lama_status_cohort_dict);
        #say join(",", map($status{$_},@fields));
    }
    return \%store;
}

sub readJson{
    my ($json_file) = @_;
    sayInfo("  Parsing input json file $json_file");
    my $json_txt = read_file( $json_file );
    my $json_obj = decode_json( $json_txt );
    return( $json_obj );
}

sub sayInfo{
    my ($msg) = @_;
    say "[INFO] " . (strftime "%y%m%d %H:%M:%S", localtime) . " - " . $msg;
}

sub sayWarn{
    my ($msg) = @_;
    warn "[WARN] " . (strftime "%y%m%d %H:%M:%S", localtime) . " - " . $msg . "\n";
}

sub parseDictFile{
    my ($file, $fileType) = @_;
    sayInfo("  Parsing input dictionary file $file");
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
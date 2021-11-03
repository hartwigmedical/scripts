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

my %WARN_IF_ABSENT_FIELDS = (
    _id=>1,
    status=>1
);

my %fields_to_merge = (
    'sample' => {
        ''
    },
    'isolate' => {},
    'prep' => {},
);

my %lama_patient_dict = (
    '_id' => 'hmf_patient_id',
    'hospitalPatientId' => 'hospital_patient_id'
);

my %lama_patient_tumor_sample_dict = (
    'refFrBarcode'          => 'ref_sample_id',
    'patientId'             => 'patient_id',
    'sampleId'              => 'sample_name',
    'hospitalPaSampleId'    => 'hospital_patient_id',
    'patientGermlineChoice' => 'germline_choice',
    'remark'                => 'remark',
    'primaryTumorType'      => 'ptum',
    'biopsySite'            => 'biopsy_site',
    'sopVersion'            => 'blood_registration_sop',
    'collectionDate'        => 'sampling_date',
    'isCUP'                 => 'is_cup',
    'arrivalHmf'            => 'arrival_date',
    'submissionNr'          => 'submission',
);

my %lama_patient_blood_sample_dict = (
    'patientId'       => 'patient_id',
    'remark'          => 'remark',
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
    'coupeBarcode'  => 'coupe_barcode',
    'sampleId'      => 'sample_name',
    'frBarcode'     => 'sample_id',
    'isolationNr'   => 'isolation_id',
    'status'        => 'isolation_status',
    'concentration' => 'conc'
);

my %lama_libraryprep_library_dict = (
    '_id'                => 'sample_id',
    'startConcentration' => 'conc',
    'isShallowSeq'       => 'shallowseq',
    'yield'              => 'yield',
    'prepType'           => 'prep_type',
    'prepNr'             => 'prep_id',
    'sampleId'           => 'sample_name',
    'status'             => 'prep_status'
);

my %lama_status_dict = (
    '_id'                  => 'received_sample_id',
    'prepStatus'           => 'lab_status',
    'registrationDateTime' => 'registration_date_epoch',
    'sampleId'             => 'sample_name',
    'frBarcodeDNA'         => 'sample_id',
    'isTissue'             => 'is_tissue',
    'shallowPurity'        => 'shallow_purity',
    'finalPurity'          => 'final_purity'
);

my %lama_status_cohort_dict = (
    'cohortCode'        => 'cohort',
    'reportViral'       => 'report_viral',
    'reportGermline'    => 'report_germline',
    'reportPGX'         => 'report_pgx',
    'isShallowStandard' => 'shallowseq',
);


my $input_dir = "./jsons";

my $lamaIsolationJson = $input_dir . '/Isolations.json';
my $lamaPatientJson = $input_dir . '/Patients.json';
my $lamaLibraryPrepJson = $input_dir . '/LibraryPreps.json';
my $lamaSampleStatusJson = $input_dir . '/SampleStatus.json';

sayInfo("Reading inputs");
my $lamaIsolation = readJson($lamaIsolationJson);
my $lamaPatient = readJson($lamaPatientJson);
my $lamaSampleStatus = readJson($lamaSampleStatusJson);
my $lamaLibraryPrep = readJson($lamaLibraryPrepJson);

my $statusObjects = parseLamaSampleStatus($lamaSampleStatus);
my $isolationObjects = parseLamaIsolation($lamaIsolation);
my $sampleObjects = parseLamaPatients($lamaPatient);
my $prepObjects = parseLamaLibraryPreps($lamaLibraryPrep);

combineAndPrintJson($sampleObjects, $statusObjects, $isolationObjects, $prepObjects, "./lama.json");
printJson($sampleObjects, $statusObjects, $isolationObjects, $prepObjects, "./lama_raw.json");

sub combineAndPrintJson{
    my ($samples, $statuses, $isolations, $preps, $output_file) = @_;
    my %samples = ();

    my $missing_sample_count = 0;
    my $missing_isolate_count = 0;
    my $missing_prep_count = 0;

    while (my ($isolate_barcode, $object) = each %$statuses){
        my %sample_to_store = %{$object};
        my $sample_barcode = $sample_to_store{received_sample_id};

        if ( exists $samples->{$sample_barcode} ){
            storeRecordByKey($samples->{$sample_barcode}, 'lama_sample', \%sample_to_store, "adding sample info for $isolate_barcode", 1);
        }else{
            $missing_sample_count++;
        }

        if ( exists $isolations->{$isolate_barcode} ){
            storeRecordByKey($isolations->{$isolate_barcode}, 'lama_isolate', \%sample_to_store, "adding isolate info for $isolate_barcode", 1);
        }
        else{
            $missing_isolate_count++;
        }

        if ( exists $preps->{$isolate_barcode} ){
            storeRecordByKey($preps->{$isolate_barcode}, 'lama_prep', \%sample_to_store, "adding prep info for $isolate_barcode", 1);
        }
        else{
            $missing_prep_count++;
        }

        storeRecordByKey(\%sample_to_store, $isolate_barcode, \%samples, "final store $isolate_barcode", 1);
    }

    my $count = scalar keys %samples;
    my %lama = ('samples' => \%samples);
    my $coder = JSON::XS->new->utf8->canonical;
    my $lama_txt = $coder->encode(\%lama);

    open my $fh, '>', $output_file or die "Unable to open output file ($output_file): $!\n";
        print $fh $lama_txt;
    close $fh;

    sayInfo("Total of $count samples written to $output_file");
    sayInfo("  $missing_sample_count without sample info in LAMA");
    sayInfo("  $missing_isolate_count without isolate info in LAMA");
    sayInfo("  $missing_prep_count without prep info in LAMA");

}

sub printJson{
    my ($samples, $statuses, $isolations, $preps, $lims_file) = @_;
    my $sample_count = scalar keys %$samples;
    my $status_count = scalar keys %$statuses;
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
    # epoch time in miliseconds
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
        foreach my $sample (@{$patient->{plasmaSamples}}) {
            # At time of writing there are no plasma sample in LAMA but future it might
            processSampleOfPatient($patient, $sample, 'blood', \%store);
        }
    }
    return \%store;
}

sub processSampleOfPatient {
    my ($patient, $sample, $sample_type, $store) = @_;
    my $sample_field_translations;
    my @sampleBarcodes;

    if( $sample_type eq 'tumor' ){
        $sample_field_translations = \%lama_patient_tumor_sample_dict;
        @sampleBarcodes = @{$sample->{sampleBarcodes}};
    }
    elsif( $sample_type eq 'blood' ){
        $sample_field_translations = \%lama_patient_blood_sample_dict;
        @sampleBarcodes = ($sample->{sampleBarcode});
    }
    else{
        die "Unknown sample type provided to processSampleOfPatient ($sample_type)\n";
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
        foreach my $isolate (@{$experiment->{isolates}}) {
            my $store_key = $isolate->{frBarcode};
            my $status = $isolate->{status};

            if ( not defined $store_key ) {
                print Dumper $isolate;
                die "Isolate store key not defined for above object";
            }
            if ( exists $store{$store_key} ) {
                my $existing_isolation_status = $store{$store_key}{'isolation_status'};
                if ( $existing_isolation_status eq 'Finished' and $status eq 'Finished' ){
                    sayWarn("SKIPPING: Encountered duplicate Finished isolate $store_key ($status)");
                    next;
                }
            }

            my %info = ();
            my $info_tag = "isolation->$store_key";
            copyFieldsFromObject($isolate, $info_tag, \%lama_isolation_isolate_dict, \%info);
            storeRecordByKey(\%info, $store_key, \%store, "isolation", 1);
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

sub parseLamaCohort{
    my ($objects) = @_;
    foreach my $object (@$objects){
        my $isShallowStandard = $object->{isShallowStandard};
        my $reportPGX = $object->{reportPGX};
        sayInfo("isShallow = " . $isShallowStandard);
        sayInfo("reportPGX = " . $reportPGX);
    }
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
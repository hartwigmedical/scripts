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

my $input_dir = "./lama_export";

my $lamaIsolationJson = $input_dir . '/Isolations.json';
my $lamaPatientJson = $input_dir . '/Patients.json';
my $lamaLibraryPrepJson = $input_dir . '/LibraryPreps.json';
my $lamaSampleStatusJson = $input_dir . '/SampleStatus.json';
#my $lamaCohortJson = $input_dir . '/Cohorts.json';
#my $lamaHospitalJson = $input_dir . '/Hospitals.json';
#my $lamaLogBookJson = $input_dir . '/LogBooks.json';
#my $lamaSnpCheckJson = $input_dir . '/SnpChecks.json';

my $lamaIsolation = readJson($lamaIsolationJson);
my $lamaPatient = readJson($lamaPatientJson);
my $lamaSampleStatus = readJson($lamaSampleStatusJson);
my $lamaLibraryPrep = readJson($lamaLibraryPrepJson);
#my $lamaCohort = ReadJson($lamaCohortJson);

my $statusObjects = parseLamaSampleStatus($lamaSampleStatus);
my $isolationObjects = parseLamaIsolation($lamaIsolation);
my $sampleObjects = parseLamaPatients($lamaPatient);

printJson($sampleObjects, "./lama.json");
#print Dumper $sampleObjects;

#dumpRecordsOfArray($lamaSampleStatus, "Status");
#dumpRecordsOfArray($lamaCohort, "Cohort");
#dumpRecordsOfArray($lamaIsolation, "Isolation");
#dumpRecordsOfArray($lamaPatient, "Patient");

sub storeRecordByKey{
    my ($record, $key, $store, $info_tag) = @_;
    warn "Store key ($key) already exists in store ($info_tag)\n" if exists $store->{$key};
    my %copy_of_record = %$record;
    $store->{$key} = \%copy_of_record;
}

sub copyFieldsFromObject{
    my ($object, $tag, $fieldsTranslationTable, $store) = @_;
    while (my ($src_key, $tgt_key) = each %$fieldsTranslationTable){
        if (exists $object->{$src_key}){
            $store->{$tgt_key} = $object->{$src_key};
        }
        else{
            warn "No $src_key field in object ($tag)\n";
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

    my %lama_patient_tumor_sample_dict = (
        'refFrBarcode'          => 'sample_id',
        'patientId'             => 'patient_id',
        'hospitalPaSampleId'    => 'hospital_patient_id',
        'patientGermlineChoice' => 'germline_choice',
        'remark'                => 'remark',
        'primaryTumorType'      => 'ptum',
        'sopVersion'            => 'sop',
        'collectionDate'        => 'collection_date',
        'isCUP'                 => 'is_cup',
        'arrivalHmf'            => 'arrival_date'
    );

    foreach my $object (@$objects) {

        my $patientId = $object->{_id};
        my $hospitalPatientId = $object->{hospitalPatientId};
        my $bloodSamples = $object->{bloodSamples};
        my $tumorSamples = $object->{tumorSamples};
        my $plasmaSamples = $object->{plasmaSamples};

        foreach my $sample (@$tumorSamples) {
            my $sampleBarcode = $sample->{refFrBarcode};
            next if $sampleBarcode eq "";

            my %info = ();
            copyFieldsFromObject($sample, $sampleBarcode, \%lama_patient_tumor_sample_dict, \%info);
            storeRecordByKey(\%info, $sampleBarcode, \%store, "patient_samples");

            #warn "Key ($sampleBarcode) already exists in PatientSamples store\n" if exists $store{$sampleBarcode};
            #$store{$sampleBarcode} = \%info;
        }
    }
    return \%store;
}


sub parseLamaIsolation{
    my ($objects) = @_;
    my %store = ();

    my %lama_isolation_isolate_dict = (
        '_id'         => 'original_container_id',
        'isolationNr' => 'isolation_id'
    );

    foreach my $experiment (@$objects) {
        my $sop = $experiment->{sopVersion};

        foreach my $isolate (@{$experiment->{isolates}}) {
            my $store_key = $isolate->{_id};
            my %info = ('sop' => $sop);
            copyFieldsFromObject($isolate, $store_key, \%lama_isolation_isolate_dict, \%info);

            # TODO: check if we need to diff dates to store most recent isolation info
            #warn "[WARN] Isolation key will be overwritten ($store_key)\n" if exists $store{$store_key};
            $store{$store_key} = \%info;
        }
    }
    return \%store;
}

sub parseLamaSampleStatus{
    my ($objects) = @_;
    my %store = ();

    my %lama_status_dict = (
        '_id' => 'received_sample_id',
        'prepStatus' => 'lab_status',
        'registrationDateTime' => 'registration_date_epoch',
        'sampleId' => 'sample_name',
        'frBarcodeDNA' => 'sample_id',
        'isTissue' => 'is_tissue'
    );

    my %lama_status_cohort_dict = (
        'cohortCode'        => 'cohort',
        'reportViral'       => 'report_viral',
        'reportGermline'    => 'report_germline',
        'reportPGX'         => 'report_pgx',
        'isShallowStandard' => 'shallowseq',
    );

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
        copyFieldsFromObject($object, $sampleBarcode, \%lama_status_dict, \%status);
        copyFieldsFromObject($object->{cohort}, $sampleBarcode, \%lama_status_cohort_dict, \%status);

        # Fix date fields
        $status{registration_date} = epochToDate($status{registration_date_epoch});

        # Output
        #my @fields = sort (values %lama_status_dict, values %lama_status_cohort_dict);
        #say join(",", map($status{$_},@fields));

        # Store
        die "Key ($sampleBarcode) already exists in Status store\n" if exists $store{$sampleBarcode};
        $store{$sampleBarcode} = \%status;
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

sub printJson{
    my ($samples, $lims_file) = @_;
    my $samp_count = scalar keys %$samples;

    my %lims = ('samples' => $samples);
    my $coder = JSON::XS->new->utf8->canonical;
    my $lims_txt = $coder->encode(\%lims);

    sayInfo("  Writing output to $lims_file ($samp_count samples)");
    open my $lims_json_fh, '>', $lims_file or die "Unable to open output file ($lims_file): $!\n";
    print $lims_json_fh $lims_txt;
    close $lims_json_fh;
}

sub sayInfo{
    my ($msg) = @_;
    say "[INFO] " . (strftime "%y%m%d %H:%M:%S", localtime) . " - " . $msg;
}

sub dumpRecordsOfArray{
    my ($arrayref, $tag) = @_;
    sayInfo("  Dump first record of $tag");
    print Dumper $arrayref->[1];
    sayInfo("  Dump -100 record of $tag");
    print Dumper $arrayref->[-100];
    sayInfo("  Dump last record of $tag");
    print Dumper $arrayref->[-1];
}
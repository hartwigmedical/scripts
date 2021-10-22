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

my %DO_NOT_EXIST_WARNING_FIELDS_TO_SKIP = (
    startConcentration=>1,
    patientGermlineChoice=>1,
    primaryTumorType=>1,
    biopsySite=>1,
    shallowPurity=>1,
    finalPurity=>1,
    coupeBarcode=>1,
    concentration=>1
);

my %lama_patient_tumor_sample_dict = (
    'refFrBarcode'          => 'ref_sample_id',
    'patientId'             => 'patient_id',
    'hospitalPaSampleId'    => 'hospital_patient_id',
    'patientGermlineChoice' => 'germline_choice',
    'remark'                => 'remark',
    'primaryTumorType'      => 'ptum',
    'biopsySite'            => 'biopsy_site',
    'sopVersion'            => 'sop',
    'collectionDate'        => 'sampling_date',
    'isCUP'                 => 'is_cup',
    'arrivalHmf'            => 'arrival_date'
);

my %lama_isolation_isolate_dict = (
    '_id'           => 'original_container_id',
    'isolationNr'   => 'isolation_id',
    'coupeBarcode'  => 'coupe_barcode',
    'sampleId'      => 'sample_name',
    'frBarcode'     => 'sample_id',
    'isolationNr'   => 'isolation_id',
    'concentration' => 'conc'
);

my %lama_libraryprep_libraries_dict = (
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
my $prepObjects = parseLamaLibraryPreps($lamaLibraryPrep);

printJson($sampleObjects, $statusObjects, $isolationObjects, $prepObjects, "./lama.json");
#print Dumper $prepObjects;

#dumpRecordsOfArray($lamaSampleStatus, "SampleStatus");
#dumpRecordsOfArray($lamaCohort, "Cohort");
#dumpRecordsOfArray($lamaIsolation, "Isolation");
#dumpRecordsOfArray($lamaPatient, "Patient");
#dumpRecordsOfArray($lamaLibraryPrep, "LibraryPrep");

sub printJson{
    my ($samples, $statuses, $isolations, $preps, $lims_file) = @_;
    my $sample_count = scalar keys %$samples;
    my $status_count = scalar keys %$statuses;
    my $isolation_count = scalar keys %$isolations;
    my $prep_count = scalar keys %$preps;

    my %lims = ('sample' => $samples, 'status' => $statuses, 'isolation' => $isolations, 'prep' => $preps);
    my $coder = JSON::XS->new->utf8->canonical;
    my $lims_txt = $coder->encode(\%lims);

    sayInfo(" Writing output to $lims_file:");
    sayInfo("   $sample_count samples");
    sayInfo("   $status_count statuses");
    sayInfo("   $isolation_count isolations");
    sayInfo("   $prep_count preps");
    sayInfo(" Investigate with:");
    print " cat $lims_file | jq '.sample' | less\n";
    print " cat $lims_file | jq '.status' | less\n";
    print " cat $lims_file | jq '.isolation' | less\n";
    print " cat $lims_file | jq '.prep' | less\n";
    print " cat $lims_file | jq '.sample.FB05298354'\n";
    print " cat $lims_file | jq '.status.FR30543060'\n";

    open my $lims_json_fh, '>', $lims_file or die "Unable to open output file ($lims_file): $!\n";
    print $lims_json_fh $lims_txt;
    close $lims_json_fh;
}

sub storeRecordByKey{
    my ($record, $key, $store, $info_tag) = @_;
    sayWarn("Store key already exists in store and will be overwritten ($key for $info_tag)") if exists $store->{$key};
    my %copy_of_record = %$record;
    $store->{$key} = \%copy_of_record;
}

sub copyFieldsFromObject{
    my ($object, $info_tag, $fieldsTranslationTable, $store) = @_;

    while (my ($src_key, $tgt_key) = each %$fieldsTranslationTable){
        if (exists $object->{$src_key}){
            $store->{$tgt_key} = $object->{$src_key};
        }
        elsif(exists $DO_NOT_EXIST_WARNING_FIELDS_TO_SKIP{$src_key}){
            # ignore the fact field is missing
        }
        else{
            sayWarn("No '$src_key' field in object ($info_tag)");
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

    foreach my $object (@$objects) {

        my $patientId = $object->{_id};
        my $hospitalPatientId = $object->{hospitalPatientId};
        my $bloodSamples = $object->{bloodSamples};
        my $tumorSamples = $object->{tumorSamples};
        my $plasmaSamples = $object->{plasmaSamples};

        foreach my $sample (@$tumorSamples) {
            my @sampleBarcodes = @{$sample->{sampleBarcodes}};
            my $info_tag = "patients->" . join("|", @sampleBarcodes);

            my %info = ();
            copyFieldsFromObject($sample, $info_tag, \%lama_patient_tumor_sample_dict, \%info);
            foreach my $sampleBarcode (@sampleBarcodes){
                storeRecordByKey(\%info, $sampleBarcode, \%store, "patient_samples");
            }
        }
    }
    return \%store;
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
            copyFieldsFromObject($object, $info_tag, \%lama_libraryprep_libraries_dict, \%info);

            # TODO: check if we need to diff dates to store most recent library prep info
            storeRecordByKey(\%info, $store_key, \%store, "libraryprep");
        }
    }
    return \%store;
}

sub parseLamaIsolation{
    my ($objects) = @_;
    my %store = ();

    foreach my $experiment (@$objects) {
        my $sop = $experiment->{sopVersion};

        foreach my $isolate (@{$experiment->{isolates}}) {
            my $store_key = $isolate->{frBarcode};
            my $status = $isolate->{status};

            # Only store info when OK
            next unless $status eq "Finished";

            my %info = ();
            my $info_tag = "isolation->$store_key";
            copyFieldsFromObject($isolate, $info_tag, \%lama_isolation_isolate_dict, \%info);

            # TODO: check if we need to diff dates to store most recent isolation info
            storeRecordByKey(\%info, $store_key, \%store, "isolation");
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

sub dumpRecordsOfArray{
    my ($arrayref, $tag) = @_;
    sayInfo("  Dump first record of $tag");
    print Dumper $arrayref->[1];
    sayInfo("  Dump -100 record of $tag");
    print Dumper $arrayref->[-100];
    sayInfo("  Dump last record of $tag");
    print Dumper $arrayref->[-1];
}
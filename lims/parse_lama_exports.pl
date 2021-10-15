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

my $lamaCohortJson = $input_dir . '/Cohorts.json';
my $lamaHospitalJson = $input_dir . '/Hospitals.json';
my $lamaIsolationJson = $input_dir . '/Isolations.json';
my $lamaLogBookJson = $input_dir . '/LogBooks.json';
my $lamaPatientJson = $input_dir . '/Patients.json';
my $lamaLibraryPrepJson = $input_dir . '/LibraryPreps.json';
my $lamaSnpCheckJson = $input_dir . '/SnpChecks.json';
my $lamaSampleStatusJson = $input_dir . '/SampleStatus.json';

my $lamaCohort = readJson($lamaCohortJson);
my $lamaIsolation = readJson($lamaIsolationJson);
my $lamaPatient = readJson($lamaPatientJson);
my $lamaLibraryPrep = readJson($lamaLibraryPrepJson);
my $lamaSampleStatus = readJson($lamaSampleStatusJson);
#print Dumper $lamaCohort;

sayInfo("First Objects:");
#printFirstOfArray($lamaCohort, "Cohort");
#printFirstOfArray($lamaIsolation, "Isolation");
#printFirstOfArray($lamaPatient, "Patient");
parseLamaSampleStatus($lamaSampleStatus);
printFirstOfArray($lamaSampleStatus, "Status");

sub epochToDate{
    # epoch time in miliseconds
    my ($epoch) = @_;
    my $registrationDate = strftime "%Y-%m-%d", localtime $epoch/1000;
    return $registrationDate;
}

sub parseLamaSampleStatus{
    my ($objects) = @_;

    my %lama_status_dict = (
        'prepStatus' => 'lab_status',
        'registrationDateTime' => 'registration_date_epoch',
        'sampleId' => 'sample_name',
        'frBarcodeDNA' => 'sample_id'
    );

    my %lama_status_cohort_dict = (
        'cohortCode' => 'cohort',
    );

    foreach my $object (@$objects){

        # there is no barcode in case a sample has not been prepped yet
        my $sampleBarcode = $object->{frBarcodeDNA};
        next unless defined $sampleBarcode;

        # debug
        #if ( $sampleBarcode eq 'FR12244536' ){
        #    print Dumper $object;
        #    die;
        #}

        my %status = ();
        while (my ($src_key, $tgt_key) = each %lama_status_dict){
            $status{$tgt_key} = $object->{$src_key};
        }
        while (my ($src_key, $tgt_key) = each %lama_status_cohort_dict){
            $status{$tgt_key} = $object->{cohort}->{$src_key};
        }
        $status{registration_date} = epochToDate($status{registration_date_epoch});
        #say join(" ", $sampleName, $sampleBarcode, $cohortCode, $regDate, $labStatus);
        print Dumper \%status;
    }
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

sub printFirstOfArray{
    my ($arrayref, $tag) = @_;
    sayInfo("  Print first object for $tag");
    print Dumper $arrayref->[0];
    #print Dumper $arrayref;
}
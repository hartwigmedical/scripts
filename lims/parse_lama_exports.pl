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
printFirstOfArray($lamaCohort, "Cohort");
printFirstOfArray($lamaIsolation, "Isolation");
printFirstOfArray($lamaPatient, "Patient");
#parseLama($lamaCohort);

sub parseLama{
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
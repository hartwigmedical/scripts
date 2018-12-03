#!/usr/bin/env perl
use strict;
use warnings;
use 5.016_000;

use Getopt::Long;
use File::Basename;
use Data::Dumper;

my $healthcheck_stdout_log;

GetOptions("healthcheck-log-file|f=s" => \$healthcheck_stdout_log)
    or die("Error in command line arguments\n");
say "Run with: --healthcheck-log-file /path/to/health-checks.txt" and exit(1) unless $healthcheck_stdout_log;
die "[ERROR] Provided healthcheck stdout log file ($healthcheck_stdout_log) does not exist\n" unless -f $healthcheck_stdout_log;

my $hcdata = parseHealthCheckerOutput($healthcheck_stdout_log);
my $sampleCount = $hcdata->{'sampleCount'};

if ($sampleCount == 1) {
    say "[INFO] Starting quality check (single sample mode)";
    doSingleSampleTests($hcdata);
} elsif ($sampleCount == 2) {
    say "[INFO] Starting quality check (somatic mode)";
    doSomaticTests($hcdata);
} else {
    die "[ERROR] Incompatible number of samples ($sampleCount)\n";
}

sub doSingleSampleTests {
    my ($hcdata) = @_;

    my $sample = $hcdata->{'refSample'};
    die "[ERROR] REF sample not determined\n" unless $sample;

    my $COV_PCT_10X_R = getValueBySample($hcdata, $sample, 'COVERAGE_10X');
    my $COV_PCT_20X_R = getValueBySample($hcdata, $sample, 'COVERAGE_20X');

    say "[INFO] Info Summary:";
    say "[INFO]   SAMPLE = $sample";
    say "[INFO] QC Tests:";

    my $fails = 0;

    ifLowerFail(\$fails, 0.90, $COV_PCT_10X_R, 'COVERAGE_10X_R');
    ifLowerFail(\$fails, 0.70, $COV_PCT_20X_R, 'COVERAGE_20X_R');

    generateFinalResultMessage($sample, $fails);
}

sub doSomaticTests {
    my ($hcdata) = @_;

    my @samples = @{$hcdata->{'sampleNames'}};
    my $refSample = $hcdata->{'refSample'};
    my $tumorSample = $hcdata->{'tumorSample'};
    die "[ERROR] REF sample not determined\n" unless $refSample;
    die "[ERROR] TUMOR sample not determined\n" unless $tumorSample;

    my $SAMPLE_NAMES = join(", ", @samples);
    my $SOM_SNP_COUNT = getValueBySample($hcdata, $tumorSample, 'SOMATIC_SNP_COUNT');
    my $SOM_SNP_DBSNP_COUNT = getValueBySample($hcdata, $tumorSample, 'SOMATIC_SNP_DBSNP_COUNT');
    my $AMBER_MEAN_BAF = getValueBySample($hcdata, $tumorSample, 'MEAN_BAF');
    my $PURPLE_QC_STATUS = getValueBySample($hcdata, $tumorSample, 'QC_STATUS');
    my $COV_PCT_10X_R = getValueBySample($hcdata, $refSample, 'COVERAGE_10X');
    my $COV_PCT_20X_R = getValueBySample($hcdata, $refSample, 'COVERAGE_20X');
    my $COV_PCT_30X_T = getValueBySample($hcdata, $tumorSample, 'COVERAGE_30X');
    my $COV_PCT_60X_T = getValueBySample($hcdata, $tumorSample, 'COVERAGE_60X');
    my $SOM_SNP_DBSNP_PCT = roundNumber($SOM_SNP_DBSNP_COUNT * 100 / $SOM_SNP_COUNT);

    say "[INFO] Info Summary:";
    say "[INFO]   SAMPLES = $SAMPLE_NAMES";
    say "[INFO]   SOMATIC_SNP_COUNT = " . commify($SOM_SNP_COUNT);
    say "[INFO]   SOMATIC_SNP_DBSNP_COUNT = " . commify($SOM_SNP_DBSNP_COUNT) . " (" . $SOM_SNP_DBSNP_PCT . "% of all SOMATIC_SNP)";
    say "[INFO] QC Tests:";

    my $fails = 0;

    my $test = 'DBSNP_CONTAMINATION';
    my $som_check1_fail = $SOM_SNP_DBSNP_COUNT > 250000;
    if ($som_check1_fail) {
        printMsg('FAIL', $test) and $fails++;
    } else {
        printMsg('INFO', "  [OK] $test");
    }

    $test = 'NON_DBSNP_CONTAMINATION';
    my $som_check2_fail = $SOM_SNP_COUNT > 1000000;
    my $som_check3_fail = $SOM_SNP_DBSNP_COUNT / $SOM_SNP_COUNT < 0.2;
    if ($som_check2_fail and $som_check3_fail) {
        printMsg('FAIL', $test) and $fails++;
    } else {
        printMsg('INFO', "  [OK] $test");
    }

    $test = 'PURPLE_QC';
    if ($PURPLE_QC_STATUS ne "PASS") {
        printMsg('FAIL', "$test ($PURPLE_QC_STATUS)") and $fails++;
    } else {
        printMsg('INFO', "  [OK] $test ($PURPLE_QC_STATUS)");
    }

    # KODU:
    # Complete sample swap  =>  <0.4
    # 10%+ contamination in tumor =>  < 0.486
    # Potential Tumor contamination in Ref  => >0.5
    ifLowerFail(\$fails, 0.48, $AMBER_MEAN_BAF, 'AMBER_MEAN_BAF');
    ifHigherFail(\$fails, 0.6, $AMBER_MEAN_BAF, 'AMBER_MEAN_BAF');

    ifLowerFail(\$fails, 0.90, $COV_PCT_10X_R, 'COVERAGE_10X_R');
    ifLowerFail(\$fails, 0.70, $COV_PCT_20X_R, 'COVERAGE_20X_R');
    ifLowerFail(\$fails, 0.80, $COV_PCT_30X_T, 'COVERAGE_30X_T');
    ifLowerFail(\$fails, 0.65, $COV_PCT_60X_T, 'COVERAGE_60X_T');

    generateFinalResultMessage($tumorSample, $fails);
}

sub generateFinalResultMessage {
    my ($sample, $fails) = @_;

    my $final_status = 'FAIL';
    $final_status = 'OK' if $fails == 0;
    my $final_msg = "TEST RESULT for $sample (fails:$fails) = $final_status";

    if ($final_status ne 'OK') {
        say "[FAIL] $final_msg";
    } else {
        say "[INFO] $final_msg";
    }
}

sub ifLowerFail {
    my ($totalFailCount, $limit, $value, $msg) = @_;
    my $fail = 0;
    $fail = $value < $limit unless $limit eq 'NA';

    if ($fail) {
        printMsg('FAIL', "$msg: $value < $limit") and $$totalFailCount++;
    } else {
        printMsg('INFO', "  [OK] $msg: $value > $limit");
    }
}

sub ifHigherFail {
    my ($totalFailCount, $limit, $value, $msg) = @_;
    my $fail = 0;
    $fail = $value > $limit unless $limit eq 'NA';

    if ($fail) {
        printMsg('FAIL', "$msg: $value > $limit") and $$totalFailCount++;
    } else {
        printMsg('INFO', "  [OK] $msg: $value < $limit");
    }
}

sub getValueBySample {
    my ($hcdata, $sample, $key, $key2) = @_;

    if (defined $hcdata->{$sample}{$key}) {
        my $return = $hcdata->{$sample}{$key};
        die "[ERROR] Key \"$key\" found for sample \"$sample\" but is ERROR\n" if $return eq 'ERROR';
        return ($return);
    } elsif ($key2 and defined $hcdata->{$sample}{$key2}) {
        my $return = $hcdata->{$sample}{$key2};
        die "[ERROR] Key \"$key2\" found for sample \"$sample\" but is ERROR\n" if $return eq 'ERROR';
        return ($return);
    } else {
        die "[ERROR] Key \"$key\" not found for sample \"$sample\"\n";
    }
}

sub printMsg {
    my ($type, $msg) = @_;
    say "[$type] $msg";
}

sub parseHealthCheckerOutput {
    my %hcdata = ();
    my ($file) = @_;
    open FILE, "<", $file or die "Unable to open \"$file\"\n";
    while (<FILE>) {
        chomp;
        if ($_ =~ /Check '([^']+)' for sample '([^']+)' has value '([^']+)'/) {
            my $key = $1;
            my $sam = $2;
            my $val = $3;
            $hcdata{$sam}{$key} = $val;
        }
    }
    close FILE;
    my @samples = keys %hcdata;
    my $sampleCount = scalar @samples;

    $hcdata{sampleNames} = \@samples;
    $hcdata{sampleCount} = $sampleCount;

    die "[ERROR] no samples detected in log file ($file)\n" if $sampleCount == 0;

    # (THHO) if somatic run: determine tum and ref sample (perhaps better to define in healthcheck output?)
    if ($sampleCount == 2) {
        my $refSample;
        my $tumorSample;

        foreach my $sample (@samples) {
            if (exists $hcdata{$sample}{'SOMATIC_SNP_COUNT'}) {
                $tumorSample = $sample;
            } else {
                $refSample = $sample;
            }
        }
        die "[ERROR] Multiple samples found so Somatic mode, but ref sample missing?\n" if not defined $refSample;
        die "[ERROR] Multiple samples found so Somatic mode, but tumor sample missing?\n" if not defined $tumorSample;

        $hcdata{refSample} = $refSample;
        $hcdata{tumorSample} = $tumorSample;
    } elsif ($sampleCount == 1) {
        $hcdata{refSample} = $samples[0];
    }

    return (\%hcdata);
}

sub roundNumber {
    my ($number) = @_;
    my $rounded = sprintf("%." . '2' . "f", $number);
    return ($rounded);
}

sub commify {
    local $_ = shift;
    1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
    return $_;
}


#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

use Data::Dumper;
use Getopt::Long;
use File::Basename;
use 5.010.000;

my $MIN_DISTANCE = 2;
my $script_name = basename $0;
my $script_dir = dirname $0;
my $out_file = "index_comparison.tsv";
my $fail_str = 'FAIL';
my $success_str = 'OK';

my $HELP =<<HELP;

  Description
    Compares two dual barcode index collections (format: name1<tab>name2<tab>seq1<tab>seq2)
    and determines the i7 and i5 distance for every possible combination. Result for each is
    '$fail_str' in case both i5 and i7 distances are lower than $MIN_DISTANCE, and otherwise '$success_str'.
  Usage
    $script_name <tsv-with-external-indexes> <tsv-with-hmf-indexes>
    $script_name <tsv-with-external-indexes>
  Example TSV input
    #i7name i5name i7seq i5seq
    IDT8_i7_1 IDT8_I5_1 CTGATCGT GCGCATAT
  Example output (written to $out_file)
    #result dist1 dist2 distance idx1_tag1 idx1_seq1 idx1_tag2 idx1_seq2 idx2_tag1 idx2_seq1 idx2_tag2 idx2_seq2
    OK 5 3 8 F498_LCH588 CTGCGGAT CZ-0379-LL_A AGATAACC IDT8_i7_384 CGACCATT IDT8_I5_384 TGATAGGC

HELP
print $HELP and exit(0) if scalar @ARGV == 0 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help';

my $file1 = "./indexes_customer.tsv";
my $file2 = "$script_dir/indexes_hmf.tsv";
$file1 = $ARGV[0] if $ARGV[0];
$file2 = $ARGV[1] if $ARGV[1];

say "[INFO] Comparing $file1 with $file2";
my $indexes1 = read_indexes_tsv($file1);
my $indexes2 = read_indexes_tsv($file2);

my @header = qw(result dist1 dist2 distance idx1_tag1 idx1_seq1 idx1_tag2 idx1_seq2 idx2_tag1 idx2_seq1 idx2_tag2 idx2_seq2);
open(my $OUT_FH, '>', $out_file) or die $!;
say $OUT_FH '##' . " Input file1: $file1";
say $OUT_FH '##' . " Input file2: $file2";
say $OUT_FH '##' . " Minimal distance: $MIN_DISTANCE";
say $OUT_FH '#' . join("\t", @header);
foreach my $index1 (@$indexes1){
    foreach my $index2 (@$indexes2){
        compare($index1, $index2, $OUT_FH);
    }
}
close $OUT_FH;
say "[INFO] Finished with $script_name (see $out_file for result).";

sub compare{
    my ($index1, $index2, $out_fh) = @_;

    my $i1_t1 = $index1->{tag1};
    my $i1_t2 = $index1->{tag2};
    my $i2_t1 = $index2->{tag1};
    my $i2_t2 = $index2->{tag2};

    my $i1_s1 = $index1->{seq1};
    my $i1_s2 = $index1->{seq2};
    my $i2_s1 = $index2->{seq1};
    my $i2_s2 = $index2->{seq2};

    my $tag = join("_", $i1_t1, $i1_t2, "vs", $i2_t1, $i2_t2);

    die "Impossible to compare $i1_s1 with $i2_s1 (reason: different size)\n" unless length($i1_s1) eq length($i2_s1);
    die "Impossible to compare $i1_s2 with $i2_s2 (reason: different size)\n" unless length($i1_s2) eq length($i2_s2);

    my $result = "";
    my $distance1 = calc_distance_between_two_seq($i1_s1, $i2_s1);
    my $distance2 = calc_distance_between_two_seq($i1_s2, $i2_s2);

    if ($distance1 < $MIN_DISTANCE and $distance2 < $MIN_DISTANCE){
        warn "[FAIL] distance too small for $tag ($distance1 and $distance2)\n";
        $result = $fail_str;
    }else{
        $result = $success_str;
    }
    my @fields = join("\t", $i1_t1, $i1_s1, $i1_t2, $i1_s2, $i2_t1, $i2_s1, $i2_t2, $i2_s2);
    say $out_fh join("\t", $result, $distance1, $distance2, $distance1+$distance2, @fields);
}

sub calc_distance_between_two_seq{
    my ($seq1, $seq2) = @_;
    my $distance = ( $seq1 ^ $seq2 ) =~ tr/\0//c;
    return $distance;
}

sub read_indexes_tsv{
    my ($file) = @_;
    say "[INFO] Reading file $file";
    my @indexes = ();
    open my $fh, '<', $file or die "Unable to open file ($file): $!\n";
    while ( <$fh> ) {
        chomp;
        next if $_ eq '' or $_ =~ /^#/;
        my ($tag1, $tag2, $seq1, $seq2) = split("\t", $_);
        my %index = (tag1 => $tag1, tag2 => $tag2, seq1 => $seq1, seq2 => $seq2);
        push(@indexes, \%index);
    }
    return \@indexes;
}
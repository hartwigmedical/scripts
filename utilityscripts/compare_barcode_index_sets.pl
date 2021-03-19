#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

use Data::Dumper;
use Getopt::Long;
use File::Basename;
use 5.010.000;

use constant MIN_DISTANCE => 2;

my $HELP =<<HELP;
  Description
    Takes two TSV dual barcode index collection TSV files (columns: name1, name2, seq1, seq2)
    and determines the i7 and i5 distance for every dual index in file1 against every in file2.
HELP
print $HELP and exit(0) if scalar @ARGV > 0 and ($ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $file1 = "index_customer.tsv";
my $file2 = "index_hmf.tsv";
my $out_file = "table.tsv";

say "Comparing $file1 with $file2";
my $indexes1 = read_indexes_tsv($file1);
my $indexes2 = read_indexes_tsv($file2);

my @header = qw(result dist_s1 dist_s2 dist_total i1_t1 i1_s1 i1_t2 i1_s2 i2_t1 i2_s1 i2_t2 i2_s2);
open(my $OUT_FH, '>', $out_file) or die $!;
say $OUT_FH '#' . join("\t", @header);
foreach my $index1 (@$indexes1){
    foreach my $index2 (@$indexes2){
        compare($index1, $index2, $OUT_FH);
    }
}
close $OUT_FH;
say "Done. See $out_file for complete result.";

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

    die "Cannot compare $i1_s1 with $i2_s1 (reason: different size)\n" unless length($i1_s1) eq length($i2_s1);
    die "Cannot compare $i1_s2 with $i2_s2 (reason: different size)\n" unless length($i1_s2) eq length($i2_s2);

    my $result = "OK";
    my $distance1 = calc_distance_between_two_seq($i1_s1, $i2_s1);
    my $distance2 = calc_distance_between_two_seq($i1_s2, $i2_s2);

    if ($distance1 < MIN_DISTANCE and $distance2 < MIN_DISTANCE){
        warn "[FAIL] distance too small for $tag ($distance1 and $distance2)\n";
        $result = "FAIL";
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
    say "Reading file $file";
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
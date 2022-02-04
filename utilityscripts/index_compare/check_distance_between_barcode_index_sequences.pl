#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

use Data::Dumper;
use Getopt::Long;
use File::Basename;
use 5.010.000;

use constant NA => "NA";

my $script_name = basename $0;
my $out_suffix = "index_comparison.tsv";
my $fail_str = 'FAIL';
my $success_str = 'OK';

my $out_name = undef;
my $max_conversion_mismatches = 1; # corresponds to setting in bcl2fastq conversion
my $help = undef;

my $HELP =<<HELP;

  Description:
    Checks all index combinations of any number of barcode index collections and determines the i7 and i5
    distance for every possible combination. Result for each is '$fail_str' in case total distance is
    lower than the allowed mismatches during conversion * 2 + 1 and otherwise '$success_str'.

  Usage:
    $script_name -out_name <output-name> <index-tsv-1> [<index-tsv-2> <index-tsv-n>]

  Example TSV input dual:
    #i7name i7seq    i5name i5seq
    BC1_i7  CTGATCGT BC1_i5 GCGCATAT

  Example output (written to <output-name>_$out_suffix):
    #result dist1 dist2 idx1_tag1 idx1_seq1 idx1_tag2 idx1_seq2 idx2_tag1 idx2_seq1 idx2_tag2 idx2_seq2
    OK 5 3 F498_LCH588 CTGCGGAT CZ-0379-LL_A AGATAACC IDT8_i7_384 CGACCATT IDT8_I5_384 TGATAGGC

  Notes:
    The type (single or dual) and length of all indexes should be identical.

HELP

GetOptions (
    "out_name=s" => \$out_name,
    "max_mismatches=i" => \$max_conversion_mismatches,
    "help|h" => \$help,
) or die "[ERROR] Issue in command line arguments\n";
my @index_input_files = @ARGV;
print $HELP and exit(0) if defined $help or scalar @ARGV == 0;

my $MIN_EDIT_DISTANCE = ($max_conversion_mismatches * 2) + 1;

die "[ERROR] Out name not defined\n" unless defined $out_name;
die "[ERROR] No index TSV files provided\n" unless scalar @index_input_files > 0;

say "[INFO] Starting with index comparison $out_name";
my $out_file = join("_", $out_name, $out_suffix);
my @all_indexes = ();
foreach my $file (@index_input_files) {
    die "[ERROR] File does not exist ($file)\n" unless -f $file;
    read_indexes_tsv($file, \@all_indexes);
}

my @header = qw(result dist1 dist2 idx1_tag1 idx1_seq1 idx1_tag2 idx1_seq2 idx2_tag1 idx2_seq1 idx2_tag2 idx2_seq2);
my @output_lines = ();

say "[INFO] Checking all possible combinations.";
check_all_index_combinations(\@all_indexes, \@output_lines);

open(my $OUT_FH, '>', $out_file) or die $!;
foreach my $file (@index_input_files) {
    say $OUT_FH '##' . " Input file: $file";
}
say $OUT_FH '##' . " Minimal distance for $success_str: $MIN_EDIT_DISTANCE";
say $OUT_FH '#' . join("\t", @header);
say $OUT_FH join("\n", @output_lines);
close $OUT_FH;

say "[INFO] Finished with $script_name (see $out_file for result).";

sub check_all_index_combinations{
    my ($indexes, $output_lines) = @_;
    my $index_count = scalar @$indexes;
    my $combination_count = 0;
    my $failure_count = 0;

    foreach my $x (1..$index_count){
        foreach my $y ($x..$index_count){
            next if $x == $y;
            $combination_count++;
            my $index1 = $indexes->[$x-1];
            my $index2 = $indexes->[$y-1];
            compare_two_indexes($index1, $index2, $output_lines, \$failure_count);
        }
    }

    if ($failure_count == 0) {
        say "[INFO] Total of $failure_count failing combinations encountered on a total of $combination_count combinations.";
    }
    else{
        warn "[WARN] Total of $failure_count failing combinations encountered on a total of $combination_count combinations.\n";
    }
}

sub compare_two_indexes{
    my ($index1, $index2, $output_array, $fail_count) = @_;

    my $i1_t1 = $index1->{tag1};
    my $i1_t2 = $index1->{tag2};
    my $i2_t1 = $index2->{tag1};
    my $i2_t2 = $index2->{tag2};

    my $i1_s1 = $index1->{seq1};
    my $i1_s2 = $index1->{seq2};
    my $i2_s1 = $index2->{seq1};
    my $i2_s2 = $index2->{seq2};

    die "[ERROR] Impossible to compare $i1_s1 with $i2_s1 due to different length\n" unless length($i1_s1) eq length($i2_s1);
    die "[ERROR] Impossible to compare $i1_s2 with $i2_s2 due to different length\n" unless length($i1_s2) eq length($i2_s2);

    # can be single or dual indexes so combine in case of dual
    my $i1_seq = NA;
    my $i2_seq = NA;
    my $tag_info = NA;
    my $seq_info = NA;
    my $is_dual_index = 0;

    if ($i1_t2 and $i2_t2) {
        $is_dual_index = 1;
        $i1_seq = join("+", $i1_s1, $i1_s2);
        $i2_seq = join("+", $i2_s1, $i2_s2);
        $tag_info = join(" ", "$i1_t1+$i1_t2", "vs", "$i2_t1+$i2_t2");
        $seq_info = join(" ", "$i1_s1+$i1_s2", "vs", "$i2_s1+$i2_s2");
    }else{
        $i1_seq = $i1_s1;
        $i2_seq = $i2_s1;
        $tag_info = join(" ", $i1_t1, "vs", $i2_t1 );
        $seq_info = join(" ", $i1_s1, "vs", $i2_s1 );
    }

    my $i1_distance = get_edit_distance_of_two_sequences($i1_s1, $i2_s1);
    my $i2_distance = get_edit_distance_of_two_sequences($i1_s2, $i2_s2);
    my $total_distance = get_edit_distance_of_two_sequences($i1_seq, $i2_seq);
    my $final_call = "";

    ## sanity check: total distance of both should always be sum of i1 and i2 distances
    if ($total_distance != $i1_distance + $i2_distance){
        die "[FAIL] Total is somehow not the sum of s1 and s2 distance ($total_distance != $i1_distance + $i2_distance) for $tag_info ($seq_info)\n";
    }

    ## determine fail or success (at least one index should meet the minimal requirement)
    if ($i1_distance < $MIN_EDIT_DISTANCE and $i2_distance < $MIN_EDIT_DISTANCE) {
        $final_call = $fail_str;
        $$fail_count++;
    }
    else{
        $final_call = $success_str;
    }

    ## output
    my @index_fields = join("\t", $i1_t1, $i1_s1, $i1_t2, $i1_s2, $i2_t1, $i2_s1, $i2_t2, $i2_s2);
    push( @$output_array, join("\t", $final_call, $i1_distance, $i2_distance, @index_fields));
}

sub get_edit_distance_of_two_sequences{
    my ($seq1, $seq2) = @_;
    my $distance = ( $seq1 ^ $seq2 ) =~ tr/\0//c;
    return $distance;
}

sub read_indexes_tsv{
    my ($file, $indexes) = @_;
    my $index_count = 0;
    say "[INFO] Reading file $file";
    open my $fh, '<', $file or die "Unable to open file ($file): $!\n";
    while ( <$fh> ) {
        chomp;
        next if $_ eq '' or $_ =~ /^#/;
        my ($tag1, $seq1, $tag2, $seq2) = split("\t", $_);

        if (not defined $tag1 or not defined $seq1){
            die "First two columns must contain tag1 and seq1";
        }
        elsif (not defined $tag2 or not defined $seq2){
            # create empty second in case of single index
            $tag2 = "";
            $seq2 = "";
        }

        # make sure all sequences are DNA
        die "[ERROR] Sequence seq1 contains chars other than ACGT ($seq1)?" unless $seq1 =~ /^[ACGT]*$/;
        die "[ERROR] Sequence seq2 contains chars other than ACGT ($seq2)?" unless $seq2 =~ /^[ACGT]*$/;

        # all ok so store result (in case of single index the tag2 and seq2 stay empty string)
        my %index = (tag1 => $tag1, tag2 => $tag2, seq1 => $seq1, seq2 => $seq2);
        push(@$indexes, \%index);
        $index_count++;
    }
    close $fh;
    say "[INFO]   Found $index_count indexes";
}
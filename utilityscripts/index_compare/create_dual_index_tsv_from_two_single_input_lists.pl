#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Text::CSV qw/csv/;
use 5.010.000;

use constant NA => "NA";

my $script_name = basename $0;
my $i7tsv;
my $i5tsv;
my $out_file;
my $help = undef;

my $HELP =<<HELP;

  Description:
    Reads in two single index lists with each two columns (name and sequence) and creates a dual index file.

  Usage:
    $script_name -out_name <output-name> -i7 <single-index-i7-tsv1> -i5 <single-index-i5-tsv2>
    $script_name -out_name ID_dual_indexes.tsv -i7 ID_i7_indexes.tsv -i5 ID_i5_indexes.tsv

  Example input i7 TSV:
    #i7name i7seq
    i7-bc-1  CTGATCGT
    i7-bc-2  CTGATCGA

  Example input i5 TSV:
    #i5name i5seq
    i5-bc-1  CTGATCGT
    i5-bc-2  CTGATCGA

  Example TSV output:
    #i7name i7seq    i5name i5seq
    i7-bc-1 CTGATCGT i5-bc-1 GCGCATAT
    i7-bc-2 CTGATCGT i5-bc-1 GCGCATAT
    i7-bc-1 CTGATCGT i5-bc-2 GCGCATAT
    i7-bc-2 CTGATCGT i5-bc-2 GCGCATAT

HELP

print $HELP and exit(0) if scalar @ARGV == 0;
GetOptions (
    "out_file=s" => \$out_file,
    "i5=s" => \$i5tsv,
    "i7=s" => \$i7tsv,
    "help|h" => \$help,
) or die "[ERROR] Issue in command line arguments\n";
print $HELP and exit(0) if defined $help;

defined $out_file || die "[ERROR] Out name not defined\n";
defined $i7tsv || die "[ERROR] i7tsv not defined\n";
defined $i5tsv || die "[ERROR] i7tsv not defined\n";

say "[INFO] Starting with $out_file";
my $i7idx = read_indexes_tsv($i7tsv, "#i7name", "i7seq");
my $i5idx = read_indexes_tsv($i5tsv, "#i5name", "i5seq");
my $combined = combine($i7idx, $i5idx, $out_file);

sub combine {
    my ($i7, $i5, $output_file) = @_;
    open my $fh, '>', $output_file or die "Unable to open output file ($output_file): $!\n";
    print $fh join("\t", '#i7name', 'i7seq', 'i5name', 'i5seq') . "\n";

    foreach my $i7seq (keys %$i7){
        my $i7name = $i7->{$i7seq};
        foreach my $i5seq (keys %$i5){
            my $i5name = $i5->{$i5seq};
            print $fh join("\t", $i7name, $i7seq, $i5name, $i5seq) . "\n";
        }
    }
    close $fh;
    say "[INFO] Written file $output_file"
}

sub read_indexes_tsv {
    my ($file, $name_header_tag, $seq_header_tag) = @_;
    my $index_count = 0;
    my %indexes = ();
    say "[INFO] Reading file $file";
    my $lines = csv(in => $file, headers => 'auto', sep_char=> "\t", skip_empty_rows => 1);
    foreach my $line (@$lines) {
        my $name = $line->{$name_header_tag} || dieAndPrintObject("Unable to get name from object", $line);
        my $seq = $line->{$seq_header_tag};
        $indexes{$seq} = $name;
        $index_count++;
    }
    say "[INFO]   Encountered $index_count indexes in file $file";
    return \%indexes;
}

sub dieAndPrintObject {
    my ($msg, $object) = @_;
    print Dumper $object;
    die "$msg\n";
}
#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
# use File::KGlob;

## ======================================
## Get options
## ======================================

my $help;
my $rundir;
my $debug;
my $javaMem = '8g';
my $slices = '/data/common/repos/scripts/dataslicer/CHPv2_OAv2_BRCA_nochr_merged.bed';
my $sambamba = '/data/common/tools/sambamba_v0.5.8/sambamba';
my $snpsift = '/data/common/tools/snpEff_v4.1h/SnpSift.jar';

die usage() if @ARGV == 0;
#get arguments
GetOptions(
    'h|help'	=> \$help,
    'slices=s'	=> \$slices,
    'rundir=s'	=> \$rundir,
    'debug'	=> \$debug,
);
die usage() if $help;
#checks
die "[ERROR] No run ID given.\n" unless $rundir;
die "Slices file $slices doesn't exist.\n" unless -e $slices;
die "Rundir $rundir doesn't exist.\n" unless -d $rundir;


## ======================================
## Slicing
## ======================================

my $run = basename($rundir);

## BAMS
my @bams = glob("$rundir/*/mapping/*dedup.realigned.bam");
foreach my $bam ( @bams ){
    my $outbam = $bam;
    $outbam =~ s/\.realigned.bam/\.realigned_sliced.bam/;
    if ( ! -e $outbam ){
	print "[INFO] Intersecting BAM files...\n";
	my $cmd = "$sambamba view -f bam -o $outbam -L $slices $bam";
	system( $cmd ) unless $debug;
    }
    my $outbai = "$outbam.bai";
    if ( ! -e $outbai ){
	print "[INFO] Indexing sliced BAM files...\n";
	my $cmd = "$sambamba index $outbam $outbai";
	system( $cmd ) unless $debug;
    }
}

## VCFs
my $annotVCF = (glob("$rundir/\*snpEff*.vcf"))[0]; #use annotated vcf
my $intersect_vcf = $annotVCF;
$intersect_vcf =~ s/\.vcf/\_sliced.vcf/;

if ( ! -e $intersect_vcf ){
    print "[INFO] Intersecting VCF file...\n";
    my $cmd = "java -Xmx$javaMem -jar $snpsift intervals $slices -i $annotVCF > $intersect_vcf";
    system( $cmd ) unless $debug;
}

## ======================================
## Subs
## ======================================


sub usage{
    print <<END;
    
  Description:
    This script slices BAM and VCF files based on a default or given BED file.
    The IAP directory structure is expected.

  Usage:
    perl HW_slicer.pl -rundir /path/to/rundir

END
  exit;
}

#!/usr/bin/perl

use strict;
use warnings;

use 5.16.0;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Find::Rule;
use Array::Diff qw(:all);
use JSON::XS;

my $NA_CHAR = '-';
my $SCRIPT = `basename $0`; chomp($SCRIPT);
my %GT_CONVERSION = ('./.'=>$NA_CHAR,'0/0'=>'A','1/0'=>'B','0/1'=>'B','1/1'=>'C','1/2'=>'D','2/1'=>'D');
my $MIN_USABLE = 18;
my $MIN_IDENTICAL = 18;
my $MAX_DISTANCE = 4;

my $snpcheck_db = "/data/snpcheck/output/";
my $delim_char = ' ';
my $input_vcf_dir = "";
my $output_file = "";
my $design_vcf = "";

my $use_tab_delim;
my $collapse_call;
my $qc_mode;
my $always_pass;

my $HELP_TEXT = <<EOF;
 ---
 Description:  
    Compares genotypes of two vcfs. Returns OK if SNPcheck criteria are met.
      criterium: minimal usable pos = $MIN_USABLE
      criterium: minimal identical calls = $MIN_IDENTICAL
      criterium: maximum sample distance = $MAX_DISTANCE
      setting: mismatch penalty HET to HOM = 1
      setting: mismatch penalty HOM to HOM = 2

    Run modus depends on input amount of input VCFs:
      * One input VCF: compare vs all vcfs in snpcheckDB ($snpcheck_db)
      * Two input VCFs: compare the two vcfs
      * Three or more input VCFs: no comparison just print table

 Usage:
    Compare one VCF file vs all in snpcheckDB
      $SCRIPT \${vcf1} -outputFile \$my_out_file
    Compare two VCF files and perform check:
      $SCRIPT \${vcf1} \${vcf2}
    Print table of calls for 3 or more VCFs:
      $SCRIPT /path/to/vcfs/*vcf
    Change output format with tab delim + single char calls
      $SCRIPT /path/to/vcfs/*vcf -tab -collapse
    Perform plate check
      $SCRIPT -qc -vcfdir /path/to/dir/with/vcfs/
   
 Options:
    -snpcheckDB    [s] path to dir with VCFs to serve as DB ($snpcheck_db)
    -outputFile    [s] path to output file with comparisons if run in genotypeDB mode
    -tabDelim          use tab as output delimiter for output table
    -collapseCall      collapse genotypes to single characters
    -designVcf     [s] vcf design file with positions to include
    -vcfDir        [s] directory with input vcf files
    -qc                quality check mode (needs vcfDir as well)
    -alwaysPass        sets QC values such that final result is always OK
 ---
EOF

my $help;
my $verbose;

print $HELP_TEXT and exit(0) if scalar(@ARGV) == 0;
GetOptions (
    "help|h" => \$help,
    "verbose" => \$verbose,
    "alwaysPass" => \$always_pass,
    "tabDelim" => \$use_tab_delim,
    "outputFile=s" => \$output_file,
    "collapseCall" => \$collapse_call,
    "snpcheckDB=s" => \$snpcheck_db,
    "designVcf=s" => \$design_vcf,
    "vcfDir=s" => \$input_vcf_dir,
    "qc" => \$qc_mode,
) or die "Error in command line arguments\n";
print $HELP_TEXT and exit(0) if $help;
$delim_char = "\t" if $use_tab_delim;

## input checks 
die "[ERROR] vcfDir is not a directory ($input_vcf_dir)\n" if $input_vcf_dir and not -d $input_vcf_dir;
die "[ERROR] Use either vcfDir or separate VCFs as input.\n" unless scalar @ARGV or $input_vcf_dir;
die "[ERROR] Use either vcfDir or separate VCFs as input, not both.\n" if scalar @ARGV and $input_vcf_dir;

## some cleanup
$input_vcf_dir =~ s/\/$//;

## set qc values such that we always pass if requested
if ($always_pass){
    say "[INFO] Running in Always pass mode!";
    $MIN_USABLE = 0;
    $MIN_IDENTICAL = 0;
    $MAX_DISTANCE = 1e6;
}

## -----
## main
## -----
my @vcfs = @ARGV;
my $vcfCount = scalar @vcfs;
my $posToKeep;
if ( $design_vcf ){
    $posToKeep = parseDesign( $design_vcf );
}

if ( $input_vcf_dir ){
   my @dir_vcfs = find->file()->name('*.vcf')->in($input_vcf_dir);
   @vcfs = ( @vcfs, @dir_vcfs );
}

my %names_seen = ();
my @names;
my %store;
foreach my $vcf ( @vcfs ){
    my $name = fileparse( $vcf );
    if ( exists $names_seen{ $name } ){
        die "[ERROR] Input name clash: VCF name appears more than once in input ($name)\n";
    }
    $names_seen{ $name } = "seen";
    push( @names, $name );
    parseVcf( $vcf, $name, \%store, $posToKeep );
}


if ( $qc_mode ){
    say "[INFO] Quality check mode: checking call count per VCF";
    die "[ERROR] Param vcfDir required with quality check mode.\n" unless $input_vcf_dir;
    my $plate_name = fileparse( $input_vcf_dir );
    say "[INFO] Working on plate: $plate_name";
    performQualityCheck( \%store, \@names );
}
elsif ( $vcfCount == 1 ){
    say "[INFO] Single vcf mode: will perform comparison to database ($snpcheck_db)";
    die "[ERROR] GenotypeDB is not a directory ($snpcheck_db)\n" if not -d $snpcheck_db;
    die "[ERROR] Param outputFile required with single vcf mode.\n" unless $output_file;
    die "[ERROR] Output file already exists ($output_file)\n" if -f $output_file;
    
    my $out = $output_file;
    my $log = $out . ".log";
    my $sample1 = $names[0];
    
    my @comparisons = ();
    foreach my $vcf2 (find->file()->name('*.vcf')->in($snpcheck_db)){
        my $sample2 = fileparse( $vcf2 );
        my %store1and2 = %store;
        parseVcf( $vcf2, $sample2, \%store1and2, $posToKeep );
        my $comparison = compareTwoSamples( \%store1and2, $sample1, $sample2, \*STDOUT );
        push( @comparisons, $comparison );
    }
    
    open my $out_fh, ">", $out or die "[ERROR] Unable to open file ($out): $!\n";
        printComparisonTable( \@comparisons, $out_fh );
    close $out_fh;
    
    say "[INFO] DONE: output in $out";
}
elsif ( $vcfCount == 2 ){
    say "[INFO] Two VCF mode: will perform pairwise comparison";
    my $comparison = compareTwoSamples( \%store, $names[0], $names[1], \*STDOUT );
    exit(0) if $comparison->{result} eq 'OK';
    exit(1);
}
else{
    say "[INFO] More than two input vcfs: will perform no comparisons";
    printSampleGenotypesTable( \%store, \@names);
}

## -----
## subroutines
## -----
sub parseDesign{
    my ( $file ) = @_;
    my %positionsToKeep = ();
    open IN, "<", $file or die "[ERROR] Unable to open $file: $!\n";
    while ( <IN> ){
        next if $_ =~ /^#/;
        my ( $chr, $pos ) = split( "\t", $_ );
        $positionsToKeep{$chr}{$pos} = 1;
    }
    close IN;
    return( \%positionsToKeep );
}

sub printComparisonTable{
    my ( $comparisons, $fh ) = @_;
    my @columns = qw( result usable failed skipped identical different sample1 sample2);
    my $matches = 0;
    
    say "[INFO] Setting up comparison table (any matches will also be printed to STDOUT)";
    say $fh '#' . join( "\t", @columns );
    foreach my $comparison ( @$comparisons ){
        say $fh join( "\t", map( $comparison->{$_}, @columns ) );
        if ( $comparison->{result} eq 'OK' ){
            say STDOUT "[INFO] Match between $comparison->{sample1} and $comparison->{sample2}" ;
            $matches++;
        }
    }
    say "[INFO] A total of $matches matches found";
}

sub performQualityCheck{
    my ( $storeObject, $vcfNames ) = @_;
    
    my @regions = sort keys %$storeObject;
    my $callField = 'gtCallsCollapsed';
    my @qc_pass = ();
    my @qc_fail = ();
    
    foreach my $vcfName ( @$vcfNames ){
        
        ## skip negative control sample
        next if $vcfName =~ /^NTC\d+(-NTC)?_[A-Z]{3}\d{2}_OpenArrayCalls\.vcf$/;
        
        my ($sampleName,$plateName) = (split( "_", $vcfName ))[0,1];
        my $usable = 0;
        my $total = 0;
        
        foreach my $region ( @regions ){
            my $call = $NA_CHAR;
            if ( exists $storeObject->{ $region }{ $vcfName }{ $callField } ){
                $total++;
                $call = $storeObject->{ $region }{ $vcfName }{ $callField };
            }            
            $usable++ unless $call =~ /$NA_CHAR/;
        }     
        my $is_qc_fail = $usable < $MIN_USABLE;
        
        if ( $is_qc_fail ){
            push( @qc_fail, "  [FAIL] $usable usable calls of $total total (Plate:$plateName Sample:$sampleName)");
        }
        else{
            push( @qc_pass, "  [PASS] $usable usable calls of $total total (Plate:$plateName Sample:$sampleName)");
        }
    }
    
    my $fail_count = scalar @qc_fail;
    say "[INFO] Sample QC Results:";
    say join( "\n", @qc_pass ) if scalar @qc_pass;
    say join( "\n", @qc_fail ) if scalar @qc_fail;
    if ( scalar @qc_fail ){
        say "[WARN] FINAL RESULT = FAIL: there were $fail_count samples failing qc!";
    }
    else{
        say "[INFO] FINAL RESULT = OK: all samples pass qc";
    }
}
    
sub printSampleGenotypesTable{
    my ( $storeObject, $vcfNames ) = @_;
    
    my @regions = sort keys %$storeObject;
    my @chroms = map( $storeObject->{ $_ }{ 'chr' }, @regions );
    say "[INFO] Genomic position order: " . join( ", ", @regions );
    
    foreach my $vcfName ( @$vcfNames ){
        my $callField = 'gtCalls';
        $callField = 'gtCallsCollapsed' if $collapse_call;
        my @calls = ();
        my $usable = 0;
        
        foreach my $region ( @regions ){
            my $call = $NA_CHAR;
            $call = $storeObject->{ $region }{ $vcfName }{ $callField } if exists $storeObject->{ $region }{ $vcfName }{ $callField };
            push( @calls, $call );
            $usable++ unless $call =~ /$NA_CHAR/;
        }
        say STDOUT join( "$delim_char", map( $_, @calls ), $usable, $vcfName );
    }
}

sub compareTwoSamples{
    my ( $storeObject, $sample1, $sample2 ) = @_;
    say "[INFO] Start comparison of $sample1 vs $sample2";
    
    my $usable    = 0;
    my $skipped   = 0;
    my $failed    = 0;
    my $identical = 0;
    my $different = 0;
    my $distance  = 0;
    
    my @print_order = sort keys %$storeObject;    
    my @toprint_chroms = ();
    my @toprint_calls1 = ();
    my @toprint_calls2 = ();
    my @toprint_result = ();
    my @toprint_distan = ();
    
    foreach my $loc ( @print_order ){
        
        ## get location object with both samples
        my $info = $storeObject->{ $loc };
        
        ## skip if not present in both samples
        unless ( exists $info->{ $sample1 } and exists $info->{ $sample2 } ){
            say "[NOTE] Ignoring location $loc: not present in both $sample1 and $sample2";
            next;
        }
        
        ## use the collapsed notation if requested
        my $callField = 'gtCalls';
        $callField = 'gtCallsCollapsed' if $collapse_call;
        
        my $chrom = $info->{ 'chr' };
        my $call1 = $info->{ $sample1 }{ $callField };
        my $call2 = $info->{ $sample2 }{ $callField };
        my $alleles1 = $info->{ $sample1 }{ 'alleles' };
        my $alleles2 = $info->{ $sample2 }{ 'alleles' };
        
        ## setup print output
        push( @toprint_chroms, $chrom );
        push( @toprint_calls1, $call1 );
        push( @toprint_calls2, $call2 );
        
        ## only use autosomal chromosomes
        if ( $chrom !~ m/^\d+$/ ){
            $skipped++;
            push( @toprint_result, $NA_CHAR );
            push( @toprint_distan, $NA_CHAR );
        }
        ## failed call
        elsif ( ($call1 =~ /$NA_CHAR/) or ($call2 =~ /$NA_CHAR/) ){
            $failed++;
            push( @toprint_result, $NA_CHAR );
            push( @toprint_distan, $NA_CHAR );
        }
        ## usable call
        else{
            ## identical call
            if ( $call1 eq $call2 ){
                $identical++;
                push( @toprint_result, 'Y' );
                push( @toprint_distan, 0 );
            }
            ## different call
            else{
                $different++;
                my $alleleDistance += calculateDistance( $alleles1, $alleles2 );
                $distance += $alleleDistance;
                push( @toprint_result, 'N' );
                push( @toprint_distan, $alleleDistance );
            }
        }
    }
    
    $usable = $identical + $different;
    my $identicalPerc = sprintf( "%d", 0 );
    my $differentPerc = sprintf( "%d", 0 );
    if ( $usable > 0 ){
        $identicalPerc = sprintf( "%d", ($identical*100 / $usable) );
        $differentPerc = sprintf( "%d", ($different*100 / $usable) );
    }
    
    ## append sample names
    push( @toprint_calls1, $sample1 );
    push( @toprint_calls2, $sample2 );
    
    ## print to stdout
    say "[INFO] Total skipped = $skipped";
    say "[INFO] Total failed = $failed";
    say "[INFO] Total usable = $usable (minimum: $MIN_USABLE)";
    say "[INFO] Total identical = $identical ($identicalPerc% of usable)";
    say "[INFO] Total different = $different ($differentPerc% of usable)";
    say "[INFO] Total distance  = $distance (maximum: $MAX_DISTANCE)";
    say join( $delim_char, "[INFO] Chrom: ", map( sprintf( "%3s", $_ ), @toprint_chroms ) );
    say join( $delim_char, "[INFO] file1: ", map( sprintf( "%3s", $_ ), @toprint_calls1 ) );
    say join( $delim_char, "[INFO] file2: ", map( sprintf( "%3s", $_ ), @toprint_calls2 ) );
    say join( $delim_char, "[INFO]  Dist: ", map( sprintf( "%3s", $_ ), @toprint_distan ) );
    say join( $delim_char, "[INFO] Match: ", map( sprintf( "%3s", $_ ), @toprint_result ) );
    
    ## determine final RESULT
    my $result;
    if ( $usable < $MIN_USABLE ){
        $result = "FAIL";
        say "[WARN] FINAL RESULT = $result (not enough usable pos: " . $usable . " < $MIN_USABLE)";
    }
    elsif( $identical < $MIN_IDENTICAL ){
        $result = "FAIL";
        say "[WARN] FINAL RESULT = $result (not enough identical calls: " . $identical . " < $MIN_IDENTICAL)";
    }
    elsif ( $distance > $MAX_DISTANCE ){
        $result = "FAIL";
        say "[WARN] FINAL RESULT = $result (distance between samples too big: " . $distance . " > $MAX_DISTANCE)";
    }
    else{
        $result = "OK";
        say "[INFO] FINAL RESULT = $result";
    }
    
    ## setup compare object (add 0 to enforce integer)
    my %comparison = ();
    $comparison{ sample1 } = $sample1;
    $comparison{ sample2 } = $sample2;
    $comparison{ skipped } = 0 + $skipped;
    $comparison{ failed } = 0 + $failed;
    $comparison{ usable } = 0 + $usable;
    $comparison{ identical } = 0 + $identical;
    $comparison{ different } = 0 + $different;
    $comparison{ distance } = 0 + $distance;
    $comparison{ result }  = $result;
    
    return( \%comparison );
}

sub calculateDistance{
    my ( $alleles1, $alleles2 ) = @_;
    my @sorted1 = sort @$alleles1;
    my @sorted2 = sort @$alleles2;
    my $diffobj = Array::Diff->diff( \@sorted1, \@sorted2 );
    my $dist = scalar( @{$diffobj->added} );
    
    return $dist;
}

sub parseVcf{
    
    my ( $file, $name, $storeObject, $posToKeep ) = @_;
    open IN, "<", $file or die "[ERROR] Unable to open file ($file): $!\n";
    
    ## read genotype calls and store them for later comparison
    while ( <IN> ){
        chomp;
        next if $_ =~ /^#/;
        my ( $chr, $pos, $id, $ref, $altString, $qual, $filter, $infoString, $formatString, $sampleString ) = split( "\t", $_ );

        ## if a list hash/dict with positions to keep was provided then filter
        next if $posToKeep and not exists $posToKeep->{ $chr }{ $pos };
        
        ## if no ref or alt info just skip position entirely
        next if $ref eq "";
        next if $altString eq "";
        
        ## gather all required info
        my @sampleInfoKeys = split( ':', $formatString );
        my @sampleInfoVals = split( ':', $sampleString );
        my @altVariants = split( ',', $altString );
        my @allAlleles = ( $ref, @altVariants );
        
        die "[ERROR] format ($formatString) and sample ($sampleString) columns diff field count...?" unless scalar @sampleInfoKeys == scalar @sampleInfoVals;
        
        my %sampleMetrics = ();
        foreach my $key ( @sampleInfoKeys ){
            my $val = shift @sampleInfoVals;
            $sampleMetrics{ $key } = $val;
        }
        
        my $genotypeCall = $sampleMetrics{ 'GT' };
        my @genotypes = split( '/', $genotypeCall );
        
        die "[ERROR] wrong number of genotyped alleles (should be 2) in $genotypeCall" unless scalar @genotypes == 2;
        
        my $genotype1 = $genotypes[0];
        my $genotype2 = $genotypes[1];
        my @calledAlleles = ($NA_CHAR, $NA_CHAR);
        
        die "[ERROR] Unknown genotype ($genotypeCall) at loc $chr:$pos" if not exists $GT_CONVERSION{ $genotypeCall };
        
        my $count1 = scalar @calledAlleles;
        
        my $collapsedGenotype = $GT_CONVERSION{ $genotypeCall };
        @calledAlleles = map( $allAlleles[ $_ ], @genotypes ) if $genotype1 =~ /^\d+/ ;        
        my $alleleCall = join( '/', @calledAlleles );
        
        if ( $verbose ){
            my $alleleString = join( '/', @calledAlleles );
            say "Found $alleleString on $chr:$pos (REF=$ref, ALT=$altString, GT=$genotypeCall)";
        }
        
        ## store info at location->sample
        my $locationString = $chr.':'.$pos;
        $storeObject->{ $locationString }{ 'chr' } = $chr;
        $storeObject->{ $locationString }{ $name }{ 'gtCalls' } = $alleleCall;
        $storeObject->{ $locationString }{ $name }{ 'gtCallsCollapsed' } = $collapsedGenotype;
        $storeObject->{ $locationString }{ $name }{ 'alleles' } = \@calledAlleles;
    }
    close IN;
    
}

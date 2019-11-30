#!/usr/bin/perl

use strict;
use warnings;
use Cwd qw(realpath);
use File::Basename qw(fileparse);
use Getopt::Long;
use Data::Dumper;
use HTML::TableExtract;
use List::Util qw(sum);
use IO::File;

my $DECIMALS = 1;

my $help;
my $searchDir ='./Data/Intensities/BaseCalls/Reports/';
my $outputName = 'conversionStats'; # ends up in output file names
my $hmfRunId = 'X00-0000'; # will be parsed out of report
my $sampleSheetFile = ''; # can be used to get extra sample info
my $outputToFiles;

## 
my $infoBlock = <<EOF;
-----
 Descr: Parses illumina bcl2fastq HTML output 
 Run with: 
    -d <ReportsDir>
 
 Options: 
    -o output base name [default: $outputName]
    -s path to SampleSheet.csv file
    -f also output to files in addition to screen
    -d path to reports dir [default: $searchDir]
 
 Example Output file names (only written with -f flag):
   HJKCJCCXX_conversionStats_perSample.tsv
   HJKCJCCXX_conversionStats_perLane.tsv

 Notes:
   1) if the internal run id cannot be found, the default is
      used in the output ($hmfRunId)
-----
EOF

printHelpMsg() unless scalar @ARGV;
GetOptions(
	"help|h" => \$help,
	"reportsDir|d=s" => \$searchDir,
	"outputName|o=s" => \$outputName,
	"hmfRunId|r=s" => \$hmfRunId,
	"sampleSheet|s=s" => \$sampleSheetFile,
	"outputToFiles|f" => \$outputToFiles,
) or die("Error in command line arguments\n");
printHelpMsg() if $help;
printHelpMsg( "Reports dir not present ($searchDir)" ) unless -d $searchDir;

## ----------
## read SampleSheet if provided
## ----------
my %sampleSheetInfo = ();
if ( $sampleSheetFile eq '' ){
	$sampleSheetFile = $searchDir;
	$sampleSheetFile =~ s/Data\/Intensities\/BaseCalls.*/SampleSheet\.csv/;
}

## get the full run name
my $shtPath = realpath $sampleSheetFile;
my ($runName, $runPath) = fileparse $shtPath;
chomp( $shtPath );
chomp( $runPath );
chomp( $runName );

if ( $sampleSheetFile ne '' ){
	open FILE, $sampleSheetFile or die "Couldn't open file ($sampleSheetFile): $!";
	while ( my $line = <FILE> ) {
		chomp($line);
		## any whitespace not allowed in csv
		$line =~ s/\s+//g;
		next if $line =~ /^,/;
		next if $line eq "";
		
		my @fields = split( ",", $line);
		
		#print "[DEBUG] FIELDS:\n";
		#print Dumper( \@fields )."\n";
		
		my ($key, $value) = @fields[0,1];
		#if ($key eq 'ExperimentName' and $value ne ''){
		if ($key eq 'ExperimentName' and $value ){
			$hmfRunId = $value;
			if ( $hmfRunId =~ /(X\d{2})(\d{4})/ ){
				$hmfRunId = $1.'-'.$2; # adding dash in case lab forgot
			}
		}
		if ($key =~ m/^[12345678]{1}$/ ){ # only read lines starting with a lane nr
			my ($id, $name, $indexseq, $project, $desc) = @fields[1,2,6,7,8];
			$sampleSheetInfo{ $name }{ sampleid } = $id;
			$sampleSheetInfo{ $name }{ indexseq } = $indexseq;
			$sampleSheetInfo{ $name }{ project } = $project;
			$sampleSheetInfo{ $name }{ desc } = $desc || 'NA';
			$sampleSheetInfo{ $name }{ fastq } = 'noFastq';
			$sampleSheetInfo{ $name }{ fastqDir } = "$runPath/Data/Intensities/BaseCalls/$project/$id/";
			
			## depending on how the sample-id and sample-name are set in SampleSheet.csv
			## the resulting fastq files can be at different dir level (2 or 3 deep)
			if ( $name eq $id ){
				$sampleSheetInfo{ $name }{ fastqDir } = "$runPath/Data/Intensities/BaseCalls/$project/";
			}
			my $checkFqExistenceCmd = "find $runPath/Data/Intensities/BaseCalls/ -mindepth 2 -maxdepth 3 -type f -iname \"*$name*fastq.gz\"";
			my $sampleHasFastqFiles = `$checkFqExistenceCmd`;
			if ( $sampleHasFastqFiles ){
				$sampleSheetInfo{ $name }{ fastq } = 'hasFastq';
			}
			
		}
	}	
	close FILE;
}

#print "[DEBUG] SampleSheetInfo\n";
#print Dumper( \%sampleSheetInfo );

my $sampleCount = scalar keys( %sampleSheetInfo );
print "## [INFO] Number of sample found in sampleSheet: ".$sampleCount."\n";

## ----------
## read laneBarcode.html
## ----------
my $find_cmd = "find $searchDir/ -maxdepth 6 -mindepth 6 -wholename '*all/all/all/laneBarcode.html'";
my $html_file = `$find_cmd`;
chomp( $html_file ); # remove newline
die "[ERROR] Required html file not found ($html_file)\n" unless (-f $html_file);
open FILE, $html_file or die "Couldn't open file ($html_file): $!"; 
	my $raw_html = join("", <FILE>); 
close FILE;

## ----------
## try to get sequence-run name from path
## ----------
my $seqrunName = 'NA';
my $seqrunDate = 'NA';
my $machine = 'NA';
my $fcid = 'NA';
## eg: 160115_ST-E00278_0021_BHJKHGCCXX
if ( $html_file =~ m/.*\/((\d{6})_(.*)_\d{4}_[AB]{1}(.{9}))\/.*/i ){
	$seqrunName = $1;
	$seqrunDate = $2;
	$machine = $3;
        $fcid = $4;
	print( "## [INFO] seqRun date found in html code: $seqrunDate ($machine)\n")
}


## ----------
## grep the flowcell name
## ----------
#my $fcid = 'NA';
#if ( $raw_html =~ m/(.+)\/all\/all/i ){
#	$fcid = $1;
#	print( "## [INFO] flowcell ID found in html code: $fcid\n")
#}else{
#	die "[ERROR] EXIT: unable to find flowcell id in html file ($html_file)\n";
#}

## ----------
## output of an earlier version of bcl2fq converter is not supported
## ----------
if ( $raw_html =~ m/<th colspan="4">Lane<\/th>/ ){
	die "[ERROR] EXIT: old bcl2fq converter output not supported (FCID: $fcid)\n";
}

## ----------
## Define all in/out column/field names
## ----------
my @fieldsRgxFlow = (
	'Clusters\s+\(Raw\)', 'Clusters\(PF\)', 'Yield'
);
my @fieldsFlow = qw( clustersRaw clustersPF yield );
my @fieldsRgxLane = ( 
	'Lane', 'Project', 'Sample', 
	'Barcode\s+sequence', 'PF\s+Clusters', 
	'Perfect\s+barcode', 'One\s+mismatch\s+barcode', 
	'Yield', 'PF', 'Q30', 'Mean\s+Quality',
);
my @fields = qw( lane project sample bc clusters mBC mmBC yield passf q30 qual );
my @outputFields = qw( yield q30 mBC mmBC qual );
my $typeSkipRegex = 'sample|bc|lane|project';
die "[ERROR] Unequal arrays (fields vs regexes)\n" unless scalar(@fields) == scalar(@fieldsRgxLane);

my $te_flow = HTML::TableExtract->new( headers => \@fieldsRgxFlow );
my $te_lane = HTML::TableExtract->new( headers => \@fieldsRgxLane );

$te_flow->parse( $raw_html );
$te_lane->parse( $raw_html );

## Examine all matching tables
my %stats = ();
my %samplesFound = ();
my %lanesFound = ();

## should have matched one table
my ($flow_table) = $te_flow->tables;
my ($lane_table) = $te_lane->tables;

## read the flowcell table
#print Dumper( $flow_table->rows );
my $rows = $flow_table->rows;
my ( $flowcellClustersRaw, $flowcellClustersPF, $flowcellYield ) = @{$rows->[0]};
$flowcellClustersRaw =~ s/,//g;
$flowcellClustersPF =~ s/,//g;
$flowcellYield =~ s/,//g;
my $flowcellPercClustersPF = sprintf("%.1f", ($flowcellClustersPF * 100 / $flowcellClustersRaw) );

## read the perLane table
foreach my $row ( $lane_table->rows ) {
	my %vals = ();	
	$vals{ $_ } = shift @$row foreach @fields;

	## init: make sucd re all present
	foreach ( @fields ){
		unless( $vals{ $_ }	){
			$vals{ $_ } = 0		
		}
	}
	
	## allow no comma's anywhere
	map( $vals{ $_ } =~ s/,//g, keys %vals);

	my $sample = $vals{ sample };
	my $lane = $vals{ lane };

	next if $sample eq "Undetermined";

	$samplesFound{ $sample } = 1;
	$lanesFound{ $lane } = 1;

	foreach my $f (@fields){
		my $v = $vals{ $f };
		$v = 0 if ($v =~ /^Nan|nan$/);
		push( @{ $stats{ $f }{ perSample }{ $sample }}, $v );
		push( @{ $stats{ $f }{ perLane }{ $lane }}, $v );	
		push( @{ $stats{ $f }{ totals }}, $v );	
	}	
}



## --------------------
## Setup output-data per Sample and Lane
## --------------------
my $samplesYield = getAggregateValue( 'yield', $stats{ 'yield' }{ 'totals' } );
my $samplesQ30 = getAggregateValue( 'q30', $stats{ 'q30' }{ 'totals' } );
my $undetYield = $flowcellYield - $samplesYield;
print "## [INFO] Flowcell Total Yield = $flowcellYield\n";
#print "## [INFO] Samples Total Yield = $samplesYield\n";
print "## [INFO] Undetermined Yield = $undetYield\n";
print "## [INFO] Samples Av Q30 = $samplesQ30\n";
print "## [INFO] Passing Filter Clusters % = $flowcellPercClustersPF\n";
my @runOverviewVals = ( $flowcellYield, $undetYield, $samplesQ30, $flowcellPercClustersPF );
print "## [INFO] Next line can be pasted into runOverview sheet:\n";
print "## [INFO] ".join( "\t", @runOverviewVals )."\n";

my %meanStatsPerSample = ();
my %meanStatsPerLane = ();
foreach my $type ( keys %stats ) {
	next if ( $type =~ /$typeSkipRegex/ ); ## these wont end up in output
			
	foreach my $sample ( keys %samplesFound ) {
		my @data = @{$stats{ $type }{ perSample }{ $sample }};
		my $aggr = getAggregateValue( $type, \@data );
		if ( $aggr =~ /NaN|nan/ ){
			warn "[WARN] Sample $sample for type $type is NaN\n";	
			$aggr = 'NA';
		}
		$meanStatsPerSample{ $sample }{ $type } = $aggr;
		#print "TYPE; $type\n";
		#print Dumper( \@data );
		#<>;
	}

	foreach my $lane ( keys %lanesFound ) {
		my @data = @{$stats{ $type }{ perLane }{ $lane }};
		my $aggr = getAggregateValue( $type, \@data );
		if ( $aggr =~ /NaN|nan/ ){
			warn "[WARN] Lane $lane for type $type is NaN\n";
			$aggr = 'NA';
		}
		$meanStatsPerLane{ $lane }{ $type } = $aggr;
	}
}
#print Dumper( \%statsPerSample );

## --------------------
## OUTPUT
## --------------------
my @outputHandlesPerSample = ( *STDOUT );
my @outputHandlesPerLane = ( *STDOUT );

if ( $outputToFiles ){
	my $sampleOutFile = $hmfRunId.'_'.$fcid.'_'.$outputName.'_perSample.tsv';
	my $laneOutFile = $hmfRunId.'_'.$fcid.'_'.$outputName.'_perLane.tsv';
	open PER_S, ">$sampleOutFile" or die "Unable to open out file ($sampleOutFile)\n";
	open PER_L, ">$laneOutFile" or die "Unable to open out file ($laneOutFile)\n";
	push( @outputHandlesPerSample, *PER_S );
	push( @outputHandlesPerLane, *PER_L );
	
	foreach my $fh ( *PER_S, *PER_L ){
		print $fh '## yield=yield in Mbases'."\n";
		print $fh '## mBC=matched barcodes without mismatch'."\n";
		print $fh '## mmBC=matched barcodes with mismatch'."\n";
		print $fh '## q30=percentage bases with qual above 30'."\n";
		print $fh '## qual=mean base quality'."\n";
		print $fh '## passf=percentage clusters passing filters'."\n";
	}
}

my @sampleSheetFields = qw( sampleid indexseq project fastq desc fastqDir );
my @outputFields2 = qw( date sequencer hmfrunid flowcell );

#print "[DEBUG] Samplesheet Hash:\n";
#print Dumper( \%sampleSheetInfo );

foreach my $fh ( @outputHandlesPerSample ){
	print $fh '## PerSampleInfo'."\n";
	print $fh '#'.join( "\t", @outputFields, @outputFields2, 'samplename', @sampleSheetFields )."\n";
	foreach my $sample ( sort keys %meanStatsPerSample ) {
		my @data = ();
		push( @data, ($meanStatsPerSample{ $sample }{ $_ } || 'NA' )) foreach @outputFields;
		push( @data, $seqrunDate, $machine, $hmfRunId, $fcid, $sample);
		#push( @data, ($sampleSheetInfo{ $sample }{ $_ } || 'NA' )) foreach @sampleSheetFields;
		push( @data, $sampleSheetInfo{ $sample }{ $_ } || 'NA' ) foreach @sampleSheetFields;
		#print Dumper( \@data )."\n";
		print $fh join( "\t", @data )."\n";
	}
}

foreach my $fh ( @outputHandlesPerLane ){
	print $fh '## PerLaneInfo'."\n";
	print $fh '#'.join( "\t", @outputFields, @outputFields2, 'lane' )."\n";
	foreach my $lane ( sort keys %meanStatsPerLane ) {
		my @data = ();
		push( @data, ($meanStatsPerLane{ $lane }{ $_ } || 'NA' )) foreach @outputFields;
		print $fh join( "\t", ( @data, $seqrunDate, $machine, $hmfRunId, $fcid, $lane) )."\n";
	}
}

if ( $outputToFiles ){
	close PER_S;
	close PER_L;
}

## ====================
## SUBROUTINES
## ====================
sub getAggregateValue{
	my ($type, $arrayRef) = @_;
	my $out;
	my $sum = sum( @$arrayRef );

	if( $type =~ /yield|clusters/ ){ ## sum of array
		$out = $sum;

	}elsif( $type =~ /mBC|q30|qual|passf/ ){ ## average of array
		$out = sprintf( '%.'.$DECIMALS.'f', $sum / @$arrayRef);
	}else{ ## undefined field?
		die "[ERROR] something wrong in getAggregateValue\n";
	}

	## Undetermined is filtered out, this could correct if still present
	#$out = 0 if ($out =~ /^NaN|nan$/);
	return( $out );
}

sub printHelpMsg{
	my @msg = @_;
	print $infoBlock;
	print "[ERROR] $_\n" foreach @msg;
	exit(1);
}

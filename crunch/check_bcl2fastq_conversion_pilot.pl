#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use File::Basename;
use Getopt::Long;
use File::Slurp;
use JSON;
use XML::Simple;
use List::Util qw/sum/;
use 5.010.000;

my $SCRIPT = basename $0;
my $OUT_SEP = "\t";
my $NA_CHAR = "NA";

my @OUT_FIELDS = qw(flowcell yld q30 pf yld_p mm0 mm1 id name submission index1 index2);
my @SUM_FIELDS = qw(submission id name yld q30 index1 index2);
my $ROUND_DECIMALS = 1;

my $RUN_PATH;
my $OUT_JSON_PATH;

my $MCSV_PATH = 'Fastq/Reports/Quality_Metrics.csv2';
my $JSON_PATH = 'Fastq/Stats/Stats.json';
my $RXML_PATH = 'RunInfo.xml';
my $SSHT_PATH = 'SampleSheet.csv';

my %SETTINGS_PER_PLATFORM = (
    'NovaSeq' => {'min_flowcell_q30' => 85, 'min_sample_yield' => 1e9, 'max_undetermined' =>  8, 'yield_factor' => 1e6},
    'NextSeq' => {'min_flowcell_q30' => 75, 'min_sample_yield' => 1e9, 'max_undetermined' => 50, 'yield_factor' => 1e6},
    'ISeq'    => {'min_flowcell_q30' => 75, 'min_sample_yield' => 1e6, 'max_undetermined' => 50, 'yield_factor' => 1},
    'HiSeq'   => {'min_flowcell_q30' => 75, 'min_sample_yield' => 1e9, 'max_undetermined' =>  8, 'yield_factor' => 1e6},
);
my @KNOWN_PLATFORMS = keys %SETTINGS_PER_PLATFORM;

my $HELP =<<HELP;

  Description
    Parses conversion output json and prints table with flowcell/lane/sample/read metrics
    
  Usage
    $SCRIPT -runDir \${run-path}
    - OR -
    $SCRIPT -sampleSheet /path/to/SampleSheet.csv -statsJson /path/to/Stats.json -metricsCsv /path/to/Quality_Metrics.csv
    
  Options
    -sep <s>           Output sep (default = <TAB>)
    -yieldFactor <i>   Factor to divide all yields with (default present per platform)
    -decimals    <i>   Decimals to keep for q30 (default: $ROUND_DECIMALS)
    -sampleSheet <s>   Path to SampleSheet.csv file
    -statsJson   <s>   Path to Stats.json file
    -runInfoXml  <s>   Path to RunInfo.xml file
    -metricsCsv  <s>   Path to Quality_Metrics.csv file
    -outputJson  <s>   Path to output json file to create (in addition to stdout output)
    -noQc              Skip QC checks and just print table
    -summary           Prints extra sample summary table
    -debug             Prints complete datastructure
    
HELP
print $HELP and exit(0) if scalar @ARGV == 0;

## -----
## Gather input
## -----
my %opt = ();
GetOptions (
  "runDir=s"      => \$RUN_PATH,
  "statsJson=s"   => \$JSON_PATH,
  "metricsCsv=s"  => \$MCSV_PATH,
  "sampleSheet=s" => \$SSHT_PATH,
  "runInfoXml=s"  => \$RXML_PATH,
  "outputJson=s"  => \$OUT_JSON_PATH,
  "sep=s"         => \$OUT_SEP,
  "decimals=i"    => \$ROUND_DECIMALS,
  "yieldFactor=i" => \$opt{ yield_factor },
  "noQc"          => \$opt{ no_qc },
  "summary"       => \$opt{ print_summary },
  "debug"         => \$opt{ debug },
  "help|h"        => \$opt{ help },
) or die "[ERROR] Issue in command line arguments\n";
print $HELP and exit(0) if $opt{ help };

my $seq_run = "FileModeSoNoRunDir";
if (defined $RUN_PATH){
    die "[ERROR] Provided run dir does not exist ($RUN_PATH)\n" if not -d $RUN_PATH;
    $RXML_PATH = "$RUN_PATH/$RXML_PATH";
    $SSHT_PATH = "$RUN_PATH/$SSHT_PATH";
    $JSON_PATH = "$RUN_PATH/$JSON_PATH";
    $MCSV_PATH = "$RUN_PATH/$MCSV_PATH";
    $seq_run = basename $RUN_PATH;
}
die "[ERROR] Provided xml does not exist ($RXML_PATH)\n" if not -f $RXML_PATH;
die "[ERROR] Provided sample sheet does not exist ($SSHT_PATH)\n" if not -f $SSHT_PATH;

my $ssht_info = readSampleSheet( $SSHT_PATH );
my $platform = determinePlatformByString($ssht_info->{'platform'}, \@KNOWN_PLATFORMS);
my $SETTINGS = $SETTINGS_PER_PLATFORM{$platform};
$SETTINGS->{yield_factor} = $opt{yield_factor} if defined $opt{yield_factor};
my $YIELD_FACTOR = $SETTINGS->{yield_factor};

my %store = ();

say "## Input File: $RXML_PATH";
my $runinfo = parseRunInfoXml($RXML_PATH);

if( -f $MCSV_PATH ){
    say "## Input File: $MCSV_PATH";
    my $metrics_info = readQualityMetricsCsv($MCSV_PATH);
    parseMetricsInfo(\%store, $metrics_info, $runinfo);
}
elsif ( -f $JSON_PATH ){
    say "## Input file: $JSON_PATH";
    my $json_info = readJson( $JSON_PATH );
    parseJsonInfo(\%store, $json_info, $runinfo);
}
else{
    die "[ERROR] Both Stats JSON and Metrics CSV not found. Need one of both!\n"
}

$store{stats}{platform} = $platform;

printJson(\%store, "test.json");

addSamplesheetInfo(\%store, $ssht_info, $seq_run);
performQC(\%store, $SETTINGS) unless $opt{no_qc};

if ( $opt{print_summary} ){
    printSummaryTable(\%store, \@SUM_FIELDS );
}
else{
    printTable(\%store, \@OUT_FIELDS);
    printJson(\%store, $OUT_JSON_PATH) if defined $OUT_JSON_PATH;
}

if ( $opt{ debug } ){
    say "[DEBUG] Samplesheet data structure";
    print Dumper $ssht_info;
    say "[DEBUG] Complete final data structure";
    print Dumper \%store;
}

## -----
## /MAIN
## -----

sub parseMetricsInfo{
    my ($store, $metrics, $runinfo) = @_;

    my $total_non_index_cycle_count = $runinfo->{total_non_index_cycle_count};
    my $cycle_string = $runinfo->{cycle_string};
    my $fid = $runinfo->{flowcell_id};

    $store->{ 'flow' }{ $fid }{ 'id' } = $fid;
    $store->{ 'flow' }{ $fid }{ 'name' } = "TODO";

    while (my ($lid, $lane) = each %$metrics) {
        $store->{lane}{ $lid }{name} = $lid;
        $store->{flow}{ $fid }{clust_raw} += 0;
        $store->{flow}{ $fid }{clust_pf} += 0;
        $store->{lane}{ $lid }{clust_raw} += 0;
        $store->{lane}{ $lid }{clust_pf} += 0;

        while (my ($sid, $sample) = each %$lane) {
            $store->{samp}{ $sid }{name} = "NoName";
            $store->{samp}{ $sid }{index1}  = "NoIndexSeq1";
            $store->{samp}{ $sid }{index2}  = "NoIndexSeq2";

            while (my ($rid, $read) = each %$sample) {
                my $yield = $read->{Yield};
                my $yield_q30 = $read->{YieldQ30};
                my $index_seq1 = $read->{index};
                my $index_seq2 = $read->{index2};
                my $seq = $indexseq1 . '+' . $$indexseq2;
                $store->{samp}{$sid}{index1} = $index_seq1;
                $store->{samp}{$sid}{index2} = $index_seq2;
                $store->{flow}{$fid}{yield} += $yield;
                $store->{lane}{$lid}{yield} += $yield;
                $store->{samp}{$sid}{yield} += $yield;
                $store->{read}{$rid}{yield} += $yield;
                $store->{indx}{$seq}{yield} += $yield;
                $store->{flow}{$fid}{yield_q30} += $yield_q30;
                $store->{lane}{$lid}{yield_q30} += $yield_q30;
                $store->{samp}{$sid}{yield_q30} += $yield_q30;
                $store->{read}{$rid}{yield_q30} += $yield_q30;
                $store->{indx}{$seq}{yield_q30} += $yield_q30;
            }
        }
    }
    addPrintFields(\%store, $fid);
    addGeneralStats(\%store, $fid, $cycle_string);

}

sub addSamplesheetInfo{
    my ($data, $ssht_info, $seq_run) = @_;

    ## add info to stats (eg hmf_runname could be X17-0001)
    my $hmf_runname = $ssht_info->{'runname'};
    $data->{ 'stats' }{ 'seq_runname' } = $seq_run;
    $data->{ 'stats' }{ 'hmf_runname' } = $hmf_runname;
    $data->{ 'stats' }{ 'submissions' } = {};
    
    ## override flowcell name with ExperimentName from SampleSheet
    my @fcids = keys %{$data->{'flow'}};
    my $flowcell_count = scalar @fcids;
    die "[ERROR] There should be exactly one flowcell in data but found $flowcell_count\n" unless $flowcell_count == 1;
    my $fcid = $fcids[0];
    $data->{ 'flow' }{ $fcid }{ 'name' } = $hmf_runname;
    $data->{ 'flow' }{ $fcid }{ 'name_print' } = $hmf_runname;
    
    ## add sample metadata
    my $samples = $data->{ 'samp' };
    foreach my $sample_id ( keys %$samples ){
        my $sample = $samples->{ $sample_id };
        
        my $submission = $NA_CHAR;
        $submission = $ssht_info->{'samples'}{$sample_id}{ 'Sample_Project' } if defined $ssht_info->{'samples'}{$sample_id}{ 'Sample_Project' };
        $sample->{ 'submission_print' } = $submission;
        $data->{ 'stats' }{ 'submissions' }{ $submission } = 1;
        
        my $description = $NA_CHAR;
        $description = $ssht_info->{'samples'}{$sample_id}{ 'Description' } if defined $ssht_info->{'samples'}{$sample_id}{ 'Description' };
        $sample->{ 'description_print' } = $description;
    }
}

sub parseRunInfoXml{
    my ($run_info_xml) = @_;
    my $content = readXml( $run_info_xml );
    my $flowcell_id = $content->{Run}{Flowcell};
    my @cycle_counts = map( $_->{NumCycles}, @{$content->{Run}{Reads}{Read}});
    my $cycle_string = join( "|", @cycle_counts);
    my @read_cycle_counts = map( $_->{IsIndexedRead} eq "N" ? $_->{NumCycles} : (), @{$content->{Run}{Reads}{Read}});
    my $total_non_index_cycle_count = sum(@read_cycle_counts);
    my %info = (
        "cycle_string" => $cycle_string,
        "total_non_index_cycle_count" => $total_non_index_cycle_count,
        "flowcell_id" => $flowcell_id
    );
    return \%info;
}

sub readXml{
    my ($file) = @_;
    my $obj = XMLin($file);
    return($obj);
}

sub readJson{
    my ($json_file) = @_;
    my $json_txt = read_file( $json_file );
    my $json_obj = decode_json( $json_txt );
    return( $json_obj );
}

sub determinePlatformByString {
    my ($platform, $known_platforms) = @_;
    my $final_platform = "";

    ## exact platform string varies too much so try to match by regex
    foreach my $known (@$known_platforms) {
        if ($platform =~ m/^$known/i) {
            $final_platform = $known;
            say "## INFO: platform configured to $known (based on input '$platform')";
        }
    }

    if ($final_platform eq "") {
        die "[ERROR] Unable to determine platform from string ($platform)\n";
    }
    return $final_platform;
}

sub performQC{
    my ($info, $qc_limits) = @_;
   
    my $stats = $info->{'stats'};
    my $samps = $info->{'samp'};
    my $lanes = $info->{'lane'};
    my $fails = 0;
    my $identifier = $stats->{'identifier'};

    ## flowcell checks
    my $undet = $stats->{'undet_perc'};
    my $max_undet = $qc_limits->{ 'max_undetermined' };
    if ( $undet > $max_undet ){
        warn "## WARNING Percentage undetermined ($undet) too high (max=$max_undet)\n";
        $fails += 1;
    }
        
    ## lane and sample checks
    $fails += checkObjectField( $lanes, 'q30',   $qc_limits->{ 'min_flowcell_q30' } );
    $fails += checkObjectField( $samps, 'yield', $qc_limits->{ 'min_sample_yield' } );
    
    ## conclusion
    my $final_result = "NoQcResult";
    if ( $fails == 0 ){
        $final_result = "PASS";
        say "## FINAL QC RESULT: OK";
    }
    else{
        $final_result = "FAIL";
        warn "## WARNING Some checks failed, inspect before proceeding (for $identifier)\n";
        say "## FINAL QC RESULT: FAIL ($fails failures for $identifier)";
    }
    $stats->{'flowcell_qc'} = $final_result;
}

sub checkObjectField{
    my ($objects, $field, $min) = @_;
    my $fails = 0;
    foreach my $obj_key ( sort { $objects->{$b}{'name'} cmp $objects->{$a}{'name'} } keys %$objects){
        my $obj = $objects->{$obj_key};
        my $name = $obj->{'name'};
        next if $name eq 'UNDETERMINED';
        next if $name =~ /^VirtualSample/;
        my $value = 0;
        $value = $obj->{$field} if exists $obj->{$field};
        if ( $value < $min ){
            warn "## WARNING $field for $name too low: $value < $min\n";
            $fails += 1;
        }
    }
    return $fails;
}

sub parseJsonInfo{
    my ($store, $stats, $runinfo) = @_;

    my $total_non_index_cycle_count = $runinfo->{total_non_index_cycle_count};
    my $cycle_string = $runinfo->{cycle_string};

    my $fid = $stats->{ 'Flowcell' };
    $store->{ 'flow' }{ $fid }{ 'id' } = $fid;
    $store->{ 'flow' }{ $fid }{ 'name' } = $stats->{ 'RunId' };

    ## First reading the "unknown barcodes" and add them to "index sequences"
    my $unknowns = $stats->{'UnknownBarcodes'};
    foreach my $lane ( @$unknowns ){
        my $lid = join( "", "lane", $lane->{ Lane } );
        my $unknown_barcodes = $lane->{'Barcodes'};
        foreach my $bc ( keys %$unknown_barcodes ){
            $store->{ indx }{ $bc }{ name } = 'IndexFromUnknown';
            my $seq1 = (split(/\+/, $bc))[0] || $NA_CHAR;
            my $seq2 = (split(/\+/, $bc))[1] || $NA_CHAR;
            $store->{ indx }{ $bc }{ index1 } = $seq1;
            $store->{ indx }{ $bc }{ index2 } = $seq2;
            # Unlike actual samples, the unknowns are reported as cluster counts instead of yield
            # So need to calculate the yield using non-index cycle counts from RunInfo.xml
            $store->{ indx }{ $bc }{ yield } += $unknown_barcodes->{ $bc } * $total_non_index_cycle_count;
        }
    }

    ## Then read samples
    my $lanes = $stats->{'ConversionResults'};
    foreach my $lane ( @$lanes ){
        my $lid = join( "", "lane", $lane->{LaneNumber} );
        $store->{lane}{ $lid }{name} = $lid;

        $store->{flow}{ $fid }{clust_raw} += $lane->{TotalClustersRaw};
        $store->{flow}{ $fid }{clust_pf} += $lane->{TotalClustersPF};
        $store->{lane}{ $lid }{clust_raw} += $lane->{TotalClustersRaw};
        $store->{lane}{ $lid }{clust_pf} += $lane->{TotalClustersPF};
        
        ## Undetermined info is stored separate from samples in json
        my $undet_id = 'UNDETERMINED';
        my $undet_obj = $lane->{Undetermined};
        my $undet_reads = $undet_obj->{ReadMetrics};
        my $undet_info = \%{$store->{undt}{ $undet_id }};
        $undet_info->{name} = $undet_id;
        foreach my $read ( @$undet_reads ){
            my $rid = join( "", "read", $read->{ReadNumber} );
            $undet_info->{yield} += $read->{Yield};
            $undet_info->{yield_q30} += $read->{YieldQ30};
            $store->{flow}{ $fid }{yield} += $read->{Yield};
            $store->{flow}{ $fid }{yield_q30} += $read->{YieldQ30};
            $store->{lane}{ $lid }{yield} += $read->{Yield};
            $store->{lane}{ $lid }{yield_q30} += $read->{YieldQ30};
            $store->{read}{ $rid }{yield} += $read->{Yield};
            $store->{read}{ $rid }{yield_q30} += $read->{YieldQ30};
        }

        my $samples = $lane->{DemuxResults};
        foreach my $sample ( @$samples ){
            
            ## Sanity checks
            die "Field for sample id not found\n" unless defined $sample->{SampleId};
            die "Field for sample name not found\n" unless defined $sample->{SampleName};

            my $sid = $sample->{SampleId};
            my $snm = $sample->{SampleName};
            my $seq = $sample->{IndexMetrics}[0]{IndexSequence} || $NA_CHAR;

            my $seq1 = (split(/\+/, $seq))[0] || $NA_CHAR; # in case of one sample there is no first index
            my $seq2 = (split(/\+/, $seq))[1] || $NA_CHAR; # in case of single index there is no second index

            ## Reset info for all real samples
            $store->{samp}{ $sid }{name} = $snm;
            $store->{samp}{ $sid }{index1}  = $seq1;
            $store->{samp}{ $sid }{index2}  = $seq2;
            
            my $reads = $sample->{ReadMetrics};
            foreach my $read ( @$reads ){
                my $rid = join( "", "read", $read->{ReadNumber} );
                $store->{read}{ $rid }{name} = $rid;
                $store->{indx}{ $seq }{name} = 'IndexFromSample';
                $store->{indx}{ $seq }{index1}  = $seq1;
                $store->{indx}{ $seq }{index2}  = $seq2;

                $store->{flow}{ $fid }{yield} += $read->{Yield};
                $store->{lane}{ $lid }{yield} += $read->{Yield};
                $store->{samp}{ $sid }{yield} += $read->{Yield};
                $store->{read}{ $rid }{yield} += $read->{Yield};
                $store->{indx}{ $seq }{yield} += $read->{Yield};

                $store->{flow}{ $fid }{yield_q30} += $read->{YieldQ30};
                $store->{lane}{ $lid }{yield_q30} += $read->{YieldQ30};
                $store->{samp}{ $sid }{yield_q30} += $read->{YieldQ30};
                $store->{read}{ $rid }{yield_q30} += $read->{YieldQ30};
                $store->{indx}{ $seq }{yield_q30} += $read->{YieldQ30};
            }
            
            my %bc_mismatch_counts = (
                'mm0' => $sample->{IndexMetrics}[0]{MismatchCounts}{0},
                'mm1' => $sample->{IndexMetrics}[0]{MismatchCounts}{1},
            );
            
            my @types = keys %bc_mismatch_counts;
            foreach my $mm ( @types ){
                my $count = 0;
                $count = $bc_mismatch_counts{ $mm } if defined $bc_mismatch_counts{ $mm };
                $store->{flow}{ $fid }{ $mm } += $count;
                $store->{lane}{ $lid }{ $mm } += $count;
                $store->{samp}{ $sid }{ $mm } += $count;
                $store->{indx}{ $seq }{ $mm } += $count;
            }
        }
    }
    addPrintFields(\%store, $fid);
    addGeneralStats(\%store, $fid, $cycle_string);
}

sub addGeneralStats{
    my ($store, $fid, $cycle_string) = @_;

    my $undet_perc = getPerc( $store->{'undt'}{'UNDETERMINED'}{'yield'}, $store->{'flow'}{$fid}{'yield'});
    my $run_overview_yield_factor = 1e6; # always report the run info in MBase
    $store->{'stats'}{'run_overview_string'} = sprintf "%s\t%s\t%s\t%s\t%s\t%s",
        round( $store->{'flow'}{$fid}{'yield'}, 0, $run_overview_yield_factor ),
        round( $store->{'undt'}{'UNDETERMINED'}{'yield'}, 0, $run_overview_yield_factor ),
        $store->{'flow'}{$fid}{'q30_print'},
        $store->{'flow'}{$fid}{'pf_print'},
        $cycle_string,
        round($undet_perc,1) . '%';

    $store->{'stats'}{'undet_perc'} = $undet_perc;
    $store->{'stats'}{'lane_count'} = scalar( keys %{$store->{'lane'}} );
    $store->{'stats'}{'samp_count'} = scalar( keys %{$store->{'samp'}} );
    $store->{'stats'}{'indx_count'} = scalar( keys %{$store->{'indx'}} );
    $store->{'stats'}{'identifier'} = join( "_", keys %{$store->{'flow'}} );
    $store->{'stats'}{'cycle_string'} = $cycle_string;
}

sub addPrintFields{
    my ($store, $fid) = @_;

    foreach my $type ( keys %$store ){
        foreach my $id ( keys %{$store->{ $type }} ){
            my $obj = $store->{ $type }{ $id };
            my $name = $obj->{ 'name' };

            $obj->{q30} = getPerc( $obj->{yield_q30}, $obj->{yield} );
            $obj->{yld_p} = getPerc( $obj->{yield}, $store->{flow}{ $fid }{yield} );
            $obj->{flowcell_print} = $fid;

            $obj->{ 'q30_print' } = round( $obj->{ 'q30' }, $ROUND_DECIMALS, 1 );
            $obj->{ 'yld_print' } = round( $obj->{ 'yield' }, 0, $YIELD_FACTOR );
            $obj->{ 'id_print' } = $id;
            $obj->{ 'name_print' } = $obj->{ 'name' };
            $obj->{ 'index1_print' } = $obj->{ 'index1' };
            $obj->{ 'index2_print' } = $obj->{ 'index2' };
            $obj->{ 'yld_p_print' } = round( $obj->{ 'yld_p' }, $ROUND_DECIMALS, 1 );

            ## percentage filtered does not exist for samples
            if ( exists $obj->{ 'clust_pf' } ){
                $obj->{'pf_print'} = round( getPerc( $obj->{'clust_pf'}, $obj->{'clust_raw'} ), $ROUND_DECIMALS );
            }

            ## some types do not have mismatch stats so skipping those from here
            next if $name =~ /read|UNDETERMINED/;
            next if $type eq "indx";

            $obj->{ 'total_reads' } = $obj->{ 'mm0' } + $obj->{ 'mm1' };
            $obj->{ 'mm0_print' } = 0;
            $obj->{ 'mm1_print' } = 0;
            if ( $obj->{ 'total_reads' } != 0 ){
                $obj->{ 'mm0_print' } = round( getPerc( $obj->{ 'mm0' }, $obj->{ 'total_reads' } ), $ROUND_DECIMALS );
                $obj->{ 'mm1_print' } = round( getPerc( $obj->{ 'mm1' }, $obj->{ 'total_reads' } ), $ROUND_DECIMALS );
            }

        }
    }
}

sub getPerc {
    my ($value, $total) = @_;
    
    if ( not defined($value) or not defined($total) ) {
        return 0;
    }
    elsif ($value < 0 or $total < 0) {
        die "[ERROR] Cannot calculate percentage if either value ($value) or total ($total) is < 0\n";
    }
    elsif ($value > $total){
        die "[ERROR] value ($value) should never be higher than total ($total)\n";
    }
    elsif ( $total == 0 and $value == 0) {
       return 0;
    }
    else {
        return $value*100/$total;
    }
}

sub printJson {
    my ($info, $output_file) = @_;
    
    my $coder = JSON->new->utf8->canonical->pretty;
    my $json_txt = $coder->encode($info);
    
    open my $fh, '>', $output_file or die "Unable to open output file ($output_file): $!\n";
        print $fh $json_txt;
    close $fh;
}

sub printTable {
    my ($info, $fields) = @_;

    say sprintf '## Settings: YieldFactor=%s, RoundDecimals=%s', 
      commify($YIELD_FACTOR),
      commify($ROUND_DECIMALS);
    say sprintf '## Flowcell: %s (%s, %d lanes, %d samples, %d indexes, %s cycles)', 
      $info->{'stats'}{'seq_runname'}, 
      $info->{'stats'}{'hmf_runname'}, 
      $info->{'stats'}{'lane_count'}, 
      $info->{'stats'}{'samp_count'},
      $info->{'stats'}{'indx_count'},
      $info->{'stats'}{'cycle_string'};
    
    say "#".join( $OUT_SEP, "level", @$fields );
    printTableForLevelSortedByName( $info->{'flow'}, $fields, 'RUN' );
    printTableForLevelSortedByName( $info->{'lane'}, $fields, 'LANE' );
    printTableForLevelSortedByYield( $info->{'indx'}, $fields, 'INDEX' );
    printTableForLevelSortedByName( $info->{'samp'}, $fields, 'SAMPLE' );
    printTableForLevelSortedByName( $info->{'read'}, $fields, 'READ' );
    printTableForLevelSortedByName( $info->{'undt'}, $fields, 'UNDET' );
}

sub printSummaryTable{
    my ($info, $fields) = @_;
    
    my @submissions = sort keys %{$info->{'stats'}{'submissions'}};
    map( $_ =~ s/HMFreg//, @submissions );
    my $submissions_string = join( ',', @submissions );
   
    say sprintf '## Settings: YieldFactor=%s, RoundDecimals=%s', 
      commify($YIELD_FACTOR),
      commify($ROUND_DECIMALS);
    say sprintf '## Flowcell: %s (%s, %d lanes, %d samples, %d indexes, %s cycles)', 
      $info->{'stats'}{'seq_runname'}, 
      $info->{'stats'}{'hmf_runname'}, 
      $info->{'stats'}{'lane_count'}, 
      $info->{'stats'}{'samp_count'},
      $info->{'stats'}{'indx_count'},
      $info->{'stats'}{'cycle_string'};
        
    say sprintf "## RunOverviewInfoLine: %s\t%s\t%s\t%s\t%s", 
      $info->{'stats'}{'hmf_runname'},
      $info->{'stats'}{'seq_runname'},
      $submissions_string,
      $info->{'stats'}{'run_overview_string'},
      $info->{'stats'}{'flowcell_qc'};
      
    say "#".join( $OUT_SEP, @$fields );
    printTableForLevelSortedByName( $info->{'samp'}, $fields );
}

sub printTableForLevelSortedByName{
    my ($info, $fields, $level) = @_;
    foreach my $id ( sort { $info->{$b}{'name'} cmp $info->{$a}{'name'} } keys %$info){
        my @output = map( $info->{ $id }{ $_."_print" } || $NA_CHAR, @$fields );
        unshift @output, $level if defined $level;
        say join( $OUT_SEP, @output );
    }
}

sub printTableForLevelSortedByYield{
    my ($info, $fields, $level) = @_;
    foreach my $id ( sort { $info->{$b}{'yld_print'} <=> $info->{$a}{'yld_print'} } keys %$info){
        my @output = map( $info->{ $id }{ $_."_print" } || $NA_CHAR, @$fields );
        unshift @output, $level if defined $level;
        say join( $OUT_SEP, @output );
    }
}

sub readSampleSheet{
    my ($csv_file) = @_;
    
    ## SampleSheet file has windows returns
    my $return_str = $/;
    $/ = "\r\n";
    
    my %output;
    $output{ 'samples' } = {};
    $output{ 'runname' } = 'NO_RUNNAME_FROM_SAMPLESHEET';
    $output{ 'platform' } = 'NO_PLATFORM_FROM_SAMPLESHEET';
    my @header;
    
    if ( ! -e $csv_file ){
        say "## WARNING skipping SampleSheet read: file not found ($csv_file)";
        return( \%output );
    }
    
    open FILE, "<", $csv_file or die "Couldn't open file ($csv_file): $!";
    while ( <FILE> ) {
        chomp($_);
        next if $_ =~ /^[\[\,]/;
        next if $_ eq "";
        my @fields = split( ",", $_);

        ## get hmf run id from config line
        if ($fields[0] =~ /Experiment(.)*Name/ ){
            my $run_name = $fields[1] || 'NA';
            $output{ 'runname' } = $run_name;
        }
        elsif ($fields[0] =~ /InstrumentType/ ){
            my $platform = $fields[1] || 'NA';
            $output{ 'platform' } = $platform;
        }

        ## find header
        elsif ( $_ =~ m/Sample_ID/ ){
            @header = @fields;
        }
        ## read sample line if header seen before
        elsif ( @header ){
            my %tmp = ();
            $tmp{ $_ } = shift @fields foreach @header;

            ## skip "empty" lines (where no sample was defined)
            my $sample_id_column = 'Sample_ID';
            next unless defined $tmp{ $sample_id_column } and $tmp{ $sample_id_column } ne '';

            my $sample_id = $tmp{ $sample_id_column };
            $output{ 'samples' }{ $sample_id } = \%tmp;
        }
    }
    close FILE;
    
    ## reset return string
    $/ = $return_str;
    
    return( \%output );
}

sub readQualityMetricsCsv{
    my ($file) = @_;

    my %qualityMetricsStore = ();
    my @header;

    open FILE, "<", $file or die "Couldn't open file ($file): $!";
    my $header_line = <FILE>;
    die "Header line has unexpected format" unless $header_line =~ m/^Lane,SampleID/;
    @header = split(",", $header_line);
    while ( <FILE> ) {
        chomp($_);
        next if $_ eq "";
        my @fields = split(",", $_);
        my %tmp = ();
        $tmp{$_} = shift @fields foreach @header;
        my $lane = $tmp{'Lane'};
        my $barcode = $tmp{'SampleID'};
        my $read = $tmp{'ReadNumber'};
        foreach my $field (@header){
            $qualityMetricsStore{$lane}{$barcode}{$read}{$field} = $tmp{$field};
        }
    }
    close FILE;
    return(\%qualityMetricsStore);
}

sub round{
    my ($number, $decimal, $factor) = @_;
    if ( not defined $number ){
        return $NA_CHAR;
    }
    $decimal = 0 unless defined $decimal;
    $factor = 1 unless defined $factor;
    my $rounded = sprintf("%.".$decimal."f", $number/$factor);
    return( $rounded );
}

## input "1000" gives output "1,000"
sub commify {
    local $_ = shift;
    $_ = int($_);
    1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
    return $_;
}

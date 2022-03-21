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

## -----
## Global variables
## -----
my $SCRIPT = basename $0;
my $OUT_SEP = "\t";
my $NA_CHAR = "NA";

my $SSHT_LOC = 'SampleSheet.csv';
my $RXML_LOC = 'RunInfo.xml';
my $JSON_LOC1 = 'Data/Intensities/BaseCalls/Stats/Stats.json';
my $JSON_LOC2 = 'Fastq/Stats/Stats.json';
my @OUT_FIELDS = qw( flowcell yld q30 pf yld_p mm0 mm1 id name submission index1 index2 );
my @SUM_FIELDS = qw( submission id name yld q30 index1 index2 );
my $ROUND_DECIMALS = 1; # Q30 with one decimal by default

## QC limits should be in sync with those from api (can be checked with cmd "api platforms")
my %SETTINGS_PER_PLATFORM = (
    'NovaSeq' => {
        'platform_name'    => "NOVASEQ",
        'min_flowcell_q30' => 85,
        'min_sample_yield' => 1e9,
        'max_undetermined' => 8,
        'yield_factor'     => 1e6,
        'run_name_regex'   => '^NO\d{2}-\d{4}$'
    },
    'NextSeq' => {
        'platform_name' => "NEXTSEQ",
        'min_flowcell_q30' => 75,
        'min_sample_yield' => 1e9,
        'max_undetermined' => 50,
        'yield_factor'     => 1e6,
        'run_name_regex'   => '^NS\d{2}-\d{4}$'
    },
    'ISeq' => {
        'platform_name' => "ISEQ",
        'min_flowcell_q30' => 75,
        'min_sample_yield' => 1e6,
        'max_undetermined' => 50,
        'yield_factor'     => 1,
        'run_name_regex'   => '^IS\d{2}-\d{4}$'
    },
    # 'HiSeq' => {
    #     'platform_name' => "HISEQ",
    #     'min_flowcell_q30' => 75,
    #     'min_sample_yield' => 1e9,
    #     'max_undetermined' => 8,
    #     'yield_factor'     => 1e6
    # },
);
my $RUN_PATH;
my $JSON_PATH;
my $RXML_PATH;
my $SSHT_PATH;
my $OUT_JSON_PATH;

my $HELP =<<HELP;

  Description
    Parses conversion output json and prints table with
    flowcell/lane/sample/read metrics
    
  Usage
    $SCRIPT -run_dir \${run-path}
    $SCRIPT -run_dir \${run-path} -samplesheet \${samplesheet-path}
    $SCRIPT -json_path \${json-path}
    
  Options
    -sep <s>               Output sep (default = <TAB>)
    -yield_factor <i>      Factor to divide all yields with (default present per platform)
    -round_decimals <i>    Decimals to keep for q30 (default: $ROUND_DECIMALS)
    -samplesheet <s>       Path to SampleSheet.csv file
    -run_info_xml <s>      Path to RunInfo.xml file
    -json_out <s>          Path to output json file (only written if provided)
    -no_qc                 Skip QC checks
    -summary               Prints extra sample summary table
    -debug                 Prints complete datastructure
    
HELP
print $HELP and exit(0) if scalar @ARGV == 0;

## -----
## Gather input
## -----
my %opt = ();
GetOptions (
  "run_dir=s"          => \$RUN_PATH,
  "json_path=s"        => \$JSON_PATH,
  "samplesheet=s"      => \$SSHT_PATH,
  "run_info_xml=s"     => \$RXML_PATH,
  "json_out=s"         => \$OUT_JSON_PATH,
  "sep=s"              => \$OUT_SEP,
  "round_decimals=i"   => \$ROUND_DECIMALS,
  "yield_factor=i"     => \$opt{ yield_factor },
  "no_qc"              => \$opt{ no_qc },
  "summary"            => \$opt{ print_summary },
  "debug"              => \$opt{ debug },
  "help|h"             => \$opt{ help },
) or die "[ERROR] Issue in command line arguments\n";
print $HELP and exit(0) if $opt{ help };
die "[ERROR] Provide either run-dir or json-path not both\n" if (defined $RUN_PATH and defined $JSON_PATH);
die "[ERROR] Provide either run-dir or json-path, see -h\n" unless (defined $RUN_PATH or defined $JSON_PATH);
die "[ERROR] Provided run dir does not exist ($RUN_PATH)\n" if defined $RUN_PATH and not -d $RUN_PATH;
die "[ERROR] Provided json does not exist ($JSON_PATH)\n" if defined $JSON_PATH and not -f $JSON_PATH;
die "[ERROR] Provided xml does not exist ($RXML_PATH)\n" if defined $RXML_PATH and not -f $RXML_PATH;

$RXML_PATH = "$RUN_PATH/$RXML_LOC" if defined $RUN_PATH;
$SSHT_PATH = "$RUN_PATH/$SSHT_LOC" if (defined $RUN_PATH and not defined $SSHT_PATH);

$JSON_PATH = "$RUN_PATH/$JSON_LOC1" if defined $RUN_PATH;
$JSON_PATH = "$RUN_PATH/$JSON_LOC2" if not -f $JSON_PATH;

## -----
## MAIN
## -----

my $seq_run = basename $RUN_PATH;
my $ssht_info = readSampleSheet( $SSHT_PATH );

## Need to make sure we know the exact platform
my $platform = determinePlatformByRunName($ssht_info->{'runname'}, \%SETTINGS_PER_PLATFORM);
my $SETTINGS = $SETTINGS_PER_PLATFORM{$platform};
$SETTINGS->{yield_factor} = $opt{yield_factor} if defined $opt{yield_factor};
my $YIELD_FACTOR = $SETTINGS->{yield_factor};

my $json_info = readJson( $JSON_PATH );
my $rxml_info = readXml( $RXML_PATH );
my $parsed_info = parseJsonInfo( $json_info, $rxml_info );
$parsed_info->{stats}{platform} = $platform;

addSamplesheetInfo( $parsed_info, $ssht_info, $seq_run );
performQC( $parsed_info, $SETTINGS) unless $opt{ no_qc };

if ( $opt{ print_summary } ){
    printSummaryTable( $parsed_info, \@SUM_FIELDS );
}
else{
    printTable( $parsed_info, \@OUT_FIELDS );
    printJson( $parsed_info, $OUT_JSON_PATH ) if defined $OUT_JSON_PATH;
}

if ( $opt{ debug } ){
    say "[DEBUG] Samplesheet data structure";
    print Dumper $ssht_info;
    say "[DEBUG] Complete final data structure";
    print Dumper $parsed_info;
}

## -----
## /MAIN
## -----

sub addSamplesheetInfo{
    my ($json_info, $ssht_info, $seq_run) = @_;

    ## add info to stats (eg hmf_runname could be X17-0001)
    my $hmf_runname = $ssht_info->{'runname'};
    $json_info->{ 'stats' }{ 'seq_runname' } = $seq_run;
    $json_info->{ 'stats' }{ 'hmf_runname' } = $hmf_runname;
    $json_info->{ 'stats' }{ 'submissions' } = {};
    
    ## override flowcell name with ExperimentName from SampleSheet
    my @fcids = keys %{$json_info->{'flow'}};
    die "[ERROR] There should only be one flowcell in json_info" unless scalar @fcids == 1;
    my $fcid = $fcids[0];
    $json_info->{ 'flow' }{ $fcid }{ 'name' } = $hmf_runname;
    $json_info->{ 'flow' }{ $fcid }{ 'name_print' } = $hmf_runname;
    
    ## add sample metadata
    my $samples = $json_info->{ 'samp' };
    foreach my $sample_id ( keys %$samples ){
        my $sample = $samples->{ $sample_id };
        
        my $submission = $NA_CHAR;
        $submission = $ssht_info->{'samples'}{$sample_id}{ 'Sample_Project' } if defined $ssht_info->{'samples'}{$sample_id}{ 'Sample_Project' };
        $sample->{ 'submission_print' } = $submission;
        $json_info->{ 'stats' }{ 'submissions' }{ $submission } = 1;
        
        my $description = $NA_CHAR;
        $description = $ssht_info->{'samples'}{$sample_id}{ 'Description' } if defined $ssht_info->{'samples'}{$sample_id}{ 'Description' };
        $sample->{ 'description_print' } = $description;
    }
}

sub readXml{
    my ($file) = @_;
    my $obj = XMLin( $file );
    return($obj);
}

sub readJson{
    my ($json_file) = @_;
    my $json_txt = read_file( $json_file );
    my $json_obj = decode_json( $json_txt );
    return( $json_obj );
}

sub determinePlatformByRunName {
    my ($runname, $settings_per_platform) = @_;
    my $final_platform = "";

    ## try to match by regex
    while(my($platform_option, $settings) = each $settings_per_platform) {
        if ($runname =~ m/$settings->{ 'run_name_regex' }/) {
            if ($final_platform ne "") {
                die "Multiple platforms match run name regex ($final_platform and $platform_option)\n"
            }
            $final_platform = $platform_option;
            say "## INFO: platform configured to $platform_option (based on input '$runname')";
        }
    }

    if ($final_platform eq "") {
        die "[ERROR] Unable to determine platform from run name ($runname)\n";
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
    my ($raw_json_info, $run_xml_info) = @_;
    my %info = ();

    ## Reading phase
    my @cycle_counts = map( $_->{NumCycles}, @{$run_xml_info->{Run}{Reads}{Read}});
    my @read_cycle_counts = map( $_->{IsIndexedRead} eq "N" ? $_->{NumCycles} : (), @{$run_xml_info->{Run}{Reads}{Read}});
    my $total_non_index_cycle_count = sum(@read_cycle_counts);

    my $cycle_string = join( "|", @cycle_counts);
    my $fid = $raw_json_info->{ 'Flowcell' };
    $info{ 'flow' }{ $fid }{ 'id' } = $fid;
    $info{ 'flow' }{ $fid }{ 'name' } = $raw_json_info->{ 'RunId' };

    ## First reading the "unknown barcodes" and add them to "index sequences"
    my $unknowns = $raw_json_info->{'UnknownBarcodes'};
    foreach my $lane ( @$unknowns ){
        my $lid = join( "", "lane", $lane->{ Lane } );
        my $unknown_barcodes = $lane->{'Barcodes'};
        foreach my $bc ( keys %$unknown_barcodes ){
            $info{ indx }{ $bc }{ name } = 'IndexFromUnknown';
            my $seq1 = (split(/\+/, $bc))[0] || $NA_CHAR;
            my $seq2 = (split(/\+/, $bc))[1] || $NA_CHAR;
            $info{ indx }{ $bc }{ index1 } = $seq1;
            $info{ indx }{ $bc }{ index2 } = $seq2;
            # Unlike actual samples, the unknowns are reported as cluster counts instead of yield
            # So need to calculate the yield using non-index cycle counts from RunInfo.xml
            $info{ indx }{ $bc }{ yield } += $unknown_barcodes->{ $bc } * $total_non_index_cycle_count;
        }
    }

    ## Then read samples
    my $lanes = $raw_json_info->{'ConversionResults'};
    foreach my $lane ( @$lanes ){
        my $lid = join( "", "lane", $lane->{LaneNumber} );
        $info{lane}{ $lid }{name} = $lid;
        
        $info{flow}{ $fid }{clust_raw} += $lane->{TotalClustersRaw};
        $info{flow}{ $fid }{clust_pf} += $lane->{TotalClustersPF};
        $info{lane}{ $lid }{clust_raw} += $lane->{TotalClustersRaw};
        $info{lane}{ $lid }{clust_pf} += $lane->{TotalClustersPF};
        
        ## Undetermined info is stored separate from samples in json
        my $undet_id = 'UNDETERMINED';
        my $undet_obj = $lane->{Undetermined};
        my $undet_reads = $undet_obj->{ReadMetrics};
        my $undet_info = \%{$info{undt}{ $undet_id }};
        $undet_info->{name} = $undet_id;
        foreach my $read ( @$undet_reads ){
            my $rid = join( "", "read", $read->{ReadNumber} );
            $undet_info->{yield} += $read->{Yield};
            $undet_info->{yield_q30} += $read->{YieldQ30};
            $info{flow}{ $fid }{yield} += $read->{Yield};
            $info{flow}{ $fid }{yield_q30} += $read->{YieldQ30};
            $info{lane}{ $lid }{yield} += $read->{Yield};
            $info{lane}{ $lid }{yield_q30} += $read->{YieldQ30};
            $info{read}{ $rid }{yield} += $read->{Yield};
            $info{read}{ $rid }{yield_q30} += $read->{YieldQ30};
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
            $info{samp}{ $sid }{name} = $snm;
            $info{samp}{ $sid }{index1}  = $seq1;
            $info{samp}{ $sid }{index2}  = $seq2;
            
            my $reads = $sample->{ReadMetrics};
            foreach my $read ( @$reads ){
                my $rid = join( "", "read", $read->{ReadNumber} );
                $info{read}{ $rid }{name} = $rid;
                $info{indx}{ $seq }{name} = 'IndexFromSample';
                $info{indx}{ $seq }{index1}  = $seq1;
                $info{indx}{ $seq }{index2}  = $seq2;
                
                $info{flow}{ $fid }{yield} += $read->{Yield};
                $info{lane}{ $lid }{yield} += $read->{Yield};
                $info{samp}{ $sid }{yield} += $read->{Yield};
                $info{read}{ $rid }{yield} += $read->{Yield};
                $info{indx}{ $seq }{yield} += $read->{Yield};
                
                $info{flow}{ $fid }{yield_q30} += $read->{YieldQ30};
                $info{lane}{ $lid }{yield_q30} += $read->{YieldQ30};
                $info{samp}{ $sid }{yield_q30} += $read->{YieldQ30};
                $info{read}{ $rid }{yield_q30} += $read->{YieldQ30};
                $info{indx}{ $seq }{yield_q30} += $read->{YieldQ30};
            }
            
            my %bc_mismatch_counts = (
                'mm0' => $sample->{IndexMetrics}[0]{MismatchCounts}{0},
                'mm1' => $sample->{IndexMetrics}[0]{MismatchCounts}{1},
            );
            
            my @types = keys %bc_mismatch_counts;
            foreach my $mm ( @types ){
                my $count = 0;
                $count = $bc_mismatch_counts{ $mm } if defined $bc_mismatch_counts{ $mm };
                $info{flow}{ $fid }{ $mm } += $count;
                $info{lane}{ $lid }{ $mm } += $count;
                $info{samp}{ $sid }{ $mm } += $count;
                $info{indx}{ $seq }{ $mm } += $count;
            }
        }
    }

    ## Create the info to print later
    foreach my $type ( keys %info ){
        foreach my $id ( keys %{$info{ $type }} ){
            my $obj = $info{ $type }{ $id };
            my $name = $obj->{ 'name' };
            
            $obj->{q30} = getPerc( $obj->{yield_q30}, $obj->{yield} );
            $obj->{yld_p} = getPerc( $obj->{yield}, $info{flow}{ $fid }{yield} );
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
    
    ## Collect some general stats/info
    my $undet_perc = getPerc( $info{'undt'}{'UNDETERMINED'}{'yield'}, $info{'flow'}{$fid}{'yield'});
    my $run_overview_yield_factor = 1e6; # always report the run info in MBase
    $info{'stats'}{'run_overview_string'} = sprintf "%s\t%s\t%s\t%s\t%s\t%s", 
      round( $info{'flow'}{$fid}{'yield'}, 0, $run_overview_yield_factor ),
      round( $info{'undt'}{'UNDETERMINED'}{'yield'}, 0, $run_overview_yield_factor ),
      $info{'flow'}{$fid}{'q30_print'},
      $info{'flow'}{$fid}{'pf_print'},
      $cycle_string,
      round($undet_perc,1) . '%';
      
    $info{'stats'}{'undet_perc'} = $undet_perc;
    $info{'stats'}{'lane_count'} = scalar( keys %{$info{'lane'}} );
    $info{'stats'}{'samp_count'} = scalar( keys %{$info{'samp'}} );
    $info{'stats'}{'indx_count'} = scalar( keys %{$info{'indx'}} );
    $info{'stats'}{'identifier'} = join( "_", keys %{$info{'flow'}} );
    $info{'stats'}{'cycle_string'} = $cycle_string;

    return \%info;
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
    
    my $coder = JSON->new->utf8->canonical;
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

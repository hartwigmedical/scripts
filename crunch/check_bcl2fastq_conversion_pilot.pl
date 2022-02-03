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
use Text::CSV qw/csv/;
use POSIX qw(strftime);
use 5.010.000;

my $SCRIPT = basename $0;
my $OUT_SEP = "\t";

use constant NA => 'NA';
use constant UNDETERMINED => 'Undetermined';

my @OUT_FIELDS = qw(flowcell yld q30 yld_p mm0 mm1 mm2 id name submission index1 index2);
my @SUM_FIELDS = qw(submission id name yld q30 index1 index2);
my $ROUND_DECIMALS = 1;

my $RUN_PATH;
my $OUT_JSON_PATH;

my $QCSV_PATH = 'Fastq/Reports/Quality_Metrics.csv';
my $DCSV_PATH = 'Fastq/Reports/Demultiplex_Stats.csv';
my $RXML_PATH = 'Fastq/Reports/RunInfo.xml';
my $SSHT_PATH = 'Fastq/Reports/SampleSheet.csv';

my %SETTINGS_PER_PLATFORM = (
    'NovaSeq' => {'min_flowcell_q30' => 85, 'min_sample_yield' => 1e9, 'max_undetermined' =>  8, 'yield_factor' => 1e6},
    'NextSeq' => {'min_flowcell_q30' => 75, 'min_sample_yield' => 1e9, 'max_undetermined' => 50, 'yield_factor' => 1e6},
    'ISeq'    => {'min_flowcell_q30' => 75, 'min_sample_yield' => 1e6, 'max_undetermined' => 50, 'yield_factor' => 1},
    'HiSeq'   => {'min_flowcell_q30' => 75, 'min_sample_yield' => 1e9, 'max_undetermined' =>  8, 'yield_factor' => 1e6},
);
my @KNOWN_PLATFORMS = keys %SETTINGS_PER_PLATFORM;

my $HELP =<<HELP;

  Description:
    Parses conversion metrics from bclconvert output files and performs QC
    
  Usage:
    $SCRIPT -runDir \${run-path}
    - OR -
    $SCRIPT -sampleSheetCsv SampleSheet.csv -runInfoXml RunInfo.xml -qualityCsv Quality_Metrics.csv -demuxCsv Demultiplex_Stats.csv
    
  Options:
    -sampleSheetCsv <s>   Path to SampleSheet.csv file
    -runInfoXml     <s>   Path to RunInfo.xml file
    -qualityCsv     <s>   Path to Quality_Metrics.csv file
    -demuxCsv       <s>   Path to Demultiplex_Stats.csv file
    -outputJson     <s>   Path to output json file to create
    -sep            <s>   Output sep [default = <TAB>]
    -yieldFactor    <i>   Factor to divide all yields with [default differs per platform]
    -decimals       <i>   Decimals to keep for q30 [$ROUND_DECIMALS]
    -noQc                 Skip QC checks and just print table
    -summary              Prints extra sample summary table

HELP
print $HELP and exit(0) if scalar @ARGV == 0;

my %opt = ();
GetOptions (
    "runDir=s"      => \$RUN_PATH,
    "qualityCsv=s"  => \$QCSV_PATH,
    "demuxCsv=s"    => \$DCSV_PATH,
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
    $QCSV_PATH = "$RUN_PATH/$QCSV_PATH";
    $DCSV_PATH = "$RUN_PATH/$DCSV_PATH";
    $seq_run = basename $RUN_PATH;
}
die "[ERROR] Provided xml does not exist ($RXML_PATH)\n" if not -f $RXML_PATH;
die "[ERROR] Provided sample sheet does not exist ($SSHT_PATH)\n" if not -f $SSHT_PATH;

my $ssht_info = readSampleSheet( $SSHT_PATH );
my $platform = determinePlatformByString($ssht_info->{'platform'}, \@KNOWN_PLATFORMS);
my $settings = $SETTINGS_PER_PLATFORM{$platform};
$settings->{yield_factor} = $opt{yield_factor} if defined $opt{yield_factor};
my $YIELD_FACTOR = $settings->{yield_factor};

my $run_info = readRunInfoXml($RXML_PATH);
my $quality_info = readBclConvertQualityMetricsCsv($QCSV_PATH, $run_info);
my $demux_info = readBclConvertDemuxCsv($DCSV_PATH, $run_info);
my $store = aggregateBclconvertInfo($quality_info, $demux_info, $run_info);

$store->{stats}{platform} = $platform;

addSamplesheetInfo($store, $ssht_info, $seq_run);
performQC($store, $settings) unless $opt{no_qc};

if ( $opt{print_summary} ){
    printSummaryTable($store, \@SUM_FIELDS );
}
else{
    printTable($store, \@OUT_FIELDS);
    printJson($store, $OUT_JSON_PATH) if defined $OUT_JSON_PATH;
}

if ( $opt{debug} ){
    printJson($ssht_info, 'tmp_debug_ssheet.json');
    printJson($quality_info, "tmp_debug_quality.json");
    printJson($demux_info, "tmp_debug_demux.json");
    printJson($store, 'tmp_debug_store.json');
}

sub aggregateBclconvertInfo{
    my ($qual, $demux, $runinfo) = @_;
    my %out = ();

    my $cycle_string = $runinfo->{cycle_string};
    my $fid = $runinfo->{flowcell_id};

    $out{flowcell}{$fid}{id} = $fid;
    $out{flowcell}{$fid}{name} = "FlowcellNameTODO";

    while (my ($type, $objects) = each %$qual) {
        while (my ($id, $object) = each %$objects) {
            while (my ($qkey, $qval) = each %$object) {
                $out{$type}{$id}{$qkey} = $qval;
                if (exists $demux->{$type}{$id}) {
                    while (my ($dkey, $dval) = each %{$demux->{$type}{$id}}) {
                        $out{$type}{$id}{$dkey} = $dval;
                    }
                }
            }
        }
    }
    addPrintFields(\%out, $fid);
    addGeneralStats(\%out, $fid, $cycle_string);
    return \%out;
}

sub addSamplesheetInfo{
    my ($data, $ssht_info, $seq_run) = @_;

    # Add info to stats (eg hmf_runname could be X17-0001)
    my $hmf_runname = $ssht_info->{runname};
    $data->{stats}{seq_runname} = $seq_run;
    $data->{stats}{hmf_runname} = $hmf_runname;
    $data->{stats}{submissions} = {};
    
    # Override flowcell name with ExperimentName from SampleSheet
    my @fcids = keys %{$data->{flowcell}};
    scalar @fcids == 1 || die "[ERROR] There should be exactly one flowcell in data\n";
    my $fcid = $fcids[0];
    $data->{flowcell}{$fcid}{name} = $hmf_runname;
    $data->{flowcell}{$fcid}{name_print} = $hmf_runname;

    # Add sample metadata
    my %submissions_encountered = ();
    my $samples = $data->{samples};
    foreach my $barcode (keys %$samples){
        my $sample = $samples->{$barcode};

        my $submission = $ssht_info->{samples}{$barcode}{Sample_Project} || NA;
        $sample->{submission} = $submission;
        $sample->{submission_print} = $submission;
        $submissions_encountered{$submission} = 1 unless $barcode eq UNDETERMINED;

        my $name = $ssht_info->{samples}{$barcode}{Sample_Name} || NA;
        $sample->{name} = $name;
        $sample->{name_print} = $name;

        my $description = $ssht_info->{samples}{$barcode}{Description} || NA;
        $sample->{description} = $description;
        $sample->{description_print} = $description;
    }

    my @submissions = sort keys %submissions_encountered;
    $data->{stats}{submissions} = \@submissions;
}

sub readRunInfoXml{
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

sub determinePlatformByString {
    my ($platform_string, $platforms) = @_;
    my $final_platform = 'UnknownPlatform';

    # Exact platform string varies too much so try to match by regex
    foreach my $known (@$platforms) {
        if ($platform_string =~ m/^$known/i) {
            $final_platform = $known;
        }
    }
    return $final_platform;
}

sub performQC{
    my ($info, $qc_limits) = @_;
   
    my $stats = $info->{stats};
    my $samples = $info->{samples};
    my $lanes = $info->{lanes};
    my $fails = 0;
    my $identifier = $stats->{identifier};

    # Flowcell checks
    my $undet = $stats->{undet_perc};
    my $max_undet = $qc_limits->{max_undetermined};
    if ( $undet > $max_undet ){
        warn "## WARNING Percentage undetermined ($undet) too high (max=$max_undet)\n";
        $fails += 1;
    }
        
    # Lane and sample checks
    $fails += checkObjectField( $lanes, 'q30',   $qc_limits->{min_flowcell_q30} );
    $fails += checkObjectField( $samples, 'yield', $qc_limits->{min_sample_yield} );
    
    # Conclusion
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
    $stats->{flowcell_qc} = $final_result;
}

sub checkObjectField{
    my ($objects, $field, $min) = @_;
    my $fails = 0;
    foreach my $obj_key ( sort { $objects->{$b}{name} cmp $objects->{$a}{name} } keys %$objects){
        my $obj = $objects->{$obj_key};
        my $name = $obj->{name};
        next if $name eq UNDETERMINED;
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

sub addGeneralStats{
    my ($store, $fid, $cycle_string) = @_;

    my $undet_perc = getPerc( $store->{samples}{Undetermined}{yield}, $store->{flowcell}{$fid}{yield});
    my $run_overview_yield_factor = 1e6; # always report the run info in MBase
    $store->{stats}{run_overview_string} = sprintf "%s\t%s\t%s\t%s\t%s\t%s",
        round( $store->{flowcell}{$fid}{yield}, 0, $run_overview_yield_factor ),
        round( $store->{samples}{Undetermined}{yield}, 0, $run_overview_yield_factor ),
        $store->{flowcell}{$fid}{q30_print},
        NA, # placeholder for "passing filter percentage" was available in bcl2fastq output but not in bclconvert output
        $cycle_string,
        round($undet_perc,1) . '%';

    $store->{stats}{undet_perc} = $undet_perc;
    $store->{stats}{lane_count} = scalar( keys %{$store->{'lanes'}} );
    $store->{stats}{samp_count} = scalar( keys %{$store->{'samples'}} );
    $store->{stats}{identifier} = join( "_", keys %{$store->{'flowcell'}} );
    $store->{stats}{cycle_string} = $cycle_string;
}

sub addPrintFields{
    my ($store, $fid) = @_;

    foreach my $type ( keys %$store ){
        foreach my $id ( keys %{$store->{$type}} ){
            my $obj = $store->{$type}{$id};

            $obj->{q30} = getPerc($obj->{yield_q30}, $obj->{yield});
            $obj->{yld_p} = getPerc($obj->{yield}, $store->{flowcell}{$fid}{yield});

            $obj->{q30_print} = round($obj->{q30}, $ROUND_DECIMALS, 1);
            $obj->{yld_print} = round($obj->{yield}, 0, $YIELD_FACTOR);
            $obj->{yld_p_print} = round( $obj->{yld_p}, $ROUND_DECIMALS, 1 );

            $obj->{flowcell_print} = $fid;
            $obj->{id_print} = $id;
            $obj->{name_print} = $obj->{name};
            $obj->{index1_print} = $obj->{index1};
            $obj->{index2_print} = $obj->{index2};

            # Some types do not have mismatch stats so skipping those from here
            next unless $type =~ /samples|lanes|flowcell/;

            $obj->{total_reads} = $obj->{mm0} + $obj->{mm1} + $obj->{mm2};
            $obj->{mm0_print} = 0;
            $obj->{mm1_print} = 0;
            $obj->{mm2_print} = 0;
            if ( $obj->{total_reads} != 0 ){
                $obj->{mm0_print} = round(getPerc($obj->{mm0}, $obj->{total_reads}), $ROUND_DECIMALS);
                $obj->{mm1_print} = round(getPerc($obj->{mm1}, $obj->{total_reads}), $ROUND_DECIMALS);
                $obj->{mm2_print} = round(getPerc($obj->{mm2}, $obj->{total_reads}), $ROUND_DECIMALS);
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
    say "[INFO] Written JSON file $output_file"
}

sub printTable {
    my ($info, $fields) = @_;

    printStandardHeaderLines($info);
    
    say "#".join( $OUT_SEP, "level", @$fields );
    printTableForLevelSortedByName($info->{flowcell}, $fields, 'RUN');
    printTableForLevelSortedByName($info->{lanes}, $fields, 'LANE');
    printTableForLevelSortedByName($info->{samples}, $fields, 'SAMPLE');
    printTableForLevelSortedByName($info->{reads}, $fields, 'READ');
    printTableForLevelSortedByName($info->{undetermined}, $fields, 'UNDET');
}

sub printSummaryTable{
    my ($info, $fields) = @_;
    
    my @submissions = sort @{$info->{stats}{submissions}};
    map( $_ =~ s/HMFreg//, @submissions );
    my $submissions_string = join( ',', @submissions );

    printStandardHeaderLines($info);
        
    say sprintf "## RunOverview INFO: %s\t%s\t%s\t%s\t%s",
      $info->{stats}{hmf_runname},
      $info->{stats}{seq_runname},
      $submissions_string,
      $info->{stats}{run_overview_string},
      $info->{stats}{flowcell_qc};
      
    say "#".join($OUT_SEP, @$fields);
    printTableForLevelSortedByName($info->{samples}, $fields);
}

sub printStandardHeaderLines{
    my ($data) = @_;
    say sprintf '## Settings INFO: Platform=%s|YieldFactor=%s|RoundDecimals=%s',
        $data->{stats}{platform},
        commify($YIELD_FACTOR),
        commify($ROUND_DECIMALS);

    say sprintf '## Flowcell INFO: %s (%s, %d lanes, %d samples, %s cycles)',
        $data->{stats}{seq_runname},
        $data->{stats}{hmf_runname},
        $data->{stats}{lane_count},
        $data->{stats}{samp_count},
        $data->{stats}{cycle_string};
}

sub printTableForLevelSortedByName{
    my ($info, $fields, $level) = @_;
    foreach my $id ( sort { $info->{$b}{name} cmp $info->{$a}{name} } keys %$info){
        my @output = map( $info->{ $id }{ $_."_print" } || NA, @$fields );
        unshift @output, $level if defined $level;
        say join( $OUT_SEP, @output );
    }
}

sub printTableForLevelSortedByYield{
    my ($info, $fields, $level) = @_;
    foreach my $id ( sort { $info->{$b}{yld_print} <=> $info->{$a}{yld_print} } keys %$info){
        my @output = map( $info->{ $id }{ $_."_print" } || NA, @$fields );
        unshift @output, $level if defined $level;
        say join( $OUT_SEP, @output );
    }
}

sub readSampleSheet{
    my ($csv_file) = @_;
    
    # SampleSheet file has windows returns
    my $return_str = $/;
    $/ = "\r\n";
    
    my %output;
    $output{samples} = {};
    $output{runname} = 'NO_RUNNAME_FROM_SAMPLESHEET';
    $output{platform} = 'NO_PLATFORM_FROM_SAMPLESHEET';
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

        # Get hmf run id from config line
        if ($fields[0] =~ /Experiment(.)*Name/ ){
            my $run_name = $fields[1] || 'NA';
            $output{runname} = $run_name;
        }
        elsif ($fields[0] =~ /InstrumentType/ ){
            my $platform = $fields[1] || 'NA';
            $output{platform} = $platform;
        }
        elsif ( $_ =~ m/Sample_ID/ ){
            # Find header
            @header = @fields;
        }
        elsif ( @header ){
            # Read sample line but only if header seen before
            my %tmp = ();
            $tmp{ $_ } = shift @fields foreach @header;

            # Skip lines without sample info
            my $sample_id_column = 'Sample_ID';
            next unless defined $tmp{ $sample_id_column } and $tmp{ $sample_id_column } ne '';

            my $sample_id = $tmp{ $sample_id_column };
            $output{samples}{ $sample_id } = \%tmp;
        }
    }
    close FILE;
    
    # Reset return string
    $/ = $return_str;
    
    return( \%output );
}

sub readBclConvertQualityMetricsCsv{
    my ($csv, $run) = @_;
    my $fcid = $run->{flowcell_id};

    my $lines = csv(in => $csv, headers => 'auto');
    my %value_keys = (
        'Yield' => 'yield',
        'YieldQ30' => 'yield_q30'
    );
    my %quality_store = ();

    foreach my $line (@$lines){
        my $barcode = $line->{SampleID};
        my $lane = 'lane' . $line->{Lane};
        my $read = 'read' . $line->{ReadNumber};

        $quality_store{lanes}{$lane}{name} = $lane;
        $quality_store{reads}{$read}{name} = $read;
        $quality_store{samples}{$barcode}{name} = $barcode;

        if ($barcode ne UNDETERMINED){
            $quality_store{samples}{$barcode}{index1} = $line->{index};
            $quality_store{samples}{$barcode}{index2} = $line->{index2};
        }

        foreach my $value_key (keys %value_keys){
            my $store_key = $value_keys{$value_key};
            $quality_store{flowcell}{$fcid}{$store_key} += $line->{$value_key};
            $quality_store{lanes}{$lane}{$store_key} += $line->{$value_key};
            $quality_store{reads}{$read}{$store_key} += $line->{$value_key};
            $quality_store{samples}{$barcode}{$store_key} += $line->{$value_key};
        }
    }
    return(\%quality_store);
}

sub readBclConvertDemuxCsv{
    my ($csv, $run) = @_;
    my $fcid = $run->{flowcell_id};

    my $lines = csv(in => $csv, headers => 'auto');
    my %value_keys = (
        '# Perfect Index Reads'      => 'mm0',
        '# One Mismatch Index Reads' => 'mm1',
        '# Two Mismatch Index Reads' => 'mm2',
        '# Reads'                    => 'reads'
    );
    my %demux_store;

    foreach my $line (@$lines){
        my $lane = 'lane' . $line->{Lane};
        my $barcode = $line->{SampleID};
        foreach my $value_key (keys %value_keys){
            my $store_key = $value_keys{$value_key};
            $demux_store{flowcell}{$fcid}{$store_key} += $line->{$value_key};
            $demux_store{lanes}{$lane}{$store_key} += $line->{$value_key};
            $demux_store{samples}{$barcode}{$store_key} += $line->{$value_key};
        }

    }
    return \%demux_store;
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

sub round{
    my ($number, $decimal, $factor) = @_;
    if ( not defined $number ){
        return NA;
    }
    $decimal = 0 unless defined $decimal;
    $factor = 1 unless defined $factor;
    my $rounded = sprintf("%.".$decimal."f", $number/$factor);
    return( $rounded );
}

sub commify {
    # Input "1000" produces output "1,000"
    local $_ = shift;
    $_ = int($_);
    1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
    return $_;
}

sub info{
    my ($msg) = @_;
    say "[INFO] " . (strftime "%y%m%d %H:%M:%S", localtime) . " - " . $msg;
}

sub warn{
    my ($msg) = @_;
    warn "[WARN] " . (strftime "%y%m%d %H:%M:%S", localtime) . " - " . $msg . "\n";
}

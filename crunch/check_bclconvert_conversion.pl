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

my @OUT_FIELDS = qw(flowcell yld q30 yld_p mm0 mm1 id name submission index1 index2);
my @SUM_FIELDS = qw(submission id name yld q30 index1 index2);
my $ROUND_DECIMALS = 1;

my $REL_QCSV_PATH = 'Fastq/Reports/Quality_Metrics.csv';
my $REL_DCSV_PATH = 'Fastq/Reports/Demultiplex_Stats.csv';
my $REL_RXML_PATH = 'Fastq/Reports/RunInfo.xml';
my $REL_SSHT_PATH = 'Fastq/Reports/SampleSheet.csv';

my $NOVASEQ = "NovaSeq";
my $NEXTSEQ = "NextSeq";
my $ISEQ = "ISeq";
my $HISEQ = "HiSeq";

# The official QC settings are in the platforms endpoint of HMFAPI (the copy here is for convenience)
my %SETTINGS_PER_PLATFORM = (
    $NOVASEQ => {'min_flowcell_q30' => 85, 'min_sample_yield' => 1e9, 'max_undetermined' =>  8, 'yield_factor' => 1e6},
    $NEXTSEQ => {'min_flowcell_q30' => 75, 'min_sample_yield' => 1e9, 'max_undetermined' => 50, 'yield_factor' => 1e6},
    $ISEQ    => {'min_flowcell_q30' => 75, 'min_sample_yield' => 5e6, 'max_undetermined' => 50, 'yield_factor' => 1},
);
my @KNOWN_PLATFORMS = keys %SETTINGS_PER_PLATFORM;

my $HELP =<<HELP;

  Description:
    Parses conversion metrics from bclconvert output and performs the flowcell QC.
    
  Usage:
    $SCRIPT -runDir \${run-path}
    - OR -
    $SCRIPT -sampleSheet SampleSheet.csv -runInfoXml RunInfo.xml -qualityCsv Quality_Metrics.csv -demuxCsv Demultiplex_Stats.csv
    
  Options:
    -sampleSheet    <s>   Path to SampleSheet.csv file
    -runInfoXml     <s>   Path to RunInfo.xml file
    -qualityCsv     <s>   Path to Quality_Metrics.csv file
    -demuxCsv       <s>   Path to Demultiplex_Stats.csv file
    -outputJson     <s>   Path to output json file to create
    -sep            <s>   Output sep [default = <TAB>]
    -yieldFactor    <i>   Factor to divide all yields with [default differs per platform]
    -decimals       <i>   Decimals to keep for q30 [$ROUND_DECIMALS]
    -noQc                 Skip QC checks and just print table

HELP
print $HELP and exit(0) if scalar @ARGV == 0;

my $run_path;
my $output_json_path;
my $quality_metrics_csv_path;
my $demultiplex_stats_csv_path;
my $run_info_xml_path;
my $sample_sheet_csv_path;

my %opt = ();
GetOptions (
    "runDir=s"      => \$run_path,
    "qualityCsv=s"  => \$quality_metrics_csv_path,
    "demuxCsv=s"    => \$demultiplex_stats_csv_path,
    "runInfoXml=s"  => \$run_info_xml_path,
    "sampleSheet=s" => \$sample_sheet_csv_path,
    "outputJson=s"  => \$output_json_path,
    "sep=s"         => \$OUT_SEP,
    "decimals=i"    => \$ROUND_DECIMALS,
    "yieldFactor=i" => \$opt{yield_factor},
    "noQc"          => \$opt{no_qc},
    "debug"         => \$opt{debug},
    "help|h"        => \$opt{help},
) or die "[ERROR] Issue in command line arguments\n";
print $HELP and exit(0) if $opt{help};

if (defined $run_path){
    -d $run_path || die "[ERROR] Provided run dir does not exist ($run_path)\n";
    $run_info_xml_path = "$run_path/$REL_RXML_PATH";
    $sample_sheet_csv_path = "$run_path/$REL_SSHT_PATH";
    $quality_metrics_csv_path = "$run_path/$REL_QCSV_PATH";
    $demultiplex_stats_csv_path = "$run_path/$REL_DCSV_PATH";
}

fileExistsOrDie($run_info_xml_path);
fileExistsOrDie($sample_sheet_csv_path);
fileExistsOrDie($quality_metrics_csv_path);
fileExistsOrDie($demultiplex_stats_csv_path);

my $ssht_info = readSampleSheet($sample_sheet_csv_path);
my $platform = determinePlatformByString($ssht_info->{platform}, \@KNOWN_PLATFORMS);
my $settings = $SETTINGS_PER_PLATFORM{$platform};
$settings->{yield_factor} = $opt{yield_factor} if defined $opt{yield_factor};
my $YIELD_FACTOR = $settings->{yield_factor};

my $run_info = readRunInfoXml($run_info_xml_path);
my $quality_info = readBclConvertQualityMetricsCsv($quality_metrics_csv_path, $run_info);
my $demux_info = readBclConvertDemuxCsv($demultiplex_stats_csv_path, $run_info);
my $store = aggregateBclconvertInfo($quality_info, $demux_info, $run_info);
$store->{stats}{platform} = $platform;
$store->{stats}{run_id} = $run_info->{run_id};
addSamplesheetInfo($store, $ssht_info);
addPrintFields($store, $run_info);
addGeneralStats($store, $run_info);

performQC($store, $settings, $run_info) unless $opt{no_qc};

printTable($store, \@OUT_FIELDS);
printJson($store, $output_json_path) if defined $output_json_path;

if ($opt{debug}){
    printJson($ssht_info, 'tmp_debug_ssheet.json');
    printJson($quality_info, "tmp_debug_quality.json");
    printJson($demux_info, "tmp_debug_demux.json");
    printJson($store, 'tmp_debug_store.json');
}

sub aggregateBclconvertInfo{
    my ($qual, $demux, $runinfo) = @_;
    my %out = ();
    my $fid = $runinfo->{flowcell_id};

    $out{flowcell}{$fid}{id} = $fid;
    $out{flowcell}{$fid}{name} = "FlowcellNameUnknown";

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
    return \%out;
}

sub fileExistsOrDie{
    my ($file) = @_;
    my $result;
    if ($file =~ "gs://"){
        $result = not system("gsutil stat $file > /dev/null") ;
    }else{
        $result = -f $file;
    }
    $result || die "[ERROR] File does not exist ($file)\n";
}

sub addSamplesheetInfo{
    my ($data, $samplesheet) = @_;

    my $run_name = $samplesheet->{run_name}; # eg NO22-0001
    my @fcids = keys %{$data->{flowcell}};
    scalar @fcids == 1 || die "[ERROR] There should be exactly one flowcell in data\n";
    my $fcid = $fcids[0];
    $data->{flowcell}{$fcid}{name} = $run_name;
    $data->{flowcell}{$fcid}{name_print} = $run_name;
    $data->{stats}{run_name} = $run_name;
    $data->{stats}{override_cycles} = $samplesheet->{override_cycles};
    $data->{stats}{submissions} = {};

    # Add sample metadata
    my %submissions = ();
    my $samples = $data->{samples};
    foreach my $barcode (keys %$samples){
        my $sample = $samples->{$barcode};
        my $submission = $samplesheet->{samples}{$barcode}{Sample_Project} || NA;
        $submissions{$submission} = 1;
        $sample->{submission} = $submission;
        $sample->{name} = $samplesheet->{samples}{$barcode}{Sample_Name} || NA;
        $sample->{description} = $samplesheet->{samples}{$barcode}{Description} || NA;;
    }
    $data->{stats}{submissions} = [sort keys %submissions];
}

sub readRunInfoXml{
    my ($run_info_xml) = @_;
    my $content = readXml($run_info_xml);

    my $flowcell_id = $content->{Run}{Flowcell};
    my $run_id = $content->{Run}{Id};
    my @cycle_counts = map($_->{NumCycles}, @{$content->{Run}{Reads}{Read}});
    my @non_index_cycle_counts = map($_->{IsIndexedRead} eq "N" ? $_->{NumCycles} : (), @{$content->{Run}{Reads}{Read}});

    my $cycle_string = join("|", @cycle_counts);
    my $total_non_index_cycle_count = sum(@non_index_cycle_counts);
    my %info = (
        "cycle_string" => $cycle_string,
        "total_non_index_cycle_count" => $total_non_index_cycle_count,
        "flowcell_id" => $flowcell_id,
        "run_id" => $run_id
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
    my ($info, $qc_limits, $run_info) = @_;
   
    my $stats = $info->{stats};
    my $samples = $info->{samples};
    my $lanes = $info->{lanes};
    my $fails = 0;
    my $identifier = $stats->{identifier};

    # Flowcell checks
    my $undet = $stats->{undet_perc};
    my $max_undet = $qc_limits->{max_undetermined};
    if ($undet > $max_undet){
        warn "## WARNING Percentage undetermined ($undet) too high (max=$max_undet)\n";
        $fails += 1;
    }
    my $fid = $run_info->{flowcell_id};
    my $q30 = $info->{flowcell}{$fid}{q30};
    my $min_flowcell_q30 = $qc_limits->{min_flowcell_q30};
    if ($q30 < $min_flowcell_q30){
        warn "## WARNING q30 ($q30) too low (min=$min_flowcell_q30)\n";
        $fails += 1;
    }
        
    # Sample checks
    $fails += checkObjectField($samples, 'yield', $qc_limits->{min_sample_yield});
    checkObjectField($lanes, 'q30',   $qc_limits->{min_flowcell_q30}); # lane check only prints warning
    
    # Conclusion
    my $final_qc_result = "NoQcResult";
    if ($fails == 0){
        $final_qc_result = "PASS";
        say "## FINAL QC RESULT: OK";
    }
    else{
        $final_qc_result = "FAIL";
        warn "## WARNING Some checks failed, inspect before proceeding (for $identifier)\n";
        say "## FINAL QC RESULT: FAIL ($fails failures for $identifier)";
    }
    $stats->{flowcell_qc} = $final_qc_result;
}

sub checkObjectField{
    my ($objects, $field, $min) = @_;
    my $fails = 0;
    foreach my $obj_key (sort {$objects->{$b}{name} cmp $objects->{$a}{name}} keys %$objects){
        my $obj = $objects->{$obj_key};
        my $name = $obj->{name};
        next if $name eq UNDETERMINED;
        next if $name =~ /^VirtualSample/;
        my $value = 0;
        $value = $obj->{$field} if exists $obj->{$field};
        if ($value < $min){
            warn "## WARNING $field for $name too low: $value < $min\n";
            $fails += 1;
        }
    }
    return $fails;
}

sub addGeneralStats{
    my ($store, $run_info) = @_;
    my $fid = $run_info->{flowcell_id};
    my $cycle_string = $run_info->{cycle_string};

    my $undet_id = UNDETERMINED;
    my $undet_perc = percentage($store->{undetermined}{$undet_id}{yield}, $store->{flowcell}{$fid}{yield});
    my $run_overview_yield_factor = 1e6; # always report the run info in MBase
    $store->{stats}{run_overview_string} = sprintf "%s\t%s\t%s\t%s\t%s\t%s",
        round($store->{flowcell}{$fid}{yield}, 0, $run_overview_yield_factor),
        round($store->{undetermined}{$undet_id}{yield}, 0, $run_overview_yield_factor),
        $store->{flowcell}{$fid}{q30_print},
        NA, # placeholder for "passing filter percentage" was available in bcl2fastq output but not in bclconvert output
        $cycle_string,
        round($undet_perc,1) . '%';

    $store->{stats}{undet_perc} = $undet_perc;
    $store->{stats}{lane_count} = scalar(keys %{$store->{lanes}});
    $store->{stats}{samp_count} = scalar(keys %{$store->{samples}});
    $store->{stats}{identifier} = join("_", keys %{$store->{flowcell}});
    $store->{stats}{cycle_string} = $cycle_string;
}

sub addPrintFields{
    my ($store, $run_info) = @_;
    my $fid = $run_info->{flowcell_id};

    foreach my $type (keys %$store){
        next if $type eq 'stats';
        foreach my $id (keys %{$store->{$type}}){
            my $obj = $store->{$type}{$id};

            $obj->{q30} = percentage($obj->{yield_q30}, $obj->{yield});
            $obj->{yld_p} = percentage($obj->{yield}, $store->{flowcell}{$fid}{yield});

            $obj->{q30_print} = round($obj->{q30}, $ROUND_DECIMALS, 1);
            $obj->{yld_print} = round($obj->{yield}, 0, $YIELD_FACTOR);
            $obj->{yld_p_print} = round($obj->{yld_p}, $ROUND_DECIMALS, 1);

            $obj->{flowcell_print} = $fid;
            $obj->{id_print} = $id;
            $obj->{name_print} = $obj->{name};
            $obj->{index1_print} = $obj->{index1};
            $obj->{index2_print} = $obj->{index2};

            if ($type =~ /samples|lanes|flowcell/) {
                $obj->{total_reads} = $obj->{mm0} + $obj->{mm1} + $obj->{mm2};
                $obj->{mm0_print} = 0;
                $obj->{mm1_print} = 0;
                $obj->{mm2_print} = 0;
                if ($obj->{total_reads} != 0) {
                    $obj->{mm0_print} = round(percentage($obj->{mm0}, $obj->{total_reads}), $ROUND_DECIMALS);
                    $obj->{mm1_print} = round(percentage($obj->{mm1}, $obj->{total_reads}), $ROUND_DECIMALS);
                    $obj->{mm2_print} = round(percentage($obj->{mm2}, $obj->{total_reads}), $ROUND_DECIMALS);
                }
            }

            if ($type =~ /samples/){
                $obj->{submission_print} = $obj->{submission};
            }
        }
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
    my @submissions = sort @{$info->{stats}{submissions}};
    map($_ =~ s/HMFreg//, @submissions);
    my $submissions_string = join(',', @submissions);

    printStandardHeaderLines($info);

    say sprintf "## RunOverview INFO: %s\t%s\t%s\t%s\t%s",
        $info->{stats}{run_name},
        $info->{stats}{run_id},
        $submissions_string,
        $info->{stats}{run_overview_string},
        $info->{stats}{flowcell_qc};
    
    say "#".join($OUT_SEP, "level", @$fields);
    printTableForLevelSortedByName($info->{flowcell}, $fields, 'RUN');
    printTableForLevelSortedByName($info->{lanes}, $fields, 'LANE');
    printTableForLevelSortedByName($info->{samples}, $fields, 'SAMPLE');
    printTableForLevelSortedByName($info->{reads}, $fields, 'READ');
    printTableForLevelSortedByName($info->{undetermined}, $fields, 'UNDET');
}

sub printStandardHeaderLines{
    my ($data) = @_;
    say sprintf '## Settings INFO: Platform=%s|OverrideCycles=%s|YieldFactor=%s|RoundDecimals=%s',
        $data->{stats}{platform},
        $data->{stats}{override_cycles},
        commify($YIELD_FACTOR),
        commify($ROUND_DECIMALS);

    say sprintf '## Flowcell INFO: %s (%s, %d lanes, %d samples, %s cycles)',
        $data->{stats}{run_id},
        $data->{stats}{run_name},
        $data->{stats}{lane_count},
        $data->{stats}{samp_count},
        $data->{stats}{cycle_string};
}

sub printTableForLevelSortedByName{
    my ($info, $fields, $level) = @_;
    foreach my $id (sort {$info->{$b}{name} cmp $info->{$a}{name}} keys %$info){
        my @output = map($info->{$id}{$_."_print"} || NA, @$fields);
        unshift @output, $level if defined $level;
        say join($OUT_SEP, @output);
    }
}

sub printTableForLevelSortedByYield{
    my ($info, $fields, $level) = @_;
    foreach my $id (sort {$info->{$b}{yld_print} <=> $info->{$a}{yld_print}} keys %$info){
        my @output = map($info->{$id}{$_."_print"} || NA, @$fields);
        unshift @output, $level if defined $level;
        say join($OUT_SEP, @output);
    }
}

sub readSampleSheet{
    my ($csv_file) = @_;
    
    # SampleSheet file has windows returns
    my $return_str = $/;
    $/ = "\r\n";
    
    my %output;
    $output{samples} = {};
    $output{run_name} = 'NO_RUNNAME_FROM_SAMPLESHEET';
    $output{override_cycles} = 'NO_OVERRIDECYCLES_FROM_SAMPLESHEET';
    $output{platform} = 'NO_PLATFORM_FROM_SAMPLESHEET';
    my @header;
    
    if (! -e $csv_file){
        say "## WARNING skipping SampleSheet read: file not found ($csv_file)";
        return(\%output);
    }
    
    open FILE, "<", $csv_file or die "Couldn't open file ($csv_file): $!";
    while (<FILE>) {
        chomp($_);
        next if $_ =~ /^[\[\,]/;
        next if $_ eq "";
        my @fields = split(",", $_);

        if ($fields[0] =~ /Experiment(.)*Name/){
            my $run_name = $fields[1] || 'NA';
            $output{run_name} = $run_name;
            if ($run_name =~ m/^NO\d{2}-/){
                $output{platform} = $NOVASEQ;
            }elsif ($run_name =~ m/^NS\d{2}-/){
                $output{platform} = $NEXTSEQ;
            }elsif ($run_name =~ m/^IS\d{2}-/){
                $output{platform} = $ISEQ;
            }elsif ($run_name =~ m/^X\d{2}-/){
                $output{platform} = $HISEQ;
            }else{
                warn "Unable to determine platform from experiment name field ($run_name)"
            }
        }
        elsif ($fields[0] =~ /OverrideCycles/){
            $output{override_cycles} = $fields[1] || 'NA';
        }
        elsif ($_ =~ m/Sample_ID/){
            @header = @fields;
        }
        elsif (@header){
            my %tmp = ();
            $tmp{$_} = shift @fields foreach @header;
            my $sample_id_column = 'Sample_ID';
            if (defined $tmp{$sample_id_column} and $tmp{$sample_id_column} ne ''){
                my $sample_id = $tmp{$sample_id_column};
                $output{samples}{$sample_id} = \%tmp;
            }
        }
    }
    close FILE;
    $/ = $return_str; # reset return string
    return(\%output);
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
        my $store_location;
        if ($barcode eq UNDETERMINED){
            $store_location = \%{$quality_store{undetermined}{$barcode}};
        }else{
            $store_location = \%{$quality_store{samples}{$barcode}};
            $store_location->{index1} = $line->{index};
            $store_location->{index2} = $line->{index2};
        }
        foreach my $value_key (keys %value_keys){
            my $store_key = $value_keys{$value_key};
            $quality_store{flowcell}{$fcid}{$store_key} += $line->{$value_key};
            $quality_store{lanes}{$lane}{$store_key} += $line->{$value_key};
            $quality_store{reads}{$read}{$store_key} += $line->{$value_key};
            $store_location->{$store_key} += $line->{$value_key};
        }
        $store_location->{name} = $barcode;
        $quality_store{lanes}{$lane}{name} = $lane;
        $quality_store{reads}{$read}{name} = $read;
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
            if ($barcode eq UNDETERMINED){
                $demux_store{undetermined}{$barcode}{$store_key} += $line->{$value_key};
            }else {
                $demux_store{flowcell}{$fcid}{$store_key} += $line->{$value_key};
                $demux_store{lanes}{$lane}{$store_key} += $line->{$value_key};
                $demux_store{samples}{$barcode}{$store_key} += $line->{$value_key};
            }
        }

    }
    return \%demux_store;
}

sub readXml{
    my ($file) = @_;
    my $entire_file = "";
    if ($file =~ "gs://"){
        $entire_file = `gsutil cat $file`;
    }else{
        $entire_file = `cat $file`;
    }
    my $obj = XMLin($entire_file);
    return($obj);
}

sub readJson{
    my ($json_file) = @_;
    my $json_txt = read_file($json_file);
    my $json_obj = decode_json($json_txt);
    return($json_obj);
}

sub percentage{
    my ($value, $total) = @_;

    if (not defined($value) or not defined($total)) {
        return 0;
    }
    elsif ($value < 0 or $total < 0) {
        die "[ERROR] Cannot calculate percentage if either value ($value) or total ($total) is < 0\n";
    }
    elsif ($value > $total){
        die "[ERROR] value ($value) should never be higher than total ($total)\n";
    }
    elsif ($total == 0 and $value == 0) {
        return 0;
    }
    else {
        return $value*100/$total;
    }
}

sub round{
    my ($number, $decimal, $factor) = @_;
    if (not defined $number){
        return NA;
    }
    $decimal = 0 unless defined $decimal;
    $factor = 1 unless defined $factor;
    my $rounded = sprintf("%.".$decimal."f", $number/$factor);
    return($rounded);
}

sub commify{
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

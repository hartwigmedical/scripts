#!/usr/bin/perl
use strict;
use warnings;
use 5.010.000;

use constant SEMICOLON => ';';
use constant NACHAR    => 'NA';
use constant COMMA     => ',';

## -----
## Description:
##  Merges COSMIC VCF records by chr/pos/ref/alt to create records 
##  with all COSM ids per unique variant "collapsed" into one
## Example new record
##  1 906267 COSM4623802;COSM4623801 C T . . COSM2ENST=${infoFields}
## -----

my $bgzip = "/data/common/tools/tabix_v0.2.6/bgzip";
my $tabix = "/data/common/tools/tabix_v0.2.6/tabix";

my $DELIM = '|';
my @COSMIC_VCF_FIELDS = qw( GENE CDS AA CNT );
my $OUT_FORMAT = join($DELIM, 'ID', @COSMIC_VCF_FIELDS);
my $OUT_INFO_KEY = "COSM2ENST";
my $OUT_INFO_DSC = sprintf("##INFO=<ID=%s,Number=.,Type=String,Description=\"COSMIC to ENST mapping %s\">", $OUT_INFO_KEY, $OUT_FORMAT);

## Script arguments
my $csm_mut_vcf = shift; # eg CosmicCodingMuts_v84.vcf.gz
my $csm_tra_tsv = shift; # eg CosmicTranscripts_v84.tsv.gz
my $output_file = shift; # optional

## Input checks
my $script_arguments = "<cosmic-mut-vcf.gz> <cosmic-transcript-tsv.gz> [<optional-output-file-path>]";
die "[EXIT] Run with $0 $script_arguments\n" unless $csm_mut_vcf and $csm_tra_tsv;
die "[EXIT] Input vcf must be gzipped VCF\n" unless $csm_mut_vcf =~ /\.vcf\.gz$/;
die "[EXIT] Input tsv must be gzipped TSV\n" unless $csm_tra_tsv =~ /\.tsv\.gz$/;
die "[EXIT] File not found ($csm_mut_vcf)\n" unless -f $csm_mut_vcf;
die "[EXIT] File not found ($csm_tra_tsv)\n" unless -f $csm_tra_tsv;
($output_file = $csm_mut_vcf ) =~ s/\.vcf\.gz$/_collapsed\.vcf/ unless $output_file;
die "[EXIT] Output VCF $output_file(.gz) already exists\n" if -f $output_file or -f "$output_file.gz";

## -----
## MAIN
say "[INFO] Parsing transcript input TSV";
my $gene2transcript = parseTranscripts($csm_tra_tsv);

say "[INFO] Opening tmp output file";
open OUT, ">", $output_file or die "[EXIT] Unable to open output filehandle: $!\n";
  printAdjustedVcfHeader($csm_mut_vcf, *OUT, $OUT_INFO_DSC);
  printCollapsedVariants($csm_mut_vcf, *OUT, \@COSMIC_VCF_FIELDS, $gene2transcript);
close OUT;

say "[INFO] Running bgzip/tabix on tmp output file";
bgzipAndTabixVcf($output_file, $bgzip, $tabix);

say "[INFO] DONE";
## /MAIN
## -----

## -----
## SUBROUTINES
## -----
sub printCollapsedVariants{
    my ($vcf, $output_filehandle, $ann_fields, $gene2tran) = @_;
    my $prev_chr = '';
    my $prev_pos = '';
    my %variants = ();
    
    ## Assume VCF is sorted by chr + pos and gzip compressed
    open(my $fh, "zcat $vcf |") or die $!;
    
    while ( <$fh> ){
      
        chomp;
        next if $_ =~ /^#/;
        next if $_ eq "";
        
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = split( "\t", $_ );
        my @annotations = ();
        
        foreach my $field ( @$ann_fields ){
            my $ann = getFieldFromInfoString($info, $field);
            $ann = addTranscript( $ann, $gene2tran ) if ( $field eq "GENE" );
            push( @annotations, $ann );
        }
        
        my $isFirstRecord  = $prev_chr eq "";
        my $isPrevLocation = ($chr eq $prev_chr and $pos eq $prev_pos);
        
        my $info_string = join( $DELIM, $id, @annotations );
        my %one_variant = ( 'ann' => $info_string, 'id'  => $id );
        
        ## Location check
        if ( $isFirstRecord ){ 
            push( @{$variants{$ref}{$alt}}, \%one_variant );
        }
        elsif( $isPrevLocation ){ 
            push( @{$variants{$ref}{$alt}}, \%one_variant );
        }
        else{ # is new location
            collapseByRefAltAndPrint( \%variants, $prev_chr, $prev_pos, $output_filehandle );
            %variants = ();
            push( @{$variants{$ref}{$alt}}, \%one_variant );
        }
        ($prev_chr,$prev_pos) = ($chr,$pos);
      
    }

    ## Print final variant
    collapseByRefAltAndPrint( \%variants, $prev_chr, $prev_pos, $output_filehandle );
    
    close $fh;
}

sub addTranscript{
    my ($gene_name, $gene2tran) = @_;
    
    if ( $gene_name =~ /_/ ) { # already has transcript id
        return( $gene_name );
    }
    elsif ( exists $gene2tran->{ $gene_name } ){
        my $tran_count = scalar @{$gene2tran->{ $gene_name }};
        die "[EXIT] Transcripts found for \"$gene_name\" but count != 1 ($tran_count)" if ( $tran_count != 1 );
        return $gene_name . '_' . $gene2tran->{ $gene_name }[0];
    }
    else{
        return( $gene_name );
        warn "[WARN] Gene $gene_name has no transcript!\n";
    }
}

sub collapseByRefAltAndPrint{
    my ($variants, $chr, $pos, $output_filehandle) = @_;
    
    foreach my $ref ( sort keys %$variants ){
        my $alts = $variants->{ $ref };
        foreach my $alt ( keys %$alts ){
            my $ref_alt_variants = $alts->{ $alt };
            my @ids = ();
            my @ann = ();
            foreach my $var ( @$ref_alt_variants ){
                push( @ids, $var->{'id'} );
                push( @ann, $var->{'ann'} );
            }
            my $ids = join( SEMICOLON, @ids );
            my $ann = join( COMMA, @ann );
            say $output_filehandle join("\t", $chr, $pos, $ids, $ref, $alt, '.', '.', join( "=", $OUT_INFO_KEY, $ann));
        }
    }
}

sub parseTranscripts{
    my ($tsv) = @_;
    my %gene2transcript = ();
    
    ## Assume 3 columns: GeneID GeneName TranscriptID
    open(my $fh, "zcat $tsv |") or die $!;
    
    while ( <$fh> ){
        chomp;
        next if $_ =~ /^#/;
        next if $_ eq "";
        next if $_ =~ /Gene ID/;
        my ($gene_id,$gene_name,$transcript_id) = split( "\t", $_ );
        if ( exists %gene2transcript{ $gene_name } ){
            die "[EXIT] Should not occur: gene \"$gene_name\" key already present\n";
        }
        push( @{%gene2transcript{ $gene_name }}, $transcript_id );
    }
    return(\%gene2transcript);
}

sub printAdjustedVcfHeader{
    my ($vcf_file, $output_fh, $additional_header_lines) = @_;
    open(my $header_fh, "zcat $vcf_file | grep '^#' |") or die $!;
    while ( <$header_fh> ){
        chomp;
        if ( $_ =~ /^#CHROM/ ){
            say $output_fh $additional_header_lines;
        }
        say $output_fh $_;
    }
    close $header_fh;
}

sub getFieldFromInfoString{
    my ($info_string, $select_field) = @_;
    my %out_map = ();
    my @fields = split( SEMICOLON, $info_string);
    foreach my $field ( @fields ){
        my ($name,$content) = split( "=", $field );
        if ( $name eq $select_field ){
            chomp($content);
            return $content;
        }
    }
    return NACHAR;
}

sub bgzipAndTabixVcf{
    my ($file, $bgzip, $tabix) = @_;
    if ( not -f $bgzip or not -f $tabix ){
        say "[INFO] Skipping bgzip/tabix because tools not available";
        return;
    }
    system( "$bgzip -f $file" );
    system( "$tabix -p vcf ${file}.gz" );
    say "[INFO] Inspect compressed output with: zcat $file.gz | less";
}

#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use File::Slurp;
use JSON::XS;
use Number::Format qw(:subs);
use 5.010.000;

## -----
## Global variables
## -----
my $DATETIME = localtime;
my $SCRIPT = basename $0;
my $JSONS_HOME = '/data/ops/api/prod/database';
my $JSONS_HOME_ACC = '/data/ops/api/acc/database';

my %OUT_FIELDS_PER_TYPE = (
  'samples'   => [ 'submission', 'barcode', 'q30', 'yld_req', 'yld', 'status', 'name', 'id' ],
  'flowcells' => [ 'name', 'createTime', 'flowcell_id', 'sequencer', 'q30', 'yld', 'undet_rds', 'status', 'undet_rds_p_pass', 'id' ],
  'fastq'     => [ 'sample_id', 'name_r1', 'qc_pass', 'q30', 'yld', 'bucket', 'id' ],
  'sets'      => [ 'name', 'entity', 'ref_sample', 'tumor_sample', 'id' ],
  'runs'      => [ 'name', 'entity', 'ref_sample', 'tumor_sample', 'bucket', 'status', 'priority', "pipeline", "ini", 'id' ],
  'entities'  => [ 'name', 'bucket', 'id' ],
  'inis'      => [ 'name', 'id' ],
  'shares'    => [ 'entity_id', 'entity', 'set_id', 'set', 'start_time', 'end_time', 'filter', 'id' ],
  'stacks'    => [ 'name', 'revision', 'enabled', 'id' ],
);
my $available_types = join( ", ", sort keys %OUT_FIELDS_PER_TYPE );

my $delim = "\t";
my $type;
my $output_as_json;
my @filters = ();
my $must_match_exact;
my $use_acceptance;
my $no_number_formatting;

my $HELP =<<HELP;

  Description
    Parses DB/API info from api jsons ($JSONS_HOME)
    (updated regularly by cronjob)

    You can manually update the jsons to latest by:
      - executing "update_api_db"
  
  Available tables/types to query:
     $available_types

  Filtering output (see "-filter" input option):
    - You can filter on any field within the selected table/type
    - Use the "-json" param to see all keys of objects
    - Negate any filter by using "!=" instead of "="
  
  Usage
    $SCRIPT -type <type> [options]

  Usage examples
    $SCRIPT -type samples | head
    $SCRIPT -type samples | grep CPCT02 | head
    $SCRIPT -type samples -json | jq '.' | head
    $SCRIPT -type samples -filter "barcode=FR13825534"
    $SCRIPT -type runs -filter "status!=Validated"
    $SCRIPT -type flowcells -filter "status=Pending"
    $SCRIPT -type flowcells -filter "status=Pending|Sequencing"

  Input options
    -filter   <s>  (filter string with = or !=)
    -exact         (match filters exact instead of regex)
    -acc           (use acceptance api instead of prod)

  Output options 
    -delim    <s>  (output delim for tabular output)
    -json          (output as json instead of tabular)
    -no_format     (do not round and format any numbers)
    
HELP
print $HELP and exit(0) if scalar @ARGV == 0 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help';

## -----
## Gather input
## -----
GetOptions (
    "type=s"    => \$type,
    "filter=s@" => \@filters,
    "delim=s"   => \$delim,
    "acc"       => \$use_acceptance,
    "json"      => \$output_as_json,
    "exact"     => \$must_match_exact,
    "no_format" => \$no_number_formatting,
) or die "Error in command line arguments\n";
warn "[ERROR] No type given?" and exit(0) unless $type;
warn "[ERROR] Type ($type) not supported" and exit(0) unless exists $OUT_FIELDS_PER_TYPE{ $type };

## -----
## MAIN
## -----
if ( $use_acceptance ){
   $JSONS_HOME = $JSONS_HOME_ACC;
}
my $objects = readJson( $type );
my $filtered_objects = filterObjects( $objects, \@filters );
my $out_fields = $OUT_FIELDS_PER_TYPE{ $type };

if ( $output_as_json ){
    printOutputAsJson( $filtered_objects );
}else{
    printOutput( $filtered_objects, $out_fields );
}

## -----
## /MAIN
## -----

## generic json reader for all types
sub readJson{
    my ($type) = @_;
    my $json_file = "$JSONS_HOME/$type.json";
    my $json_obj = jsonFileToObject($json_file);
    
    ## optimize some fields for viewing
    foreach my $obj ( @$json_obj ){
        $obj->{ 'id' } = $obj->{ 'id' } if defined $obj->{ 'id' };
        unless ( $no_number_formatting ){
            $obj->{ 'q30' } = sprintf( "%.1f", $obj->{ 'q30' } ) if defined $obj->{ 'q30' };
            $obj->{ 'yld' } = format_number( $obj->{ 'yld' } / 1000000, 0 ) if defined $obj->{ 'yld' };
            $obj->{ 'yld_req' } = format_number( $obj->{ 'yld_req' } / 1000000, 0 ) if defined $obj->{ 'yld_req' };
            $obj->{ 'undet_rds' } = format_number( $obj->{ 'undet_rds' } / 1000000, 0 ) if defined $obj->{ 'undet_rds' };
        }
    }
    
    ## Add set info for runs (set info is situated one level deeper)
    if ( $type eq 'runs' ){
        foreach my $obj ( @$json_obj ){
            my @keys_to_move = qw( name ref_sample tumor_sample entity_id );
            foreach my $key ( @keys_to_move ){
                $obj->{ $key } = $obj->{ 'set' }{ $key } if defined $obj->{ 'set' }{ $key };
            }
        }
    }
    
    ## Add extra information for sets/runs from other jsons
    if ( $type eq 'sets' or $type eq 'runs' or $type eq 'shares' ){
        my $entities = jsonFileToObject( "$JSONS_HOME/entities.json" );
        my $inis = jsonFileToObject( "$JSONS_HOME/inis.json" );
        my $stacks = jsonFileToObject( "$JSONS_HOME/stacks.json" );
        foreach my $obj ( @$json_obj ){
            $obj->{ 'entity' } = getFieldValueById( $entities, 'name', $obj->{ 'entity_id' } ) if defined $obj->{ 'entity_id' };
            $obj->{ 'ini' } = getFieldValueById( $inis, 'name', $obj->{ 'ini_id' } ) if defined $obj->{ 'ini_id' };
            $obj->{ 'pipeline' } = getFieldValueById( $stacks, 'revision', $obj->{ 'stack_id' } ) if defined $obj->{ 'stack_id' };
        }
    }

    if ( $type eq 'shares' ){
        my $sets = jsonFileToObject( "$JSONS_HOME/sets.json" );
        foreach my $obj ( @$json_obj ){
            $obj->{ 'set' } = getFieldValueById( $sets, 'name', $obj->{ 'set_id' } ) if defined $obj->{ 'set_id' };
        }
    }

    return( $json_obj );
}

sub jsonFileToObject{
    my ($json_file_path) = @_;
    my $json_txt = read_file( $json_file_path );
    my $json_obj = decode_json( $json_txt );
    return( $json_obj );
}

sub getFieldValueById{
    my ($search_objects, $request_field, $search_id) = @_;
    my $return = "NA";
    foreach my $obj ( @$search_objects ){
        if ( $obj->{ 'id' } == $search_id ){
            $return = $obj->{ $request_field } if defined $obj->{ $request_field };
        }
    }
    
    return( $return );
}

sub filterObjects{
    my ($objects, $filters) = @_;
    
    my @out = ();
    my %filter_counts = ();
    
    foreach my $obj ( @$objects ){
        my $do_skip_object = applyFiltersOnObject( $obj, $filters, \%filter_counts );
        push( @out, $obj ) unless $do_skip_object;        
    }
    
    foreach my $filter ( keys %filter_counts ){
        my $count = $filter_counts{ $filter };
        #say "## FILTER: $count filtered away by filter \"$filter\"";
    }
    
    return \@out;
}

sub applyFiltersOnObject{
    my ($object, $filters, $counts) = @_;
        
    foreach my $filter ( @$filters ){
        my $negate_search;
        my $field;
        my $match_str;

        if ($filter =~ /\!\=/){
            $negate_search = 1;
            ($field, $match_str) = split( "!=", $filter, 2 );
        }
        elsif ($filter =~ /=/){
            $negate_search = 0;
            ($field, $match_str) = split( "=", $filter, 2 );
        }
        else{
            die "[ERROR] Incorrect filter format ($filter). Requires a = or != in string.\n"
        }

        if ( not exists $object->{ $field } ){
            $counts->{ $field.'=KeyNotExists' }++;
            return(1);
        }
        elsif ( not defined $object->{ $field } ){
            $counts->{ $field.'=KeyNotDefined' }++;
            return(1);
        }
        else{
            my $exact_match = $object->{ $field } eq $match_str;
            my $regex_match = $object->{ $field } =~ m/$match_str/i;

            if ( $negate_search ){
                $exact_match = !$exact_match;
                $regex_match = !$regex_match;
            }
            
            if ( $must_match_exact ){
                if ( not $exact_match ){
                    $counts->{ $filter }++;
                    return(1);
                }
            }
            elsif( not $regex_match ){
                $counts->{ $filter }++;
                return(1);
            }
        }
    }
    
    ## all filters applied and still here so return OK response
    return 0;
}

## print output
sub printOutputAsJson{
    my ($objects) = @_;
    my $json_text = encode_json( $objects );
    say $json_text;
}

## print output
sub printOutput{
    my ($objects, $out_fields) = @_;
    my $object_count = scalar @$objects;
    say "## Script: $SCRIPT";
    say "## DateTime: $DATETIME";
    say '## TotalCount: '.$object_count;
    say '#'.join( $delim, @$out_fields );
    foreach my $obj ( @$objects ){
        my @out_values = ();
        foreach my $field ( @$out_fields ){
            if ( defined $obj->{$field} ){
                push( @out_values, $obj->{$field} );
            }
            else{
                push( @out_values, "NA" );
            }
        }
        say join( $delim, @out_values );
    }
}

## trims whitespace on both ends of a string
sub trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

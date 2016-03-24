#!/usr/bin/perl -w
use strict;

open IN, $ARGV[0]; #list of jobids
my $command = $ARGV[1];
while (<IN>) {
    chomp();
    system "$command $_";
}
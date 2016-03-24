#!/usr/bin/perl -w
use strict;

open IN, $ARGV[0]; # submission log

while (my $line = <IN>) {
    chomp($line);
    if($line =~  /Your job (\d{4,6})/) {
	my $jid = $1;
	#print $line,"\ $jid\n";
	system "qdel $jid";
    }
}
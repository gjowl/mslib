#!/usr/bin/perl

my $infile = shift;

open INFILE, $infile or die "ERROR: DID NOT OPEN $infile\n";

while (my $line = <INFILE>) {
	my @matches = $line =~ m/\[(.+?)\]/;
	foreach $match(@matches) {
		chomp $match;
		print("$match\t");
	}
	print("\n");
}

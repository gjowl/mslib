#!/usr/bin/perl
use strict;
use Chart::Gnuplot;
use File::Basename;
use List::Util qw(min max);

# Parse text file
my $inputfile = shift;
open INPUT, $inputfile or die "ERROR: Could not open INPUT!\n";
my $outputDir = shift;
#my $axialRotateA = shift;

my ($fname, $fpath, $fsuffix) = fileparse($inputfile);

my @data;

my $header = <INPUT>;
chomp($header);
my @labels = split("\t", $header);

while (my $line = <INPUT>) {
	chomp($line);
	my @vals = split ("\t", $line);
	for (my $i = 0; $i < @vals; $i++) {
		push (@{$data[$i]}, $vals[$i]);
	}
	
}

my $min = 999999999999;
my $max = -999999999999;

for (my $i = 1; $i < @data; $i++) {
	my $posmin = min(@{$data[$i]});
	my $posmax = max(@{$data[$i]});
	if ($posmin <= $min) {
		$min = $posmin;
		print "new min $min\n";
	}
	if ($posmax >= $max) {
		$max = $posmax;
		print "new max $max\n";
	}
}

print "$min $max\n";


#Chart object
for (my $i = 1; $i < @data; $i++) {
	my $chart = Chart::Gnuplot->new (
		output	=> "$outputDir/$fname-$labels[$i].png",
		title	=> "$labels[$i]",
		#view	=> 'map',
		#size	=> 'square',
		yrange	=> [$min,$max],
	);
	
	#Data set object
	my $dataSet = Chart::Gnuplot::DataSet->new (
		xdata	=> \@{$data[0]},
		ydata	=> \@{$data[$i]},
		style	=> 'lines',
	);
	
	#plot chart
	$chart->plot2d($dataSet);
}

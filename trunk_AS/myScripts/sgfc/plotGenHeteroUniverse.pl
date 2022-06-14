#!/usr/bin/perl
use warnings;
use strict;
use Chart::Gnuplot;

# Parse text file
open INPUT, shift or die "ERROR: Could not open INPUT!\n";
my $output = shift;
my $Za = shift;
my $Zb = shift;
#my $axialRotateA = shift;


my (@theta, @omegaB, @hBond);

while (my $line = <INPUT>) {
	if ($line =~ /^[+-]?\d+/) {
		chomp($line);
		my @data = split ("\t", $line);
		if ($data[3] == $Za && $data[4] == $Zb) {
			my $omegaB = $data[1];
			my $theta = $data[2];
			my $hBond = $data[12];
			#Add blank line if x value changes--required to change scan in pm3d plot
			if ($omegaB != $omegaB[-1]) {
				push (@omegaB, "");
				push (@theta, "");
				push (@hBond, "");
			}
	
			#my $omegaPrime = (10 * $data[1] / 9) + (200 * $data[2] / 27);
			my $omegaPrime = $omegaB;
			push (@omegaB, $omegaPrime);
			push (@theta, $theta);
			push (@hBond, $hBond);
		}
	}
	
}

# CODE FOR PLOTTING MODEL FILES FROM SAB'S ORIGINAL GENHOMOUNIVERSE PROGRAM
=pod
while (my $line = <INPUT>) {
	if ($line =~ /^model/) {
		chomp($line);
		my @data = split (" ", $line);
		my @geometry = split ("_", $data[0]);
		my $omega = $geometry[1];
		my $theta = $geometry[2];
		my $hBond = $data[-1];
	
		#Add blank line if x value changes--required to change scan in pm3d plot
		if ($omega != $omega[-1]) {
			push (@omega, "");
			push (@theta, "");
			push (@hBond, "");
		}
	
		#my $omegaPrime = (10 * $data[1] / 9) + (200 * $data[2] / 27);
		my $omegaPrime = $omega;
		push (@omega, $omegaPrime);
		push (@theta, $theta);
		push (@hBond, $hBond);
	}
}
=cut

print scalar @theta;
print " " . scalar @omegaB;
print " " . scalar @hBond;
#Chart object
my $chart = Chart::Gnuplot->new (
	output	=> $output,
	title	=> "HBond energy",
	view	=> 'map',
	xlabel	=> 'omegaBPrime',
	ylabel	=> 'theta',
	zlabel	=> 'Hydrogen Bonding',
	cbrange => '[-18:0]',
	yrange  => [-55,55],
	size	=> 'square',
	palette => 'negative',
	palette	=> 'defined (-16 "red", \
		-12 "blue",\
		-8 "green",\
		-4 "yellow",\
		0 "white")' 
);

#Data set object
my $dataSet = Chart::Gnuplot::DataSet->new (
	xdata	=> \@omegaB,
	ydata	=> \@theta,
	zdata	=> \@hBond,
	style	=> 'pm3d',
);

#plot chart
$chart->plot3d($dataSet);

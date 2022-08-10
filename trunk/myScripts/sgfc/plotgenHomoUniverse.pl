#!usr/bin/perl
use warnings;
use strict;
use Chart::Gnuplot;

# Parse text file
open INPUT, shift or die "ERROR: Could not open INPUT!\n";
my $output = shift;
my $Z = shift;
my (@theta, @omega, @hBond);
while (my $line = <INPUT>) {
	if ($line =~ /^[+-]?\d+/) {
		chomp($line);
		my @data = split ("\t", $line);
		if ($data[2] == $Z) {
			my $omega = $data[0];
			my $theta = $data[1];
			my $hBond = $data[8];
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
print " " . scalar @omega;
print " " . scalar @hBond;
print "\n";
#Chart object
my $chart = Chart::Gnuplot->new (
	output	=> $output,
	title	=> "HBond energy",
	view	=> 'map',
	xlabel	=> 'omegaPrime',
	ylabel	=> 'theta',
	zlabel	=> 'Hydrogen Bonding',
	size	=> 'square',
	palette	=> 'defined (-16 "red", -12 "blue",-8 "green",-4 "yellow", 0 "white")',
	xrange  => '[0:100]',
	yrange  => '[-55:55]',
	cbrange => '[-16:0]',
	#palette => 
);

#Data set object
my $dataSet = Chart::Gnuplot::DataSet->new (
	xdata	=> \@omega,
	ydata	=> \@theta,
	zdata	=> \@hBond,
	style	=> 'pm3d',
);

#plot chart
$chart->plot3d($dataSet);

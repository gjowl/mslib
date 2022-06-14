#!usr/bin/perl
use strict;
use Chart::Gnuplot;


my $datafile = shift;
my $outputfile = shift;


open DATAFILE, $datafile or die "ERROR: COULD NOT OPEN INPUT!\n";

my $header = <DATAFILE>;
my @colNames = split("\t", $header);



my $chart = Chart::Gnuplot->new (
	file	=> $datafile,
	output	=> $outputfile,
	title	=> "Boxplots",
	ylabel	=> "Cell Length"
);

$chart->set(style	=> 'boxplot');
$chart->set(factor	=> 1.0);
$chart->set(xtics	=> "", 1);
$chart->set(yrange	=> '[*:*]');

for (my $i = 1; $i < @colNames; $i++) {
	my $command = "xtics add(word($colNames[$i],$i) $i)";
	$chart->command($command);
	$chart->using("$i:$i");
	$chart->plot();
}


#!usr/bin/perl
#
#This program will read in an upper-triangle pairwise distance matrix
#of RMSD values and cluster them based on a given RMSD threshold.
#Starting with the energy of the best model, it will find all models
#that are within the given threshold. This process will repeat with the
#next best energy, and so on until no more PDBs are left. At the end it 
#will print out a list of the clusters along with useful statistics, like
#a histogram of clusters with number of models, best energy per cluster,
#average energy for each cluster, etc.


use File::Basename;
use File::Find;
use Getopt::Long;
use strict;


#TODO Read in CSV and put into matrix
#TODO Check to make sure matrix is upper triangle
#TODO Read in file of models and energies
#
#########################################################################
#	Options
########################################################################
my $pairwiseMatrix = '';	#CSV File containing the pairwise RMSDs (must be upper traingle or symmetric)
my $headersFile = '';		#File containing the names of each model, along with energy if appropriate (NAME,ENERGY), same order as matrix
my $outputFile = '';		#File to store the output
my $threshold = '';		#RMSD Threshold for clustering
GetOptions (
	'pairwiseMatrix=s' => \$pairwiseMatrix,
	'headersFile=s' => \$headersFile,
	'outputFile=s' => \$outputFile,
	'threshold=f'	=> \$threshold
);


my @matrix;
open MATRIX, "<$pairwiseMatrix" or die "ERROR: DID NOT OPEN MATRIX FILE\n";
while (my $line = <MATRIX>) {
	chomp $line;
	if ($line =~ /==RMSD==/../^=+$/) {
		next if ($line =~ /==RMSD==/ || $line =~ /^=+$/);
		my (@matrixLine) = split (/\s/, $line);
		push (@matrix, \@matrixLine);
	}
}
close MATRIX;
my @models;
my @energies;
my @modelInfo;
open HEADERS, "<$headersFile" or die "ERROR: DID NOT OPEN HEADERS FILE\n";

print "Matrix read! Converting to symmetric:\n";
#Convert matrix to fully symmetric matrix: cell ji = cell ij
for (my $i = 0; $i < @matrix; $i++) {
	for (my $j = $i+1; $j < @{$matrix[$i]}; $j++) {
		my $rmsd = $matrix[$i][$j];
		$matrix[$j][$i] = $rmsd;
	}
}


print "Reading header file...\n";
my $index = 0;
while (my $line = <HEADERS>) {
	chomp $line;
	my %info;
	my ($model, $energy) = split(",", $line);
	push (@models, $model);
	push (@energies, $energy);
	$info{'index'} = $index;
	$info{'model'} = $model;
	$info{'energy'} = $energy;
	push (@modelInfo, \%info);
	$index ++;
}
close HEADERS;

#Use clustering algorithm to cluster models
#Sort models by energy
#Pick best model and find all columns with rmsd <= threshold
#Push model info of clustered models to separate array
#Iterate through sorted energies and check if they are in a cluster
#Iterate through model rmsds and check if they are in a cluster
my @energySorted = sort { $a->{'energy'} <=> $b->{'energy'} } @modelInfo;

my %clusteredModels;
my @clusters;
print "Clustering models...\n";
for (my $i = 0; $i < @energySorted; $i++) {
	my %centroid = %{$energySorted[$i]};
	my $index = $centroid{'index'};
	#Check if the centroid has already been clustered 
	if (!exists $clusteredModels{$index}) {
		$clusteredModels{$index} = 'clustered!';
		my @cluster;
		push (@cluster, \%centroid);
		for (my $j = 0; $j < @{$matrix[$index]}; $j++) {
			if(!exists $clusteredModels{$j} && $matrix[$index][$j] <= $threshold) {
				$modelInfo[$j]->{'rmsd'} = $matrix[$index][$j];
				push (@cluster, $modelInfo[$j]);
				$clusteredModels{$j} = 'clustered!';
			}
		}
		push (@clusters, \@cluster);
	}
}


#Print out clustered models
for (my $i = 0; $i < @clusters; $i++) {
	print "=====Cluster $i==========\n";
	my @sortedCluster = sort{ $a->{'energy'} <=> $b->{'energy'} } @{$clusters[$i]};
	for (my $j = 0; $j < @sortedCluster; $j++) {
		print "$sortedCluster[$j]->{'model'}\n";
	}
	print "\n";
}

#Print out cluster statistics
print "Number of clusters: " . scalar(@clusters) . "\n";
for (my $i = 0; $i < @clusters; $i++) {
	print "=====Cluster $i==========\n";
	print "Size:	" . scalar (@{$clusters[$i]}) . "\n";
	print "Best Model:	${$clusters[$i]}[0]->{'model'}\n";
	print "Best Energy:	${$clusters[$i]}[0]->{'energy'}\n";

	my $clusterRMSD = 0;
	my $clusterEnergy = 0;
	for (my $j = 0; $j < @{$clusters[$i]}; $j++) {
		$clusterRMSD += ${$clusters[$i]}[$j]->{'rmsd'};
		$clusterEnergy += ${$clusters[$i]}[$j]->{'energy'};
	}
	$clusterRMSD /= scalar(@{$clusters[$i]});
	$clusterEnergy /= scalar(@{$clusters[$i]});
	print "Average RMSD to centroid:	$clusterRMSD\n";
	print "Average cluster energy: $clusterEnergy\n";
}
for (my $i = 0; $i < @clusters; $i++) {
	print "$i " . scalar(@{$clusters[$i]}) . "\n";
}

				


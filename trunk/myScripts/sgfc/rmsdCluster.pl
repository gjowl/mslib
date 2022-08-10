#!usr/bin/perl
use strict;

# This program takes a pairwise distance matrix (in this case we will assume
# that it's an RMSD matrix) and group them into clusters.
#
my $clusterThreshold = 0.5;

my @proteins;		# 1-D array containing the names of the proteins used in the pairwise alignment
my @proteinMask;	# Mask to indicate if the protein has been assigned to a cluster; 1 is assigned and 0 is unassigned

my @rmsdMatrix;
my @clusters;

my $file = shift;
open FILE, $file or die "Did not open csv file!";
my @fileArray = <FILE>;

my $proteinNames = shift (@fileArray);
chomp $proteinNames;
@proteins = split(",", $proteinNames);
shift @proteins;

foreach my $line (@fileArray) {
	chomp $line;
	my (@pairwiseRMSD) = split (",", $line);
	shift @pairwiseRMSD;
	push (@rmsdMatrix, \@pairwiseRMSD);
}
close FILE;
while (@proteins) {
	my $int = int(rand(@proteins));

	my @cluster;
	print "New cluster beginning with $proteins[$int]\n";
	push (@cluster, $proteins[$int]);
	my $bait = splice (@proteins, $int, 1);


	my @alignment = @{$rmsdMatrix[$int]};
	#Remove the matching pairwise alignment score, which will be 0
	splice (@alignment, $int, 1);
	#Remove the column and row corresponding to the protein that's forming the new cluster
	splice (@rmsdMatrix, $int, 1);
	for (my $i = 0; $i < @rmsdMatrix; $i++) {
		splice (@{$rmsdMatrix[$i]}, $int-1, 1);
	}

	my @splice;
	for (my $i = 0; $i < @alignment; $i++) {
		if ($alignment[$i] < $clusterThreshold) {
			print "$bait to $proteins[$i]: $alignment[$i] is less than $clusterThreshold\n";
			push (@cluster, $proteins[$i]);
			push (@splice, $i);	#We will remove these clustered entries at the end of the loop
		}
	}

	#Remove the clustered columns and rows as well as the name of the protein from consideration
	#thisissohacky
	my $counter;
	foreach my $i (@splice) {
		splice (@proteins, $i-$counter, 1);
		splice (@rmsdMatrix, $i-$counter, 1);
		for (my $j = 0; $j < @rmsdMatrix; $j++) {
			splice (@{$rmsdMatrix[$j]}, $i-$counter, 1);
		}
		$counter++;
	}

	push (@clusters, \@cluster);
}

for (my $i = 0; $i < @clusters; $i++) {
	print "Cluster " . ($i+1) . ": \n";
	for (my $j = 0; $j < @{$clusters[$i]}; $j++) {
		print "$clusters[$i][$j]\n";
	}
	print "\n";
}

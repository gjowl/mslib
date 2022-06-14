#!usr/bin/perl

#This program will recursively search for models to concatenate
#into a large multi-model PDB file
#
#
#
#




use File::Basename;
use File::Find;
use Getopt::Long;
use Cwd;
use strict;

#TODO Find all models in a directory recursively--maybe sort?
#TODO Cat all models into a single file
#	TODO Make option to just do CAs
#TODO send that file to pairwiseRMSD
#TODO Take output from pairwiseRMSD and perform clustering analysis
#TODO Write multi-model PDBs for each cluster
#

my $baseDir = '';
my $modelName = '';		#Assume model names are all the same
my $outputDir = '';		#Where to output multimodel PDB files and results
my $allAtoms = '';		#Only pull out CA atoms unless otherwise stated
my $rmsdThreshold = '';
my $pairwiseRMSD = "/exports/home/scondon/mslib/trunk_AS/bin/pairwiseRMSD";

GetOptions (
	'baseDir=s'		=> \$baseDir,
	'modelName=s'		=> \$modelName,
	'outputDir=s'		=> \$outputDir,
	'allAtoms'		=> \$allAtoms,
	'rmsdThreshold=f'	=> \$rmsdThreshold
);

if ($baseDir eq '') {
	$baseDir = getcwd();
}
if ($outputDir eq '') {
	$outputDir = getcwd();
}
if ($modelName eq '') {
	print "ERROR: Model name not specified!\n";
	die;
}


opendir ( BASEDIR, $baseDir ) or die "ERROR: CANNOT OPEN DIRECTORY $baseDir!\n";

my @subDirs = sort { $a <=> $b } readdir(BASEDIR);
my @modelDirs;

my ($modelFile, $modelDir, $modelSuffix) = fileparse($modelName, qr/\.[^.]*/);
print "$modelName $modelFile $modelDir $modelSuffix\n";
my $multiPDB = "$outputDir/$modelFile-all.pdb";
my @modelNames;
my @modelEnergies;
open MULTIPDB, ">$multiPDB" or die "ERROR: DID NOT OPEN MULTI PDB FILE $multiPDB!\n";

foreach my $sub (@subDirs) {
	my $model = "$baseDir/$sub/$modelName";
	if (-f $model) {
		print "$model\n";
		push (@modelDirs, $sub);
		my $modelEnergy = '';
		open MODEL, $model or die "ERROR: DID NOT OPEN MODEL FILE $model!\n";

		print MULTIPDB "REMARK PDB FILE $model\n";
		print MULTIPDB "MODEL\n";
		while (my $line = <MODEL> ) {
			chomp $line;
			if ($allAtoms) {
				print MULTIPDB "$line\n";
			} elsif (!$allAtoms && $line =~ /^ATOM.*CA/) {
				print MULTIPDB "$line\n";
			}
	
			if ($line =~ /Total/ ) { #Contains the energy of the model we are interested in
				my ($remark, $total, $energy, $interactions) = split (/\s+/, $line);
				$modelEnergy = $energy;
			}
				
		}
		print MULTIPDB "ENDMDL\n";
		push (@modelNames, $model);
		push (@modelEnergies, $modelEnergy);
	}

}
close MULTIPDB;

for (my $i = 0; $i < @modelNames; $i++) {
	print "$modelDirs[$i]	$modelNames[$i]	$modelEnergies[$i]\n";
}
#Now that we have a file containing all our models, perform pairwise alignments on all of them

my $pairwiseCmd = "$pairwiseRMSD --pdb $multiPDB --outputDir $outputDir";
system($pairwiseCmd);

my $matrixFile = "$outputDir/pairwiseRMSD.tsv";
my @matrix;


open MATRIX, "<$matrixFile" or die "ERROR: DID NOT OPEN MATRIX FILE\n";
while (my $line = <MATRIX>) {
	chomp $line;
	print "$line\n";
	my (@matrixLine) = split (/\s/, $line);
	push (@matrix, \@matrixLine);
}
close MATRIX;

print "Matrix read! Converting to symmetric:\n";
#Convert matrix to fully symmetric matrix: cell ji = cell ij
for (my $i = 0; $i < @matrix; $i++) {
	for (my $j = $i+1; $j < @{$matrix[$i]}; $j++) {
		my $rmsd = $matrix[$i][$j];
		$matrix[$j][$i] = $rmsd;
	}
}

my @modelInfo;
for (my $i = 0; $i < @modelNames; $i++) {
	my %info;
	$info{'index'} = $i;
	$info{'subDir'} = $modelDirs[$i];
	$info{'model'} = $modelNames[$i];
	$info{'energy'} = $modelEnergies[$i];
	push (@modelInfo, \%info);


}


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
			if(!exists $clusteredModels{$j} && $matrix[$index][$j] <= $rmsdThreshold) {
				$modelInfo[$j]->{'rmsd'} = $matrix[$index][$j];
				push (@cluster, $modelInfo[$j]);
				$clusteredModels{$j} = 'clustered!';
			}
		}
		push (@clusters, \@cluster);
	}
}


#Print out clustered models
open LOGFILE, ">$outputDir/clusterModels.log" or die "ERROR: DID NOT OPEN LOGFILE!\n";
for (my $i = 0; $i < @clusters; $i++) {
	print LOGFILE "=====Cluster $i==========\n";
	my @sortedCluster = sort{ $a->{'energy'} <=> $b->{'energy'} } @{$clusters[$i]};
	for (my $j = 0; $j < @sortedCluster; $j++) {
		print LOGFILE "$sortedCluster[$j]->{'model'}\n";
	}
	print LOGFILE "\n";
}

#Print out cluster statistics
print LOGFILE "Clustering $modelName models from $baseDir at a RMSD threshold of $rmsdThreshold:\n";
print LOGFILE "Number of clusters: " . scalar(@clusters) . "\n";
for (my $i = 0; $i < @clusters; $i++) {
	print LOGFILE "=====Cluster $i==========\n";
	print LOGFILE "Size:	" . scalar (@{$clusters[$i]}) . "\n";
	print LOGFILE "Best Model:	${$clusters[$i]}[0]->{'model'}\n";
	print LOGFILE "Best Energy:	${$clusters[$i]}[0]->{'energy'}\n";

	my $clusterRMSD = 0;
	my $clusterEnergy = 0;
	for (my $j = 0; $j < @{$clusters[$i]}; $j++) {
		$clusterRMSD += ${$clusters[$i]}[$j]->{'rmsd'};
		$clusterEnergy += ${$clusters[$i]}[$j]->{'energy'};
	}
	$clusterRMSD /= scalar(@{$clusters[$i]});
	$clusterEnergy /= scalar(@{$clusters[$i]});
	print LOGFILE "Average RMSD to centroid:	$clusterRMSD\n";
	print LOGFILE "Average cluster energy: $clusterEnergy\n";
}
for (my $i = 0; $i < @clusters; $i++) {
	print "$i " . scalar(@{$clusters[$i]}) . "\n";
	my $clusterFile = "$outputDir/$modelFile-cluster$i.pdb";
	open CLUSTERPDB, ">$clusterFile" or die "ERROR: DID NOT OPEN CLUSTER FILE $clusterFile!\n";
	my @sortedCluster = sort{ $a->{'energy'} <=> $b->{'energy'} } @{$clusters[$i]};
	for (my $j = 0; $j < @sortedCluster; $j++) {
		
		my $modelFile = $sortedCluster[$j]->{'model'};
		open CLUSTERMODEL, $modelFile or die "ERROR: DID NOT OPEN $modelFile\n";
		print CLUSTERPDB "MODEL\n";		
		print CLUSTERPDB "REMARK $modelFile\n";

		while (my $line = <CLUSTERMODEL>) {
			chomp $line;
			print CLUSTERPDB "$line\n";
		}
		print CLUSTERPDB "ENDMDL\n";
		close CLUSTERMODEL;
	}
	print CLUSTERPDB "END\n";
	close CLUSTERPDB;

}

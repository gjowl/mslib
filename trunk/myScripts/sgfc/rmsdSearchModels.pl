#!/usr/bin/perl
# This program identifies similar models found in CATM runs compared to a reference.
# For example, in a mutagenesis experiment it can find structure models that are closest
# to the wild-type structure.  This way you can examine the energy differences between
# similar structures.
#
# For each protein predicted by CATM in the directory, we will align each of the models in that
# run to the reference PDB using alignMolecules.  The model with the lowest RMSD will be
# selected and its name will be printed.
use strict;
use Getopt::Long;	# Option parser for Perl


#The main directory holds all of the individual CATM runs.  We will loop over each 
#subdirectory to get the names of the files we need and the information about
#the interface that we will align over.
#
##############################################################################
#	Options
##############################################################################

my $allCa = '';		# Use all alpha carbons in the alignment (useful for SDM experiments)
my $mainDir = '';	# Directory where the CATM runs are located
my $refPdb = '';	# PDB That we will align models to
my $threshold = 0;	# Main 


GetOptions (
	'allCa' => \$allCa,
	'mainDir=s' => \$mainDir,
	'refPdb=s' => \$refPdb,
);


opendir(DIR, $mainDir);
my @subDirs;
while (my $subDir = readdir (DIR)){
	if (-d "$mainDir/$subDir") {
	push (@subDirs, $subDir);
	}
}
closedir DIR;



#for (my $i = 0; $i < @subDirs; $i++) {
#	print $subDirs[$i] . "\n";
#}

my @outputs;
foreach my $subDir (@subDirs) {
	my $path = "$mainDir/$subDir/";


	opendir SUBDIR, $path or next;
	while (my $file = readdir(SUBDIR)) {
		if ($file =~ /\d\d\.txt$/) {
			my %outputHash;
			$file =~ /(\d\d)(\.txt$)/;	#Get the digits coming right before the extension to get the model number
			my $model = $1;
			#print "$file is model number $model\n";
			$outputHash{"model"} = $model;
			open FILE, $file;
			while (my $line = <FILE>) {
				chomp $line;
				my ($key, $val) = split (" ", $line);
				$outputHash{$key} = $val;
			}
			#print "TEXT FILE BEFORE: $file\n";
			my $pdbFile = $file;
			$pdbFile =~ s/\.txt$/\.pdb/g;
			#print "TEXT FILE: $file\n";
			#print "PDB FILE: $pdbFile\n";
			$outputHash{"accession"} = $subDir;
			$outputHash{"pdbFile"} = "$path/$pdbFile";
			close FILE;

			# Calculate range of positions to align between structures based on the interface
			# Notation: 00100110122012200100010000000
			# 0 is not interface, 1 is interfacial, and 2 is a crossing point residue
			# We want to align the crossing point residues
			my $interface = $outputHash{"interfacial"};
			my ($istart, $alignStart, $alignEnd);
			if ($interface =~ /22..22/g) {
				$istart = pos ($interface) - 6 + $outputHash{"resStart"};
				$alignStart = $istart;
				$alignEnd = $alignStart + 6;
		
				$outputHash{"istart"} = $istart;
				$outputHash{"alignStart"} = $alignStart;
				$outputHash{"alignEnd"} = $alignEnd;
			}
			push (@outputs, \%outputHash);
		}
			close SUBDIR;
	}
			

}

for (my $i = 0; $i < 100; $i++) {
	my $command;
	if ($allCa) {
		$command = "alignMolecules --pdb1 $refPdb --pdb2 $outputs[$i]->{pdbFile} --sele1 Chain A and name ca --sele2 Chain A and name ca --sele1 Chain B and name ca --sele2 Chain B and name ca --noOutputPdb\n";
	} 
	my @alignment = qx($command);
	foreach my $line (@alignment) {
		if ($line =~ /^RMSD/g) {
			my($junk, $rmsd) = split (" ", $line);
			print "$outputs[$i]->{'accession'} model $outputs[$i]->{'model'} to reference: $rmsd\n";

			}
	}
	
}

=pod
my @rmsdArray;
for (my $i = 0; $i < @outputs; $i++) {
	for (my $j = $i+1; $j < @outputs; $j++) {

		# Print the command that will align the two pdb files
		my $pdbKey = "pdbFile";
		my $alignStartKey = "alignStart";
		my $alignEndKey = "alignEnd";
		my $command;
		if ($useAll) {
			$command = "alignMolecules --pdb1 $outputs[$i]->{$pdbKey} --pdb2 $outputs[$j]->{$pdbKey} --sele1 Chain A and name ca --sele2 Chain A and name ca --sele1 Chain B and name ca --sele2 Chain B and name ca --noOutputPdb\n";
		} else {
		$command = "alignMolecules --pdb1 $outputs[$i]->{$pdbKey} --pdb2 $outputs[$j]->{$pdbKey} --sele1 Chain A and name ca and resi $outputs[$i]->{$alignStartKey}-$outputs[$i]->{$alignEndKey} --sele2 Chain A and name ca and resi $outputs[$j]->{$alignStartKey}-$outputs[$j]->{$alignEndKey} --sele1 Chain B and name ca and resi $outputs[$i]->{$alignStartKey}-$outputs[$i]->{$alignEndKey} --sele2 Chain B and name ca and resi $outputs[$j]->{$alignStartKey}-$outputs[$j]->{$alignEndKey} --noOutputPdb\n";
		}
		print $command;
		chomp ($command);
		my @alignment = qx($command);
		foreach my $line (@alignment) {
			if ($line =~ /^RMSD/g) {
				my($junk, $rmsd) = split (" ", $line);
				# print "$outputs[$i]->{accession} to $outputs[$j]->{accession}: $rmsd\n";
				$rmsdArray[$i][$j] = $rmsd;
				}
		}

	}
}

open (MATRIXFILE, ">./clusterCATM_Matrix.csv");
print MATRIXFILE "Protein,";
foreach my $prot (@outputs) {
	print MATRIXFILE "$prot->{accession},";
}
print MATRIXFILE "\n";
for (my $i = 0; $i < @outputs; $i++) {
	print MATRIXFILE "$outputs[$i]->{accession},";
	for (my $j = 0; $j < @outputs; $j++) {
		print MATRIXFILE "$rmsdArray[$i][$j],";
	}
	print MATRIXFILE "\n";
}
close MATRIXFILE;


for (my $i=0; $i<@dir) {
	if ($dir[$i] =~ /^[A-Z0-9]{6}$/) {
		open FILE "$dir[$i]/$dir[$i]_01.txt" || die ("Error\n";)
	if $dir[$i] has the correct format
		open file

open(OUTPUT, $outputFile);
while (my $line = <OUTPUT>) {
	if ($line =~ m/interfacial /) {
		my $interface =~ /{d}*/;
	}
my $interface = "00000010011002201220111011000";
my $resStart = 23;
if ($interface =~ /22..22/g) {
	my $istart = pos ($interface) - 6 + $resStart;
	my $alignStart = $istart -4;
	my $alignEnd = $istart + 10;
#	print "alignMolecules --pdb1 $pdb1 --pdb2 $pdb2 --sele1 \"Chain A & name ca & resi $alignStart-$alignEnd\" --sele2 \"Chain B & name ca & resi $alignStart-$alignEnd\";
}
=cut

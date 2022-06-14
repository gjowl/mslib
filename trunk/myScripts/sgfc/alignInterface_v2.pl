#!/usr/bin/perl
# This program clusters a set of CATM dimers by RMSD in a pairwise fashion.
# It finds the dimer interface provided in the CATM output file and uses the
# precompiled program AlignMolecules to align the two dimers around the 
# interfaces.  
#
# interfacial 00000010011002201220111011000

use strict;
use List::Util qw[min max];
my $R0 = 3.0;
#The main directory holds all of the individual CATM runs.  We will loop over each 
#subdirectory to get the names of the files we need and the information about
#the interface that we will align over.
my $mainDir = shift;
my $useAll;
my $outputFile = shift;
if (!$outputFile) {
	$outputFile = ">./CATM_RMSD.csv";
}
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
	open FILE, "$path/$subDir"."_01.txt" or next;
	my %outputHash;
	while (my $line = <FILE>) {
		chomp $line;
		my ($key, $val) = split (" ", $line);
		$outputHash{$key} = $val;
	}

	close FILE;	
	$outputHash{"accession"} = $subDir;
	$outputHash{"pdbFile"} = "$path/$subDir"."_01.pdb";

	
#	my $interface = $outputHash{"interfacial"};
#	my ($istart, $alignStart, $alignEnd);
#	if ($interface =~ /22..22/g) {
#		$istart = pos ($interface) - 6 + $outputHash{"resStart"};
#		$alignStart = $istart;
#		$alignEnd = $alignStart + 6;
#
#		$outputHash{"istart"} = $istart;
#		$outputHash{"alignStart"} = $alignStart;
#		$outputHash{"alignEnd"} = $alignEnd;
#	}
	push (@outputs, \%outputHash);
}

my @rmsdArray;
for (my $i = 0; $i < @outputs; $i++) {
	for (my $j = $i+1; $j < @outputs; $j++) {

		# Print the command that will align the two pdb files
		my $pdbKey = "pdbFile";
		my $alignStartKey = "alignStart";
		my $alignEndKey = "alignEnd";
		my $command;

		my ($alignStart, $alignEnd) = &max_Align($outputs[$i]->{"interfacial"}, $outputs[$j]->{"interfacial"});

		my $resiStart_i = $alignStart + $outputs[$i]->{"resStart"};
		my $resiEnd_i   = $alignEnd + $outputs[$i]->{"resEnd"};
		my $resi_i = "$resiStart_i-$resiEnd_i";

		my $resiStart_j = $alignStart + $outputs[$j]->{"resStart"};
		my $resiEnd_j   = $alignEnd + $outputs[$j]->{"resEnd"};
		my $resi_j = "$resiStart_j-$resiEnd_j";

		if ($useAll) {
			$command = "alignMolecules --pdb1 $outputs[$i]->{$pdbKey} --pdb2 $outputs[$j]->{$pdbKey} --sele1 Chain A and name ca --sele2 Chain A and name ca --sele1 Chain B and name ca --sele2 Chain B and name ca --noOutputPdb\n";
		} else {
		$command = "alignMolecules --pdb1 $outputs[$i]->{$pdbKey} --pdb2 $outputs[$j]->{$pdbKey} --sele1 Chain A and name ca and resi $resi_i --sele2 Chain A and name ca and resi $resi_j --sele1 Chain B and name ca and resi $resi_i --sele2 Chain B and name ca and resi $resi_j --noOutputPdb\n";
		}
		print $command;
		chomp ($command);
		my @alignment = qx($command);
		foreach my $line (@alignment) {
			if ($line =~ /^RMSD/g) {
				my($junk, $rmsd) = split (" ", $line);
				my $qScore = (($alignEnd-$alignStart)^2)/((1+(($rmsd/$R0)^2))*length($outputs[$i]->{"interfacial"})*length($outputs[$j]->{"interfacial"}));
				print "$outputs[$i]->{accession} to $outputs[$j]->{accession}: $rmsd \n";
				$rmsdArray[$i][$j] = $rmsd;
				}
		}

	}
}


# Determines the maximum number of residues we can align from the central interface of the two structures
# For this particular problem the structures are idealized alpha helix dimers, so we can simply look at
# the same number of residues starting from the center.
	# Notation: 00100110122012200100010000000
	# 0 is not interface, 1 is interfacial, and 2 is a crossing point residue

sub max_Align () {
	my ($interface1, $interface2) = @_;
	my ($intLeft1, $intLeft2, $intRight1, $intRight2);
	if ($interface1 =~ /22..22/) {
		$intLeft1 = pos ($interface1);
		$intRight1 = length($interface1) - $intLeft1;
	}
	if ($interface2 =~ /22.22/) {
		$intLeft2 = pos ($interface2);
		$intRight2 = length($interface2) - $intLeft2;
	}
	my $alignStart = min ($intLeft1, $intLeft2);
	my $alignEnd = min ($intRight1, $intRight2);
	return ($alignStart, $alignEnd);
}







open (MATRIXFILE, '>', $outputFile) or die "Could not make $outputFile!\n";
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
=pod
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

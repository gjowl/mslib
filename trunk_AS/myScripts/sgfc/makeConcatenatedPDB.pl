# This script concatenates a bunch of pdb files in ordered job directories for clustering
#!/usr/bin/perl
use strict;
use Cwd;
#
my $baseDir = getcwd;
my $fileName = shift;
my $outputDir = shift;
my $exec = "/scratch/senesgrp/sgcondon/mslib/trunk_AS/bin/trajectoryRMSD";
my $sele = "name CA";

my $grep = 'grep  -e "[0123456789] \+N \+" -e "[0123456789] \+C \+" -e "[0123456789] \+CA \+" -e "[0123456789] \+O \+"';

open OUTPUT, ">$outputDir/all_backbone.pdb" or die "ERROR: DID NOT OPEN OUTPUT FILE!\n";
for (my $i = 1; $i <= 999; $i++) {
	my $folderIndex = sprintf "%04s", $i;
	my $subdir = $baseDir . "/" . $folderIndex;
	my $model = "$subdir/$fileName";
	if (-f $model) {
		print OUTPUT "MODEL\n";
		my $command = $grep . " $model";
		my $backbone .= qx/$command/;
		print OUTPUT "$backbone\n";
		print OUTPUT "ENDMDL\n";
	}
}


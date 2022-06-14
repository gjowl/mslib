#!usr/bin/perl
#This script will concatenate a set of files from a base directory
#
#usage: perl catFiles.pl /exports/home/scondon/FtsB-FtsL_MCR/ MC_Final.pdb allModels.pdb true;
use File::Slurp;


my $baseDir = shift;
my $fileName = shift;
my $outFile = shift;
my $isPDB = shift;


opendir ( DATA_DIR, $baseDir ) or die "ERROR: CANNOT OPEN $baseDir!\n";

my @subDirs = sort { $a <=> $b } readdir(DATA_DIR);

open (OUTFILE, ">$outFile") or die "ERROR: CANNOT OPEN $outFile!\n";
while (my $subDir = shift @subDirs) {
	my $file = "$baseDir/$subDir/$fileName";
	if (-f $file) {
		if ($isPDB) {
			print OUTFILE "MODEL\n";
		}
		my $fileText = read_file($file);
		print OUTFILE "$fileText\n";
		if ($isPDB) {
			print OUTFILE "ENDMDL\n";
		}
		print "$baseDir/$subDir/$fileName\n";
	}
}
close OUTFILE;

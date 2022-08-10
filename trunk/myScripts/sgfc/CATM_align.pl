#This program takes the .txt files generated from CATM and combines them into a single table for further analysis.

#!usr/bin/perl

use File::Basename;
use File::Find;
use Getopt::Long;
use strict;



########################################################################
#	Options
########################################################################
my $inputDirectory = '';	#Directory where CATM location is stored
my $energyterms = '';		#Boolean to search for the individual energy terms
my $outputFile = '';		#File to store the output
GetOptions (
	'inputDirectory=s' => \$inputDirectory,
	'energyterms' => \$energyterms,
	'outputFile=s' => \$outputFile
);


my @subDirs;
my @outputs;

-d $inputDirectory or die "Input directory $inputDirectory is invalid: $!\n";

opendir(DIR, $inputDirectory);
while (my $subDir = readdir (DIR)){
	if (-d "$inputDirectory/$subDir") {
		push (@subDirs, $subDir);
	}
}
closedir DIR;

foreach my $subDir (@subDirs) {
	my $path = "$inputDirectory/$subDir/";
	open FILE, "$path/$subDir"."_01.txt" or next;
	my %outputHash;
	while (my $line = <FILE>) {
		chomp $line;
		my ($key, $val) = split (" ", $line);
		$outputHash{$key} = $val;
	}

	close FILE;
	open ENERGY, "$path/$subDir"."_01.energy" or next;
	while (my $line = <ENERGY>) {
		chomp $line;
		my ($key, $val) = split (" ", $line);
		if ($key != "MODEL") {
			$outputHash{$key} = $val;
		}
	}
	$outputHash{"accession"} = $subDir;
	$outputHash{"pdbFile"} = "$path/$subDir"."_01.pdb";

	
	push (@outputs, \%outputHash);
}
=pod
find(\&read_protein_files, $inputDirectory);


sub read_protein_files
{
	my $file = $_;
	if ($file =~ /^	# looking for a file name like "A2RRL7_100_04_01.txt"
	.{6}		# match exactly six characters
	_		# match "_" and stuff in the middle
	01		# match the first model
	\.txt		# match .txt /x) {		# /x to allow whitespace for commenting

			
		my $i = 0;
		if (-r $file) {	
			open(TXT, "$file");
			print "$file\n";
			LINE: while (my $line = <TXT>) {
				if ($line =~ /hbonds/) {	#ignore the multiple models and h bonds in each text file and add the name and data tags to the row


					# first, open up the energy file and add these last two terms to the row
					my $fileName = basename($file, ".txt");
					my $energyFile = "$fileName" . ".energy";
					open(ENERGY, "$energyFile") or die "Can't find energy file: $!\n";

					while (my $energyLine = <ENERGY>) {
						(my $energyType, my $energyValue) = split (" ", $energyLine);
						push (@{ $proteinArray[$i] }, $energyValue);
						$i+=1;
					}
					close ENERGY;
					# add the filename's metadata to the row
					(my @nameArray) = split ("_", $fileName);
					foreach (@nameArray) {
						push (@{ $proteinArray[$i] }, $_);
						$i+=1;
					}
					$i = 0;
					last LINE;
				}
				(my $name, my $value) = split (" ", $line);
				push (@{ $proteinArray[$i] }, $value);
				$i+=1;
			}
			close TXT;

		} else {
			print "Can't open file $file: $!\n";
		}

	}
}
=cut


# read text file and poulate the next row in the array with the relevant info
# The two print statements below give you the dimensions of the table so you get an idea of whether it worked.
#print "The scalar atproteinArray is length " . scalar(@proteinArray) . "\n";
#print "The subarray is length " . scalar(@{$proteinArray[0]}) . "\n";

# make a new file to save your data to.  This file can be copied and pasted into Excel using space delimiters

open(OUTPUT, ">$outputFile") or die "Can't make file: $!\n";
foreach my $key (sort(keys $outputs[0])) {
	print OUTPUT $key . "	";
}
print OUTPUT "\n";

for (my $i = 0; $i < @outputs; $i++) {
	
	foreach my $key (sort(keys $outputs[$i])) {
		print OUTPUT "$outputs[$i]->{$key}	";
	}
	print OUTPUT "\n";
}


#This script will make a submit file for stockholmMSA.pl
#by reading the output of CATM_align and identifying the 
#position of glycines within the analyzed domain.
#

#!usr/bin/perl
use Getopt::Long;
use Bio::SeqIO;
use strict;



########################################################################
#	Options
########################################################################
my $input = '';			#File containing the CATM sequence information
my $inputDir = '';		#Directory of MSA files in Stockholm Format
my $submitFile = 'stockholmMSAsubmit.sh';	#Shell script for submission to stockholmMSA.pl
my $interface = -1;		#Option to filter residues by interface type (0, 1, 2, n1, c1, or c5)
my $output = '';		#Where to write the results of stockholmMSA to

my $commandline = "perl " . $0 . " " . (join " ", @ARGV);	#Stores the command line input for your records


GetOptions (
	'input=s' => \$input,
	'submitFile=s' => \$submitFile,
	'inputDir=s' => \$inputDir,
	'interface=s' => \$interface,
	'output=s' => \$output,
) or exit (1);



=pod
my $data = Bio::SeqIO->new(	-input => $input,
				-format => 'table'
				#-header => 1,
				#-delimiter => ','
			);
=cut
open OUTPUT, ">$submitFile" or die "ERROR: Could not open $submitFile!\n";
print OUTPUT "#$commandline\n";
close OUTPUT;

open CATMSEQ, $input or die "Could not open $input\n";
while (my $line = <CATMSEQ>) {
	chomp $line;
	my ($accession, $seq, $interfacial, $start, $end, $pfam) = split (",", $line);
	if ($pfam =~ /^P/) {
		my $msaFile = $inputDir . "/" . $pfam . ".stk";
		my @positions;
		my @interfaces;
		my $offset = 0;
		#Get the position of the C1 residue (3rd "2" on the interface list)
		#and check for glycines at C1 + 4 and C1 - 4
		
		my @crossing;
		my $crossingres;
		while (1) {
			$crossingres = index ($interfacial, "2", $crossingres + 1);
			last if ($crossingres < 0);
			my $crossPos = $crossingres + $start;
			push (@crossing, $crossPos);
		}
		my $n1 = $crossing[0];
		my $c1 = $crossing[2];
		my $c5 = $crossing[2] + 4;
		#print "$accession $crossing[0] $crossing[1] $crossing[2] $crossing[3]\n";

		while (1) {
			my $position = index($seq, "G", $offset);
			last if ($position < 0);
			my $resInterface = substr($interfacial, $position, 1);
			my $residue = $position + $start;
			if ($residue == $n1) {
				$resInterface = "n1";
			} elsif ($residue == $c1) {
				$resInterface = "c1";
			} elsif ($residue == $c5) {
				$resInterface = "c5";
			}
			
			
			if ($interface==-1) {	#If no interface type selected, default to all residues
				push (@positions, $residue);
				push (@interfaces, $resInterface);
			} elsif ( $resInterface eq $interface) {
				push (@positions, $residue);
				push (@interfaces, $resInterface)
			}
			$offset=$position + 1;
		}

		if (@positions) {
			open OUTPUT, ">>$submitFile" or die "ERROR: Could not open $submitFile!\n";
			print OUTPUT "perl /exports/home/scondon/mslib/trunk/myScripts/sgfc/stockholmMSA.pl --file $msaFile --protein $accession --output $output --residue ";
			foreach my $pos (@positions) {
				print OUTPUT "$pos ";
			}
			print OUTPUT "--interface ";
			foreach my $int (@interfaces) {
				print OUTPUT "$int ";
			}
			print OUTPUT "\n";
			close OUTPUT;
		}


	}
}

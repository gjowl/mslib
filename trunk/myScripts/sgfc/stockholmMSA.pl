#!usr/bin/perl
# This program will read a multiple sequence alignment in Stockholm (Pfam)
# format and search for particular residues of one protein.  For the 
# columns that those residues are in, the program will return the number
# and ratio of each amino acid in that column.
#
#
# Example usage:
# perl myScripts/sgfc/stockholmMSA.pl --file /exports/home/scondon//Desktop/pfam_alignments/PF01102_full.txt --protein P02724 --residue 98 102 113 
#
#
# Steps in the program:
# TODO Read stockholm MSA
# TODO Find reference protein in MSA
# TODO Find the specified residue in that protein
# TODO Count the amino acids in that column
# TODO Return the amino acid and gap counts

use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Getopt::Long;
use File::Basename;


########################################################################
#	Options
########################################################################
my $file = '';	#Multiple sequence alignment
my $refProtein = '';	#UNIPROT code that we will use as reference
my @residues;		#Residues in reference protein whose MSA columns we will examine	
my $keepgaps = '';	#Keep sequences with gaps in the columns, default false
my $output = '';	#Where to write the residue counts
#my $residue = '';		#Residue to count for conservation

my $commandline = "perl " . $0 . " " . (join " ", @ARGV);	#Stores the command line input for your records

GetOptions (
	'file=s' => \$file,
	'protein=s' => \$refProtein,
	'residue=i{1,}'=> \@residues,
	'output=s' => \$output

);

my $str = Bio::AlignIO->new (	-file => $file,
				-format => "stockholm" );
my $alignment = $str->next_aln();
my $msa = basename($file, ".stk");

print "\n======================= " . $msa . " =======================\n";
print $alignment->length . "\n";
print $alignment->num_residues . "\n";
print $alignment->num_sequences . "\n";
print $alignment->is_flush . "\n";
print $alignment->consensus_string(50) . "\n";


my @columns;
foreach my $seq ($alignment->each_seq) {
	if ($seq->accession_number =~ /$refProtein/) {
		
		print $seq->accession_number . "\n";
		# print $seq->start . "-" . $seq->end . " Length: " . $seq->length . "\n";
		# Hacky way to find my residue!
		# Iterate through sequence and ignore gaps, counting up until
		# you find the proper residue
		my $posCounter = $seq->start-1;	#Set poscounter to "0th" position of MSA
		for (my $i = 1; $i < $seq->length; $i++) {	#Bio::Seq is a 1-based object
			my $pos = $seq->subseq($i,$i);
			if ($pos =~ /[a-zA-Z]/) {
				$posCounter++;
				for (my $j = 0; $j < @residues; $j++) {
					if ($posCounter == $residues[$j]) {
						my %positions;
						$positions{'residue'} = $pos;
						$positions{'position'} = $residues[$j];
						$positions{'column'} = $i;
						print $pos . $posCounter . " is column $i in the MSA\n";
						push @columns, \%positions;
					}
	
				}
			}
			
		}

	}
}

my @conservation;

=pod

foreach my $i (@columns) {
	my $aln = $alignment->slice($i,$i,$keepgaps);
	my @consensus = $aln->consensus_conservation();
	my $num_consensus = $consensus[0] * $aln->num_sequences / 100;
	print "$consensus[0] percent, $num_consensus out of " . $aln->num_sequences . "\n";
	push @conservation, $aln->consensus_conservation();
}

foreach my $i (@conservation) {
	print "$i \n";
}

=cut

#print $alignment->num_sequences. "\n";
foreach my $pos (@columns) {
	my %rescounts;
	my $column = $pos->{'column'};
	#Count the number of each residue in the selected columns
	foreach my $seq ($alignment->each_seq) {
		my $res = $seq->subseq($column, $column);
		$rescounts{$res}++;
	}
	my $count = $rescounts{ $pos->{'residue'} };
	my $gaps = $rescounts{'-'} + $rescounts{'.'};
	my $totalNoGaps = $alignment->num_sequences - $gaps;


	#Calculate relative entropy in the column

	my $entropy;
	my @aminoAcids = {"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"};
	foreach my $aminoacid(@aminoAcids) {
		my $freq = $rescounts{$aminoacid}/$totalNoGaps;
		$entropy = $entropy + $freq*log($freq);
	}

	open OUTPUT, ">>$output" or die "Could not open $output!\n";
	print OUTPUT "$refProtein,$msa,$pos->{'residue'},$pos->{'position'},$pos->{'column'},$count/" . $alignment->num_sequences . ",$count/$totalNoGaps" . "\n";
	close OUTPUT;


#
#	foreach my $res (keys %rescounts) {
#		printf "Res: %s  Count: %2d\n", $res, $rescounts{$res};
#	}
#
#	print "=============================================\n"
}


	

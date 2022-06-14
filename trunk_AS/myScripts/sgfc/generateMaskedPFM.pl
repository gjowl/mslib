#!/usr/bin/perl
# This script will read in a multiple sequence alignment and
# output a position frequency matrix corresponding to just the first sequence


#TODO read in multiple sequence alignment
#TODO accept options for gap mask threshold, etc
#TODO Iterate over alignment columns
#TODO Calculate gap frequency
#TODO Switch case of columns exceeding threshold
#TODO Return masked alignment


use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Align::Utilities qw(:all);
use strict;

########################################################################
#	Options
########################################################################
my $msa = shift;
my $startTime = time;
my $fast = 1;
my $str = Bio::AlignIO->new(	-file => $msa,
				);
my $aln = $str->next_aln();

print $aln->length;
print "\n";
print $aln->num_sequences;
print "\n";
print $aln->num_residues;
print "\n";
print $aln->is_flush;
print "\n";

my $seqLength = $aln->length;
my $numSeqs = $aln->num_sequences;
if (!$aln->is_flush) {
	die "ERROR: ALIGNMENT NOT FLUSH!\n";
}

my %charList;
my @positions;

if ($fast) {
	#Make a two dimensional array of the multiple sequence alignment
	my @msa_chars;
	foreach my $seq($aln->each_seq) {
		my @seq_chars = split('', $seq->seq);
		push (@msa_chars, \@seq_chars);
	}
	
	#Loop over all characters in a column and count gaps
	for (my $col = 0; $col < $seqLength; $col++) {
		my %count;
		for (my $row = 0; $row < $numSeqs; $row++) {
			my $char = $msa_chars[$row][$col];
			$count{$char} += 1.0/$numSeqs;
			$charList{$char} = "Present";
		}
		push (@positions, \%count);


		#for (my $row = 0; $row < $numSeqs; $row++) {
		#	my $mask = $msa_chars[$row][$col];
		#	if ($mask eq "-") {
		#		$mask = ".";
		#	} else {
		#		$mask = lc($mask);
		#	}
		#	$msa_chars[$row][$col] = $mask;
		#}
		#}
	}
	print "\n";

	#Re-concatenate sequences
	#my $row = 0;
	#foreach my $seq($aln->each_seq) {
	#	my $maskedSequence = join ('', @{$msa_chars[$row]});
	#	$seq->seq($maskedSequence);

	#	print ">" . $seq->display_id . "\n";
	#	print $seq->seq . "\n";
	#	$row++;
	#}
	#
	
	#Loop over set of characters in the MSA and print out frequency of that char 
	#for each position
	
	print "POS\t";
	foreach my $character (sort keys %charList) {
		print "$character\t";
	}
	print "\n";

	for (my $i = 0; $i < @positions; $i++) {
		print "$i\t";
		my %posCount = %{$positions[$i]};
		foreach my $character (sort keys %charList) {
			if ($msa_chars[0][$i] eq $character) { #Check if this character is the wild-type residue
				print "1\t";
			} else {
				print "0\t";
			}
		}
		print "\n";
	}

			

}
my $endTime = time;
my $totalTime = $endTime - $startTime;
print "Total time: $totalTime\n";

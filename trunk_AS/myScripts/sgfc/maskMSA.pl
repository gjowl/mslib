#!/usr/bin/perl
# This script will read in a multiple sequence alignment and mask
# columns with a gap count above a certain threshold. A masked 
# column in the MSA will consist of all lower-case letters and .
# for gaps instead of -
# Example input: columns 8-10 will be masked
# ATTSGVVAACG--GTEW
# SSSAGVVW-RGWAGTEW
# SS-AGVV---GWAGTEW
# GTS-GVW---MWAATEW

# Corresponding output:
# ATTSGVVaacG--GTEW
# SSSAGVVw.rGWAGTEW
# SS-AGVV...GWAGTEW
# GTS-GVW...MWAATEW



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
my $threshold = shift;
my $outputFormat = shift;
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
		}
		if ($msa_chars[0][$col] ne '-') {
			printf ("%.2f", $count{'-'});
			print " ";
		}

		if ($count{'-'} > $threshold) {
			#mask the column if there are too many gaps

			for (my $row = 0; $row < $numSeqs; $row++) {
				my $mask = $msa_chars[$row][$col];
				if ($mask eq "-") {
					$mask = ".";
				} else {
					$mask = lc($mask);
				}
				$msa_chars[$row][$col] = $mask;
			}
		}
	}
	print "\n";

	#Re-concatenate sequences
	my $row = 0;
	foreach my $seq($aln->each_seq) {
		my $maskedSequence = join ('', @{$msa_chars[$row]});
		$seq->seq($maskedSequence);

		print ">" . $seq->display_id . " " . $seq->desc . "\n";
		print $seq->seq . "\n";
		$row++;
	}

			

} else {
	my $maskaln = $aln->slice(1,1,1);
	
	#This bit of code suppresses warnings about not being
	#able to guess the sequence alphabet when conatenating
	#the column alignments. We'll get rid of the dummy 
	#column at the end.
	foreach my $seq ($maskaln->each_seq) {
		my $dummySeq = "A";
		my $id = $seq->display_id . " " . $seq->desc;
		$seq->seq($dummySeq, 'protein');
		#correct the end value for LocatableSeq
		$seq->end(1);
	}
	
	#Loop over input alignment
	for (my $i = 1; $i <= $seqLength; $i++) {
		#final 1 keeps the gaps in the column
		my $column = $aln->slice($i, $i, 1);
		
		#Count occurences of each character in the column
		my %count;
		foreach my $seq ($column->each_seq) {
			#$seq->alphabet('protein');
			#print $seq->seq();
			#print "\n";
			my $res = $seq->seq;
			$count{$res} += 1.0/$numSeqs;
		}
		#print "$count{'-'}\n";
		my $fracGaps = $count{'-'} + $count{'.'} + $count{'_'} + $count{'~'};	#Count up all common gap characters
		if ($fracGaps >= $threshold) {
			#mask the column if there are too many gaps
			foreach my $char ($column->each_seq) {
				my $lower = lc($char->seq);
				if ($lower eq "-") {
					$lower = ".";
				}
				$char->seq($lower, 'protein');
	
			}
			#$column->map_chars('\-', '.');
		}
	
		
		$maskaln = cat($maskaln, $column);
	}
	#Get rid of the dummy column we added before
	$maskaln = $maskaln->remove_columns([0,0]);
	
	foreach my $seq ($maskaln->each_seq) {
		print ">" . $seq->display_id . " " . $seq->desc;
		print "\n";
		print $seq->seq;
		print "\n";
	}
}
my $endTime = time;
my $totalTime = $endTime - $startTime;
print "Total time: $totalTime\n";

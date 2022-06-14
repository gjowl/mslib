#!usr/bin/perl
#TODO Parse Data File
#TODO put single and pairwise freqs into matirx
#TODO RANK column pairs
#TODO OUTPUT pairwise freqs
#TODO OUTPUT all pairs
#	'===================' separates pairs above threshold from rest

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Align::Utilities qw(:all);
use strict;
use Getopt::Long;
#CalcMI <dataset> <bins> <th> <out> <tpn>
#Takes as input the expression data for a set of genes over a number of
#trials and reconstructs pairwise gene-gene dependencies using mutual
#information. The program will output the list of gene-gene dependencies 
#and annotate those that pass a specified mutual information threshold as
#positive edge predictions. It should only consider the dependencies 
#between unique genes, not the mutual information of a gene and itself
#(the #entropy of that gene).

#Data file: tab-delimited,  first line shows labels of dataaset columns
#Each column after the first represents expression of a gene during
#different time courses. The data for each trial is separated by a blank
#line.

my $msa = shift;

my @avg_expr;

my $str = Bio::AlignIO->new(	
	-file => $msa,
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

my @aas = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-');
my %aa2Index = (
	"A" => 0,
	"C" => 1,
	"D" => 2,
	"E" => 3, 
	"F" => 4,
	"G" => 5,
	"H" => 6,
	"I" => 7,
	"K" => 8,
	"L" => 9,
	"M" => 10,
	"N" => 11,
	"P" => 12,
	"Q" => 13,
	"R" => 14,
	"S" => 15,
	"T" => 16,
	"V" => 17,
	"W" => 18,
	"Y" => 19,
	"-" => 20,
);
















#Make a two dimensional array of the multiple sequence alignment
my @msa_chars;
foreach my $seq($aln->each_seq) {
	my @seq_chars = split('', $seq->seq);
	push (@msa_chars, \@seq_chars);
}


#Calculate individual frequencies for each column and pairwise for each pair of columns
my @individual_freq;	#n by 20

for (my $col = 0; $col < $seqLength; $col++) {
	
	for (my $row = 0; $row < $numSeqs; $row++) {
		my $char = $msa_chars[$row][$col];
		#my $charIndex = $aa2Index{$char};
		#print "$char $charIndex\n";
		$individual_freq[$col][$aa2Index{$char}] += 1/$numSeqs;

	}
}


my @pairwise_freq;	#

for (my $col_i = 0; $col_i < $seqLength; $col_i++) {
	for (my $col_j = $col_i + 1; $col_j < $seqLength; $col_j++) {
		if ($col_i == $col_j) {
			next;
		}

		for (my $row = 0; $row < $numSeqs; $row++) {
			my $char_i = $msa_chars[$row][$col_i];
			my $char_j = $msa_chars[$row][$col_j];

			$pairwise_freq[$col_i][$col_j][$aa2Index{$char_i}][$aa2Index{$char_j}] += 1/$numSeqs;



		}




	}
}


#calculate mututal information for each pair of columns

for (my $col_i = 0; $col_i < $seqLength; $col_i++) {
	for (my $col_j = $col_i + 1; $col_j < $seqLength; $col_j++) {
		
		if ($col_i == $col_j) {
			next;
		}

		#Loop over all amino acids here and calculate MI
		my $MI = 0;
		for (my $aa_a = 0; $aa_a < @aas; $aa_a++) {
			for (my $aa_b = 0; $aa_b < @aas; $aa_b++) {

				my $freq_ia = $individual_freq[$col_i][$aa_a];
				my $freq_jb = $individual_freq[$col_j][$aa_b];
				my $freq_pair = $pairwise_freq[$col_i][$col_j][$aa_a][$aa_b];
				if ($freq_ia > 0 and $freq_jb > 0 and $freq_pair > 0) {
					$MI += $freq_pair * log($freq_pair/$freq_ia/$freq_jb);
				}
			}
		}
		my $columni = $col_i+1;
		my $columnj = $col_j+1;
		print "$columni $columnj $MI\n";
	}
}


#Print out a sample freq table

#for (my $i = 0; $i < @aas; $i++) {
#	print "$aas[$i]	$individual_freq[3][$i]\n";
#}

#for (my $i = 0; $i < @aas; $i++) {
#	for (my $j = 0; $j < @aas; $j++) {
#		my $freq = $pairwise_freq[18][74][$i][$j];
#		if (!$freq) {
#			print "0	";
#		} else {
#			print "$freq	";
#		}
#
#	}
#	print "\n";
#}
#

















#open DATA, $dataset or die "ERROR: DID NOT OPEN DATASET $dataset!\n";
#my $firstLine = <DATA>;
#my ($time, @genes) = split ('\t', $firstLine);
#my $index = 0;
#my $num_replicates = 1;
#while (my $line = <DATA>) {
#	chomp $line;
#	if ($line =~ /^\s*$/) {
#		$index = 0;
#		$num_replicates++;
#		next;
#	}
#	my ($timepoint, @exprTimePoints) = split ('\t', $line);
#
#	for (my $g = 0; $g < @exprTimePoints; $g++) {
#		$avg_expr[$index][$g] += $exprTimePoints[$g];
#	}
#
#	$index++;
#
#}
#
##Normalize the summed expression values
#for (my $i = 0; $i < @avg_expr; $i++) {
#	for (my $j = 0; $j < @{$avg_expr[$i]}; $j++) {
#		$avg_expr[$i][$j] /= $num_replicates;
#		#print $avg_expr[$i][$j] . "\t";
#	}
#	#print "\n";
#}
#
#
#
#my @bins;
#for (my $i = 1; $i <= $bins; $i++) {
#	my $binUpper = $i/$bins;
#	push (@bins, $binUpper);
#	#print "$binUpper\t";
#}
##print "\n";
#
#
#
##Calculate Marginal and Joint probabilities for the binned data
#
##MI:	
##	--   --              P(m,c)
##	>    >   P(m,c) log P(m)P(c)
##	--m  --c
##	
#
#my $pseudocount = 0.1;
#
#my %pairedMI;
#for (my $i = 0; $i < @genes; $i++) {
#	for (my $j = $i+1; $j < @genes; $j++) {
#		
#		#Don't care about MI of a gene to itself...
#		if ($i == $j) {
#			next;
#		}
#		#print "$i $j\n";
#		my @countMatrix;
#		#now that we looped over all gene pairs, go through the 
#		#expression data and bin it appropriately
#		for (my $t = 0; $t < @avg_expr; $t++) {
#			my $bindex_i = &bin_data($avg_expr[$t][$i], \@bins);
#			my $bindex_j = &bin_data($avg_expr[$t][$j], \@bins);
#			#print "$avg_expr[$t][$i] $avg_expr[$t][$j]\t";
#			#print "$bindex_i $bindex_j\n";
#			$countMatrix[$bindex_i][$bindex_j] += 1;
#		}
#		
#		#Normalize the counts to get probabilities
#		my $countTotal;
#		for (my $x = 0; $x < @bins; $x++) {
#			for (my $y = 0; $y < @bins; $y++) {
#				$countMatrix[$x][$y] += $pseudocount;
#				$countTotal += $countMatrix[$x][$y];
#			}
#		}
#		
#		my @marginal_row = 0;	#Marginal probability for row
#		my @marginal_col = 0;	#Marginal probability for column
#
#		for (my $x = 0; $x < @bins; $x++) {
#			for (my $y = 0; $y < @bins; $y++) {
#				$countMatrix[$x][$y] /= $countTotal;
#				$marginal_row[$x] += $countMatrix[$x][$y];
#				$marginal_col[$y] += $countMatrix[$x][$y];
#				#print "$countMatrix[$x][$y]	";
#			}
#			#print "$marginal_row[$x]\n";
#		}                              
#		foreach my $col (@marginal_col) {
#			#print $col . "\t";
#		}
#		#print "\n";
#
#		#Calculate Mutual Information for this matrix of bin probabilities
#		
#		my $MI;
#		for (my $x = 0; $x < @countMatrix; $x++) {
#			
#			for (my $y = 0; $y < @{$countMatrix[$x]}; $y++) {
#				my $logTerm = $countMatrix[$x][$y]/$marginal_row[$x]/$marginal_col[$y];
#				$logTerm = &log_base(2, $logTerm);
#				$MI += $countMatrix[$x][$y] * $logTerm;
#			}
#		}
#		
#		my $gene_i = $i+1;
#		my $gene_j = $j+1;
#		my $pairID = "($gene_i, $gene_j)";
#		$pairedMI{$pairID} = $MI;
#
#	}
#}
#
#
##Sort the MI Pairs and print them out appropriately
#
#open OUT, ">$out" or die "ERROR: DID NOT OPEN OUTPUT FILE $out\n";
#open TPN, ">$tpn" or die "ERROR: DID NOT OPEN TOTAL POSITIVE NUMBER OUTPUT FILE $tpn\n";
#
#
#my $thresholdFlag = 0;
#my $totalPositive = 0;
#foreach my $genepair (sort { $pairedMI{$b} <=> $pairedMI{$a} } keys %pairedMI) {
#	if ($pairedMI{$genepair} < $th) {
#		if (!$thresholdFlag) {
#			$thresholdFlag = 1;
#			print OUT "===========================================\n";
#		}
#	} else {
#		$totalPositive ++;
#	}
#	print OUT "$genepair	$pairedMI{$genepair}\n";
#}
#print TPN "$totalPositive\n";
#
#
#
#
#
##Function determines if a data point belongs in a particular bin
##Accepts the data point and an array of bin thresholds
##Returns the index of the bin it belongs to
#sub bin_data {
#	my $_point = $_[0];
#	my @_bins = @{$_[1]};
#
#	my $binnedIndex = -1;
#	for (my $i = 0; $i < @_bins; $i++) {
#		if ($_point <= $_bins[$i]) {
#			$binnedIndex = $i;
#			last;
#		}
#	}
#
#	return $binnedIndex;
#}
#
#sub log_base {
#	my ($base, $value) = @_;
#	return log($value)/log($base);
#}

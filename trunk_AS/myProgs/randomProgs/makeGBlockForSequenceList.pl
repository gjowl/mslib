#!/usr/bin/perl
#Edited version of /data07/sanderson/2017JACSPaper/gBlockProject/gBlockSynthesis/gBlockSynthesisProgramGeneral.pl
#Edited by GJL to accept an input file of sequences rather than inputting each sequence as a command line argument
#Date: 2021_10_25

#Run with: perl gBlockSynthesisProgram.pl filename

my $seqFile = shift;
open(FILE,$seqFile);
my @data = <FILE>;
close(FILE);

my @sequenceArray;

for (my $i=0;$i<@data;$i++) {
	$sequence = substr($data[$i],0,length($data[$i])-2);
	push(@sequenceArray, $sequence);
}

my $inputLength = @data+1;

my @STP = ("TAA","TGA","TAG");
$aaToDNAHash{"STP"} = \@STP;
my @ALA = ("GCT","GCC","GCA","GCG");
$aaToDNAHash{"A"} = \@ALA;
my @CYS = ("TGT","TGC");
$aaToDNAHash{"C"} = \@CYS;
my @ASP = ("GAT","GAC");
$aaToDNAHash{"D"} = \@ASP;
my @GLU = ("GAA","GAG");
$aaToDNAHash{"E"} = \@GLU;
my @PHE = ("TTT","TTC");
$aaToDNAHash{"F"} = \@PHE;
my @GLY = ("GGT","GGC","GGA","GGG");
$aaToDNAHash{"G"} = \@GLY;
my @HIS = ("CAT","CAC");
$aaToDNAHash{"H"} = \@HIS;
my @ILE = ("ATT","ATC","ATA");
$aaToDNAHash{"I"} = \@ILE;
my @LYS = ("AAA","AAG");
$aaToDNAHash{"K"} = \@LYS;
my @LEU = ("TTA","TTG","CTT","CTC","CTG");
$aaToDNAHash{"L"} = \@LEU;
my @MET = ("ATG");
$aaToDNAHash{"M"} = \@MET;
my @ASN = ("AAT","AAC");
$aaToDNAHash{"N"} = \@ASN;
my @PRO = ("CCT","CCC","CCA","CCG");
$aaToDNAHash{"P"} = \@PRO;
my @GLN = ("CAA","CAG");
$aaToDNAHash{"Q"} = \@GLN;
my @ARG = ("CGT","CGC","CGG");
$aaToDNAHash{"R"} = \@ARG;
my @SER = ("TCT","TCC","TCA","TCG","AGT","AGC");
$aaToDNAHash{"S"} = \@SER;
my @THR = ("ACT","ACC","ACA","ACG");
$aaToDNAHash{"T"} = \@THR;
my @VAL = ("GTT","GTC","GTA","GTG");
$aaToDNAHash{"V"} = \@VAL;
my @TRP = ("TGG");
$aaToDNAHash{"W"} = \@TRP;
my @TYR = ("TAT","TAC");
$aaToDNAHash{"Y"} = \@TYR;

my %endLeuHash;
my @ALA = ("GCG");
$endLeuHash{"A"} = \@ALA;
my @GLU = ("GAG");
$endLeuHash{"E"} = \@GLU;
my @GLY = ("GGG");
$endLeuHash{"G"} = \@GLY;
my @LYS = ("AAG");
$endLeuHash{"K"} = \@LYS;
my @LEU = ("TTG","CTG");
$endLeuHash{"L"} = \@LEU;
my @MET = ("ATG");
$endLeuHash{"M"} = \@MET;
my @PRO = ("CCG");
$endLeuHash{"P"} = \@PRO;
my @GLN = ("CAG");
$endLeuHash{"Q"} = \@GLN;
my @ARG = ("CGG");
$endLeuHash{"R"} = \@ARG;
my @SER = ("TCG");
$endLeuHash{"S"} = \@SER;
my @THR = ("ACG");
$endLeuHash{"T"} = \@THR;
my @VAL = ("GTG");
$endLeuHash{"V"} = \@VAL;
my @TRP = ("TGG");
$endLeuHash{"W"} = \@TRP;

my @randomBP = ("G","C","T","A");

my $errors = 0;
for(my $i=0; $i<scalar(@sequenceArray); $i++) {
	my $seq = $sequenceArray[$i];

	for(my $j=0; $j<length($seq); $j++) {
		$aa = substr($seq, $j, 1);
		$aa = uc($aa);
		if ($aa !~ m/[ACDEFGHIKLMNPQRSTVWY]/) {
			print "$aa is not a recognized amino acid.\n";
			$errors++;
		}
	}
	my $endAminoAcid = substr($seq, length($seq)-1);
	if ($endAminoAcid =~ m/[CDFHIMYN]/) {
		print "amino acid \"$endAminoAcid\" does not have a codon that ends with G(uanosine)\n";
		$errors++;
	}
}
if ($errors > 0) {
	print "There are currently $errors errors, please correct and resubmit\n";
	exit;
}


my $counter = 0;
my $annotatedSequence = "|5' Random Sequence                                  ||5' amp site     |";
my $gBlockSequence =    "TGTGACGGTTTCCTCGCGCCGGTTGTCCTGTACTGCAAAGGGTCGGGATGGCCCCGAGTCTTCAGCACGTAC";
my $translatedSeq =     "                                                                        ";

for(my $i=0; $i<scalar(@sequenceArray); $i++) {
	$counter++;
	if ($counter > 1000) {
		print "fail, cannot find sequence without repeats.\n";
		exit;
	}
	
	my @inputSeq;
	my $seq = $sequenceArray[$i];
	for(my $j=0; $j<length($seq); $j++) {
		$aa = substr($seq, $j, 1);
		$aa = uc($aa);
		push(@inputSeq, $aa);
	}

	my $currentSeq = "";
	for(my $j=0; $j<scalar(@inputSeq); $j++) {
		my $seqAA = $inputSeq[$j];

		if ($j == scalar(@inputSeq)-1) {
			my $codonSize = scalar(@{$endLeuHash{$seqAA}});
			my $randomNum = int(rand($codonSize));
			my $randomCodon = $endLeuHash{$seqAA}[$randomNum];
			$currentSeq .= $randomCodon;
		}
		else {
		       my $codonSize = scalar(@{$aaToDNAHash{$seqAA}});
		       my $randomNum = int(rand($codonSize));
		       my $randomCodon = $aaToDNAHash{$seqAA}[$randomNum];
		       $currentSeq .= $randomCodon;
		       #print $randomCodon." ";
		}
	}
	if ($currentSeq =~ m/GATC/) {
		#print "$currentSeq\n";
		print "dpnII cut site found, redo sequence\n";
		$i--;
		next;
	}
	if ($currentSeq =~ m/GCTAGC/) {
		#print "$currentSeq\n";
		print "nheI cut site found, redo sequence\n";
		$i--;
		next;
	}

	# Check current sequence for repeats of 8 characters
	#print "$currentSeq\n";
	my $windowSize = 8;
	my $redoSeq = 0;
	my $checklength = length($currentSeq)- $windowSize;	
	if ($checklength > $windowSize) {
		for (my $j=0; $j<$checklength; $j++) {
			for (my $k=$j+$windowSize+1; $k<$checklength; $k++) {
				my $currentFragment = substr($currentSeq, $j, $windowSize);
				my $newFragment = substr($currentSeq, $k, $windowSize);
				if ($newFragment eq $currentFragment) {
					print "seq match found: $newFragment, repeating\n";
					$redoSeq = 1;
					last;
				}
			}
			if ($redoSeq == 1) {
				last;
			}
		}
	}
	if ($redoSeq > 0) {
		$i--;
		next;
	}
	# Check against previously written sequences for repeats
	$redoSeq = 0;
	if (length($gBlockSequence) >= $windowSize) {
		for (my $j=0; $j<$checklength; $j++) {
			my $gBlockLength = length($gBlockSequence) - $windowSize;
			for (my $k=0; $k<$gBlockLength; $k++) {
				my $newFragment = substr($currentSeq, $j, $windowSize);
				my $currentFragment = substr($gBlockSequence, $k, $windowSize);
				if ($newFragment eq $currentFragment) {
					print "seq match found: $newFragment, repeating\n";
					$redoSeq = 1;
					last;
				}
			}
			if ($redoSeq == 1) {
				last;
			}
		}
	}
	if ($redoSeq > 0) {
		$i--;
		next;
	}
	$annotatedSequence .= "|NheI|";
	my $seqNum = $i + 1;
	my $seqNumString = "|seq $seqNum";
	$annotatedSequence .= $seqNumString;
	$annotatedSequence .= (" " x (length($currentSeq) - length($seqNumString) - 2));
	$annotatedSequence .= "|";
	$annotatedSequence .= "|d2|";
	$annotatedSequence .= "|gap |";

	$translatedSeq .= "      "; #NheI cutsite
	for (my $j=0; $j < scalar(@inputSeq); $j++) { $translatedSeq .= "$inputSeq[$j]  "; } #sequence	
	$translatedSeq .= "         "; #dpnII cut site (needs preceeding G) + spacer

		
	$gBlockSequence .= "GCTAGC"; #NheI cutsite
	$gBlockSequence .= $currentSeq; #sequence 
	$gBlockSequence .= "ATC"; #dpnII cut site (needs preceeding G)
	for (my $j=1; $j<=6; $j++) {
		my $randomNum = int(rand(4));
		$gBlockSequence .= $randomBP[$randomNum]; # random spacer
	}
}
$annotatedSequence .= "|3' amp site         ||3' Random Sequence                                 |";
$gBlockSequence .=    "TAAATACTTGGCAGCACCCGAGAGCTCCGCGAGGCGGAATGAGAACGACCAGTTAACCTGTAACCCTCCGGCGTG";

print "\nSequence (no spaces):\n";
print $gBlockSequence."\n";

print "\nAnnotated Sequence:\n";
my $n = 0;
my $lineWidth = 120;
while ($n <= length($gBlockSequence)) {
	print substr($annotatedSequence, $n, $lineWidth)."\n";
	print substr($gBlockSequence, $n, $lineWidth)."\n";
	print substr($translatedSeq, $n, $lineWidth)."\n";
	print "\n";
	$n += $lineWidth;
}


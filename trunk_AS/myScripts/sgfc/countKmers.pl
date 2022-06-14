#!usr/bin/perl
# This script will count the number of unique k-mers present in a set of
# DNA sequences.  To do so, it chops up each sequence into all of its
# k-mers, and then looks them up in a hash table.  At the end, it will
# print out each k-mer and the number of times it appears in a sequence.
# Right now, this program assumes each line in the DNA sequence file
# is one unique sequence.
# usage:
# 	$ perl countKmers.pl <dnaFile> <k>
# 	$ perl countKmers.pl dnaSeqs.txt 12


my $dnaFile = shift;
my $k = shift;

my %kmerCounts;
open SEQFILE, $dnaFile or die "ERROR: did not open $dnaFile!\n";
while (my $dnaSeq = <SEQFILE>) {
	my $seqlength = length($dnaSeq);
	for (my $i = 0; $i < $seqlength - $k; $i++) {
		my $subSeq = substr($dnaSeq, $i, $k);
		# Check if the subsequence is unique
		if ($kmerCounts{$subSeq}) {
			$kmerCounts{$subSeq} += 1;
		} else {
			$kmerCounts{$subSeq} = 1;
		}
	}
}


foreach my $kmer (keys %kmerCounts) {
	print "$kmer	$kmerCounts{$kmer}\n";
}

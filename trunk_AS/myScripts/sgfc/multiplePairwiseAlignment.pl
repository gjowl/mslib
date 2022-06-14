#This script compares a set of sequences to a main sequence using pairwise alignment
#This is useful to compare a bunch of sequences from SDM experiments to the wild-type
#or intended mutation to rapidly screen for the ones you need

use strict;

use Bio::Tools::dpAlign;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Matrix::IO;


my $fastafile = shift;
#my $seqToCompare = shift;

#TODO if there is no seq to compare, use the first sequence in the fasta file



my $bioseq = Bio::SeqIO->new(	-file => $fastaA,
				-format => 'Fasta',
				-delimiter => '\n',
				-alphabet => 'protein');


my $seqToCompare = $bioseq->next_seq;

# create a dpAlign object
# to do global alignment, specify DPALIGN_GLOBAL_MILLER_MYERS
# to do ends-free alignment, specify DPALIGN_ENDSFREE_MILLER_MYERS
$parser = Bio::Matrix::IO->new(-format => 'scoring', -file => 'blosum50.mat');
  $matrix = $parser->next_matrix;

my $factory = new dpAlign(-matrix => $matrix,
		   -match => 3,
                   -mismatch => -1,
                   -gap => 3,
                   -ext => 1,
                   -alg => Bio::Tools::dpAlign::DPALIGN_LOCAL_MILLER_MYERS);

# actually do the alignment
#


# Iterate over sequences
while (my $sequence = $bioseq->next_seq) {
	$alignment = $factory->pairwise_alignment($seqToCompare, $sequence);
	$alnout = Bio::AlignIO->new(-format => 'pfam', -fh => \*STDOUT);
	$alnout->write_aln($out);
}

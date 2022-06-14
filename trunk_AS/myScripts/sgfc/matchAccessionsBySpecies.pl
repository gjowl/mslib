#This script will take a record from a blast FASTA file and parse the species
#information from the record. It will then find the same record in two files
#and pull out those fasta records. This will solve any issues with programs
#stripping out some vital information in the sequence manipulation.
#
#



#!/usr/bin/perl
use Bio::SeqIO;
use strict;

########################################################################
#	Options
########################################################################
my $fasta_merged = shift;
my $fasta_A = shift;
my $fasta_B = shift;

my $source = shift;


########################################################################
#	Subroutines
########################################################################

sub parse_blastFasta_gi {
	#2016-09-26 BLAST CHANGED HEADER FORMAT, this function depreciated and
	#renamed to parse_blastFasta_gi from parse_blastFasta!
	#
	#
	# This function will parse out the species present 
	# in a nr-database record from BLAST. For convenience,
	# only the first species will be used from a 
	# MULTISPECIES record.
	# >gi|xx|dbsrc|accession.version|gene name [organism]
	my $id = $_[0];
	#my ($firstId) = $id =~ /^(>.+?)>?$/;	#Only examine first id record
	my $firstId = $id;
	#my ($gi, $geneNumber, $db, $accession, $description) = split(/\|/, $firstId);
	my ($gi, $geneNumber, $db, $accession, $description) = $firstId =~
		/^
		(.+?)\|
		(.+?)\|
		(.+?)\|
		(.+?)\|
		(.+?)$/x;

	my ($geneName, $organism) = $description =~
		/(.+?)
		\[(.+?)\]/x;
	#print "$id\n";
	#print "$firstId\n";
	#print "$gi, $geneNumber, $db, $accession\n";
	#print "$description\n";
	#print "$geneName\n";
	#print "$print "$organism\n";
	my %parsedId = (
		'fullid' => $id,
		'id' => $firstId,
		'db' => $db,
		'accession' => $accession,
		'geneNumber' => $geneNumber,
		'geneName' => $geneName,
		'organism' => $organism
	);
	#print $parsedId{'organism'} . "\n";

	return \%parsedId;


}
########################################################################
#	Main
########################################################################

### Open up the fasta files

my $bioseq_merged = Bio::SeqIO->new ( 	-file => $fasta_merged,
					-format => 'Fasta',
					-delimiter => '\n');
my @merged_ids;
my @merged_seqs;

my $bioseq_A = Bio::SeqIO->new(	-file => $fasta_A,
				-format => 'Fasta',
				-delimiter => '\n');
my @fasta_A_ids;
my @fasta_A_seqs;

my $bioseq_B = Bio::SeqIO->new(	-file => $fasta_B,
				-format => 'Fasta',
				-delimiter => '\n');

my @fasta_B_ids;
my @fasta_B_seqs;


### Read the files using BIOSEQ and place sequences into arrays
while ( my $seq = $bioseq_merged->next_seq() ) {
	push(@merged_seqs, $seq);
}
while ( my $seq = $bioseq_A->next_seq() ) {
	push(@fasta_A_seqs, $seq);
}
while ( my $seq = $bioseq_B->next_seq() ) {
	push(@fasta_B_seqs, $seq);
}

### Read the files again to parse the sequence id info
open MERGED, $fasta_merged or die "ERROR: DID NOT OPEN $fasta_merged\n";
open FA_A, $fasta_A or die "ERROR: DID NOT OPEN $fasta_A\n";
open FA_B, $fasta_B or die "ERROR: DID NOT OPEN $fasta_B\n";

while (my $line = <MERGED>) {
	if ($line =~ /^>/) {
		chomp $line;
		my $parsedId = &parse_blastFasta_gi($line);
		push (@merged_ids, $parsedId);
	}
}
while (my $line = <FA_A>) {
	if ($line =~ /^>/) {
		chomp $line;
		my $parsedId = &parse_blastFasta_gi($line);
		push (@fasta_A_ids, $parsedId);
	}
}
while (my $line = <FA_B>) {
	if ($line =~ /^>/) {
		chomp $line;
		my $parsedId = &parse_blastFasta_gi($line);
		push (@fasta_B_ids, $parsedId);
	}
}


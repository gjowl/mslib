# This script accepts two fasta files corresponding to two gene families
# It will generate two ordered lists corresponding to the set of species
# that have genes in both files.  If there are multiple instances of a
# species in a fasta file, it will only use the first one in the list.
# The fasta file must be obtained from UNIPROT or in the same format:
#
#
# >tr|A0A084C6S7|A0A084C6S7_9PSED Cell division protein FtsL OS=Pseudomonas sp. WCS358 GN=ftsL PE=3 SV=1
#
# The OS=<SPECIES> will be used to explicitly match species between the two files
# 
# TODO Read files
# TODO Parse files (FASTA FORMAT)
# 	TODO Capture species identification
# TODO Match species between both sets
# TODO print out set of species found in both sets
# TODO write a fasta file for each set of species
#

#!usr/bin/perl
use Bio::SeqIO;
use strict;

########################################################################
#	Options
########################################################################
my $fastaA = shift;
my $fastaB = shift;
my $source = shift;


# Open file and put sequences into BioPerl object
# May have to separately parse species ID from gene ID...

 
sub parse_uniprotFasta {
	# Subroutine to parse important info from the uniprot fasta ID
	# http://www.uniprot.org/help/fasta-header	
	# >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName[ GN=GeneName]PE=ProteinExistence SV=SequenceVersion
	# Where:
	# *db is 'sp' for UniProtKB/Swiss-Prot and 'tr' for UniProtKB/TrEMBL.
	# *UniqueIdentifier is the primary accession number of the UniProtKB entry.
	# *EntryName is the entry name of the UniProtKB entry.
	# *ProteinName is the recommended name of the UniProtKB entry as annotated in the RecName field. For UniProtKB/TrEMBL entries without a RecName field, the SubName field is used. In case of multiple SubNames, the first one is used. The 'precursor' attribute is excluded, 'Fragment' is included with the name if applicable.
	# *OrganismName is the scientific name of the organism of the UniProtKB entry.
	# *GeneName is the first gene name of the UniProtKB entry. If there is no gene name, OrderedLocusName or ORFname, the GN field is not listed.
	# *ProteinExistence is the numerical value describing the evidence for the existence of the protein.
	# *SequenceVersion is the version number of the sequence.
	
	my $id = $_[0];
	my ($db) = $id =~ />(.+?)\|/;
	my ($accession) = $id =~ /\|(.+?)\|/;
	my ($entryName, $proteinName) = $id =~ /\|.*\|(\w+) (.+?) OS=/;
	my ($organism, $geneName, $proteinExistence, $seqVersion) = $id =~ 
		/OS=(.+?)
		 GN=(.+?)
		 PE=(.+?)
		 SV=(.+?)
		/x;
	
	my %parsedId = (
		'id' => $id,
		'db' => $db,
		'accession' => $accession,
		'entryName' => $entryName,
		'proteinName' => $proteinName,
		'organism' => $organism,
		'geneName' => $geneName,
		'PE' => $proteinExistence,
		'SV' => $seqVersion
	);
	return \%parsedId;
}
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

sub parse_blastFasta {
	#Parser for BLAST header in FASTA line
	#Updated to accomodate removal of GI information
	#https://www.ncbi.nlm.nih.gov/news/03-02-2016-phase-out-of-GI-numbers/
	#Example:
	#	>XP_017825711.1 PREDICTED: glycophorin-A isoform X2 [Callithrix jacchus]
	my $id = $_[0];

	my $firstId = $id;
	my ($accession, $geneName, $organism) = $firstId =~
		/^
		(.+?)\s
		(.+?)\s
		\[(.+?)\]$/x;

	my ($db, $geneNumber) = split (/_/, $accession);

	my %parsedId = (
		'fullid' => $id,
		'id' => $firstId,
		'db' => $db,
		'accession' => $accession,
		'geneNumber' => $geneNumber,
		'geneName' => $geneName,
		'organism' => $organism
	);
	#print "$id\n";
	#print "$accession \t $geneName \t $organism\n";

	return \%parsedId;
}


sub parse_multiSpecies {
	my $_idList = $_[0];
	my @ids = split('>', $_idList);
	
	my @parsedIds;
	foreach my $id (@ids) {
		#print "$id\n";
		my $parsedId = &parse_blastFasta($id);
		push(@parsedIds, $parsedId);
		print $parsedIds[0]->{'organism'} . "\n";
	}
	#print "$parsedIds[0]->{'organism'}" . "\n";
	return $parsedIds[0];
}

sub find_uniqueSpecies {
	# Subroutine to find set of genes with unique species ID.
	# If multiple genes with same species present, only first
	# is retained.
	my $sequences = $_[0];
	my $ids = $_[1];
	
	my %unique;
	my $counter = 0;
	for (my $i = 0; $i < @{$sequences}; $i++) {
		my $fullID = $ids->[$i]->{'id'};
		my $organism = $ids->[$i]->{'organism'};
		if (!exists $unique{$organism}) {
			$counter++;
			$unique{$organism} = $sequences->[$i];
			$sequences->[$i]->display_id($fullID);
			my $seq = $sequences->[$i]->seq();
		}
	}
	return \%unique;

}

open FASTA_A, $fastaA or die "ERROR: DID NOT OPEN THE 1st FASTA FILE!\n";
open FASTA_B, $fastaB or die "ERROR: DID NOT OPEN THE 2nd FASTA FILE!\n";

my $bioseq1 = Bio::SeqIO->new(	-file => $fastaA,
				-format => 'Fasta',
				-delimiter => '\n');
my $bioseq2 = Bio::SeqIO->new(-file => $fastaB);

my $seqOut = Bio::SeqIO->new(
				-format => 'fasta',
				-fh => \*STDOUT);

my @seq1_array;
my @seq1_ids;

my @seq2_array;
my @seq2_ids;

while ( my $seq = $bioseq1->next_seq() ) {
	push(@seq1_array, $seq);
}

while ( my $seq = $bioseq2->next_seq() ) {
	push(@seq2_array, $seq);
}

while (my $line = <FASTA_A>) {
	if ($line =~ /^>/) {
		
		chomp $line;
		$line =~ s/^>//;
		my @ids = split('>', $line);
		my $parsedId;
		if ($source eq 'uniprot') {
			$parsedId = &parse_uniprotFasta($line);
		} elsif ($source eq 'blast') {
			$parsedId = &parse_blastFasta($ids[0]);
			#$parsedId = &parse_multiSpecies($line);
		} elsif ($source eq 'blast_gi') {
			$parsedId = &parse_blastFasta_gi($ids[0]);
		} else {
			die "ERROR: SOURCE NOT RECOGNIZED\n";
		}
		push (@seq1_ids, $parsedId);
	}
}

while (my $line = <FASTA_B>) {
	if ($line =~ /^>/) {
		chomp $line;
		$line =~ s/^>//;
		my @ids = split('>', $line);		
		my $parsedId;
		if ($source eq 'uniprot') {
			$parsedId = &parse_uniprotFasta($line);
		} elsif ($source eq 'blast_gi') {
			$parsedId = &parse_blastFasta_gi($ids[0]);
		} elsif ($source eq 'blast') {
			$parsedId = &parse_blastFasta($ids[0]);
		} else {
			die "ERROR: SOURCE NOT RECOGNIZED\n";
		}
		push (@seq2_ids, $parsedId);
	}
}

my $unique1 = &find_uniqueSpecies(\@seq1_array, \@seq1_ids);
# print "$_\n" for keys %{$unique1};
my $unique2 = &find_uniqueSpecies(\@seq2_array, \@seq2_ids);
#print "$_\n" for keys %{$unique2};

my %sharedUnique;
my $sharedcounter = 0;
my $leftovercounter = 0;

open SHARED1, ">/tmp/shared1.fasta" or die "ERROR: DID NOT OPEN FILE\n";
open SHARED2, ">/tmp/shared2.fasta" or die "ERROR: DID NOT OPEN FILE\n";
open CONCATENATED, ">/tmp/concatenated.fasta" or die "ERROR: DID NOT OPEN FILE\n";
foreach my $key (keys %{$unique1}) {
	if (exists $unique2->{$key}) {
		print "$key\n";
		$sharedcounter++;
		my $seq1id = $unique1->{$key}->display_id;
		my $sequence1 = $unique1->{$key}->seq;

		my $seq2id = $unique2->{$key}->display_id;
		my $sequence2 = $unique2->{$key}->seq;
		print SHARED1 ">$seq1id\n";
		print SHARED1 "$sequence1\n";
		print SHARED2 ">$seq2id\n";
		print SHARED2 "$sequence2\n";
		print CONCATENATED ">$seq1id >$seq2id\n";
		print CONCATENATED "$sequence1$sequence2\n";
		#print "\n";
	} else {
		#print "NOPE	$key\n";
		$leftovercounter++;
	}
}
my $hash1Size = keys %{$unique1};
my $hash2Size = keys %{$unique2};
my $seq1Size = scalar (@seq1_array);
my $seq2Size = scalar (@seq2_array);
print "Total sequences in FASTA file 1: $seq1Size\n";
print "Total sequences in FASTA file 2: $seq2Size\n";
print "Unique Species in FASTA file 1: $hash1Size\n";
print "Unique Species in FASTA file 2: $hash2Size\n";
print "Shared Unique: $sharedcounter\n";
print "Leftover in FASTA file 1 but not 2: $leftovercounter\n";

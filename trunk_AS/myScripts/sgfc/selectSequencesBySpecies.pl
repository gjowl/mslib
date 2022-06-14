#!usr/bin/perl
# This script will accept a file of blast-formatted
# FASTA sequences and a list of accepted species
# (such as that obtained from a text tree result of
# a CommonTree search).  It will return the subset of
# sequences that match the list of species.  This is
# useful for limiting the results of a BLAST search
# to specific phyla, like splitting a BLAST search
# of FtsB homologues between the proteobacteria
# and the Firmicutes (like B. subtilis).
#


use Bio::SeqIO;
use strict;

########################################################################
#	Options
########################################################################
my $fastafile = shift;
my $speciesfile = shift;


########################################################################
#	Subroutines
########################################################################
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
		\[(.+?)\]\s+$/x;

	#print "$firstId\n";
	#print "Parsed: $accession \t $geneName \t $organism\n";
	
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

	return \%parsedId;
}

sub parse_blastFasta_gi {
	# This function will parse out the species present 
	# in a nr-database record from BLAST. For convenience,
	# only the first species will be used from a 
	# MULTISPECIES record.
	# >gi|xx|dbsrc|accession.version|gene name [organism]
	my $id = $_[0];
	my ($firstId) = $id =~ /^(>.+?)>?$/;	#Only examine first id record
	#my ($gi, $geneNumber, $db, $accession, $description) = split(/\|/, $firstId);
	my ($gi, $geneNumber, $db, $accession, $description) = $firstId =~
		/^>
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

	return \%parsedId;


}
sub parse_species {
	# This function will parse the species list.
	# It may be either a simple text file or a
	# tree format as shown below.  Removing the
	# leading spaces and "+" signs will take care
	# of both cases.
	# Example tree output:
	# + Bacteria
	# + + Caldiserica
	# + + + Caldisericum exile
	# + + Synergistetes
	# + + + Anaerobaculum hydrogeniformans
	# + + Chlamydiae/Verrucomicrobia group
	# + + + Verrucomicrobia
	# + + + + Verrucomicrobiae
	my $species = $_[0];
	chomp $species;
	$species =~ s/^[\s\+]+//g; #Remove all leading whitespace or + signs
	#$species =~ s/^[\s]+//g;
	return $species;
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
		my $fullID = $ids->[$i]->{'fullid'};
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


########################################################################
#	Main
########################################################################
my @seqIds;
my @sequences;
my %speciesSubset;
open FASTA, $fastafile or die "ERROR: did not open $fastafile\n";
my $bioseq = Bio::SeqIO->new(	-file => $fastafile,
				-format => 'Fasta',
				-delimiter => '\n');
			

while ( my $seq = $bioseq->next_seq() ) {
	push(@sequences, $seq);
}



while (my $line = <FASTA>) {
	if ($line =~ /^>/) {
		chomp $line;
		$line =~ s/^>//;
		my @ids = split('>', $line);
		my $parsedId = &parse_blastFasta($ids[0]);
		$parsedId->{'fullid'} = $line;
		push (@seqIds, $parsedId);
	}
}
open SPECIES, $speciesfile or die "ERROR: did not open $speciesfile\n";
while (my $line = <SPECIES>) {
	chomp($line);
	my $species = &parse_species($line);
	$speciesSubset{$species} = "Subset";
}

for (my $i = 0; $i < @seqIds; $i++) {
	my $seqOrganism = $seqIds[$i]->{'organism'};
	my $seqFullId = $seqIds[$i]->{'fullid'};
	if (exists $speciesSubset{$seqOrganism}) {
		$sequences[$i]->display_id($seqFullId);
		my $seq = $sequences[$i]->seq();
		print ">$seqFullId\n";
		print "$seq\n";
	}
}


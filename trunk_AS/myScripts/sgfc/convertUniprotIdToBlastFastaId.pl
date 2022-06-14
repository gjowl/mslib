# This program will accept a set of sequences in FASTA format
# It will take the sequence ID and search the uniprot database
# for that particular accession number.
# It will then fetch the RefSeq Accession, the gene name, and 
# the organism. With this information it will then reprint the
# sequence file with a sequence ID like what you would find from
# a BLAST search of the RefSeq Database.


#!/usr/bin/perl

use strict;
use LWP::Simple;
use Bio::SeqIO;
use Getopt::Long;

########################################################################
#	Options
########################################################################
my $fastafile = shift;
my $outfile = shift;

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

sub parse_UniprotFasta {
	#Parser for FASTA output of HMMR {REF}
	#Example: 
	#	>Q32AA6_SHIDS
	my $id = $_[0];
	my $firstId = $id;
	my ($uniprot, $speciesAbbreviation) = $firstId =~
		/^
		(.+?)_
		(.+?)$
		/x;
	print "PARSED: $uniprot \t $speciesAbbreviation\n";
	return $uniprot;


}


########################################################################
#	Main
########################################################################
my @uniprotIds;
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
		my $uniprotId = &parse_UniprotFasta($ids[0]);
		push (@uniprotIds, $uniprotId);
	}
}


#Set up LWP
my $ua = LWP::UserAgent->new;
$ua->show_progress(1);
$ua->env_proxy;

#@uniprotIds[0] = "B8GX61";
#@uniprotIds[1] = "P0AEN4";
for (my $i = 0; $i < @uniprotIds; $i++) {
	
	my $uniprotAccession = $uniprotIds[$i];
	my $website = "http://www.uniprot.org/uniprot/$uniprotAccession.txt";
	
	my $response = $ua->get( $website );
	die $response->status_line, "\n" unless $response->is_success;

	my $content = $response->content;
	my ($dbId)     = $content =~ /DR   RefSeq; (.+);/;
	my ($geneName) = $content =~ /DE   RecName: Full=(.+)/;
	my ($organism) = $content =~ /OS\s+(.+)/;
	
	
	print "$dbId\t$geneName\t$organism\n";

	#Fetch the website
	#Look for the species
	#Look for the database
	#Look for the gene name
}

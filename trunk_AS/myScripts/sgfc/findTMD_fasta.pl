#! /usr/bin/perl
# This program reads a FASTA file and prints out a submit file for CATM.
# To get the membrane start and end numbers, it uses a second file where each line
# is of the format:
#
# <NAME> <tmStartNum> <tmEndNum>
#
#This version accepts custom TM start and end sites.
#NOTE: Much of the information in this program is hard coded, and works off of the parameters 
#specified in the CATM HUMAN_UNIPROT_20131120 run.
#
use lib "/exports/home/scondon/perl5/lib/perl5/";
use LWP::Simple;
use Bio::SeqIO;
use Bio::SearchIO;

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

my $fastaFile = shift;
my $tmFile = shift;
my ($tmFileName, $tmExtension) = split (/\./, $tmFile);
my ($fastaTM, $extension) = split(/\./, $fastaFile);
$fastaTM = $fastaTM . "_TMH.fasta";
# Load the fasta file and parse the sequences using bioperl.  See 
# http://www.bioperl.org/wiki/HOWTO:SeqIO

my $seq_in = Bio::SeqIO->new(
				-format => 'fasta',
				-file 	=> $fastaFile,
			    );


my @seq_array;
while ( my $seq = $seq_in->next_seq() ) {
	push(@seq_array, $seq);
	my $seqId = $seq->display_id;
}

my @seq_ids;
open FASTA, "$fastaFile" or die "ERROR: DID NOT OPEN FASTA FILE\n";
while (my $line = <FASTA>) {
	if ($line =~ /^>/) {
		chomp $line;
		my $parsedId = &parse_uniprotFasta($line);
		push (@seq_ids, $parsedId);
	}
}

open PHOBIUS, "$tmFile" or die "ERROR: DID NOT OPEN PHOBIUS FILE\n";
my $counter = 0;

open TMDFILE, ">/tmp/TMD.fasta" or die "ERROR: DID NOT OPEN TMD FILE\n";
while (my $line = <PHOBIUS>) {
	chomp $line;
	my @tmData = split(/\s+/, $line);
	my $tmTopology = $tmData[3];
	my $Id = $tmData[0];
	my $fastaId = $seq_array[$counter]->display_id;
	#Includes signal peptides
	if ($Id eq $fastaId) {
		if ($tmTopology =~ m/[onci][0-9]+-[0-9]+[onci]/) {
	
			my @tmRange = $tmTopology =~ /[onci][0-9]+-[0-9]+[onci]/g;		
			my ($tmStartNum, $tmEndNum) = $tmRange[0] =~ /[0-9]+/g;
			
	
			#Capture full TM +- 8 residues, if possible	
			$tmStartNum -= 8;
			$tmEndNum += 8;
	
			if ($tmStartNum < 1) {
				$tmStartNum = 1;
			}
			if ($tmEndNum > $seq_array[$counter]->length()) {
				$tmEndNum = $seq_array[$counter]->length();
			}
			my $tmdseq = $seq_array[$counter]->subseq($tmStartNum, $tmEndNum);
	
			print TMDFILE "$seq_ids[$counter]->{'id'}\n";
			print TMDFILE "$tmdseq\n";
		}
		else {
		print "$fastaId\t $Id\n";
		}	
	}
	$counter++;


	

}

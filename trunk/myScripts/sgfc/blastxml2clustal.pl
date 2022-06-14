#!usr/bin/perl
# This program will parse an xml file produced from a BLAST search and send
# the sequences to the Multiple Sequence Alignment program clustalo to 
# generate an MSA.
#

use strict;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::AlignIO;
#use option parser
use Getopt::Long;

########################################################################
#	Options
########################################################################
my $file = '';
my $tmStart = '';
my $tmEnd = '';

GetOptions (
	'file' => \$file,
	'tmStart' => \$tmStart,
	'$tmEnd' => \$tmEnd
);
#Program steps:
#read the xml file
#convert xml to list of fasta sequences
#send sequences to clustalo
#capture clustalo output



my $blastin = new Bio::SearchIO (	-format => 'blastxml',
					-file => $file
				);

my $seqout = Bio::SeqIO->new (		-format => 'fasta');

my $alignIO = Bio::AlignIO->new	(	-format => 'fasta',
				);

my @alignSeq;	
open FASTA, ">blastxml2clustal.fas";
while (my $result = $blastin->next_result) {
	## $result is a Bio::Search::Result::ResultI compliant object
	while( my $hit = $result->next_hit ) {
	## $hit is a Bio::Search::Hit::HitI compliant object
		while( my $hsp = $hit->next_hsp ) {
		## $hsp is a Bio::Search::HSP::HSPI compliant object
			if( $hsp->length('total') > 50 ) {
				 if ( $hsp->start('query') < $tmStart && $hsp->end('query') > $tmEnd) {
				 	#print "Query=",   $result->query_name,
					# " Hit=",        $hit->name,
					# " Length=",     $hsp->length('total'),
					# " Percent_id=", $hsp->percent_identity, "\n",
					# "Sequence =" , $hsp->hit_string, "\n\n";
					# print "> ", $hit->name, "Range = ", 
					# $hsp->start('hit'),"-", $hsp->end('hit'),"\n",
					# $hsp->hit_string, "\n";

					print FASTA ">" . $hit->name . $hsp->start('hit') . "-" . $hsp->end('hit') . "\n" . 
					$hsp->hit_string . "\n";

					}
			}
		}  
	 }
}
close FASTA;

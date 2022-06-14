#! /usr/bin/perl
#This program accepts a UNIPROT ID and returns the sequence in FASTA format


use strict;
use LWP::Simple;
use Bio::SeqIO;
use Getopt::Long;
########################################################################
#	Options
########################################################################

my $uniprotList = shift;

########################################################################
#	Subroutines
########################################################################
sub getUniprotSequence {
	my $uniprotId = $_[0];
	#Set up LWP
	my $ua = LWP::UserAgent->new;
	$ua->show_progress(1);
	$ua->env_proxy;
	
	
	my $website = "http://www.uniprot.org/uniprot/$uniprotId.txt";
	my $response = $ua->get( $website );
	die $response->status_line, "\n" unless $response->is_success;
	
	my $content = $response->content;
	
	my ($dbId)     = $content =~ /DR   RefSeq; (.+);/;
	my ($geneName) = $content =~ /DE   ...Name: Full=(.+)/;
	my ($organism) = $content =~ /OS\s+(.+)/;
	my $seqId = ">$dbId\t$geneName\t$organism";
	
	#Get the sequence
	my @data = split(/\n/, $content);
	my $sequence = "";
	for(my($i) = 0; $i < @data; $i++) {
		if(substr($data[$i],0,2) eq "SQ") {
			while(1) {
				$i++;
				if(substr($data[$i],0,2) eq "//") {
					last;
				}
				chomp($data[$i]);
				$data[$i] =~ s/\s+//g;
				$sequence = $sequence . $data[$i];
			}
		}
	}
	
	my $fasta =  "$seqId\n$sequence";
	#print "$fasta\n";
	return $fasta;
}

########################################################################
#	Main
########################################################################

open UNIPROTLIST, $uniprotList or die "ERROR: DID NOT OPEN $uniprotList!\n";
while (my $line = <UNIPROTLIST>) {
	chomp $line;
	my $uniprotSeq = &getUniprotSequence($line);
	print "$uniprotSeq\n";
}

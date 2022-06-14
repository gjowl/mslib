#!usr/bin/perl
#
#This program will accept a UNIPROT ID along with start and end positions.  With
#this information, it will access the Pfam webpage for that protein.  Then it will
#fetch the Pfam stockholm MSAs of any matching domains that contain the start and
#end sites of interest.
#
# Read in ID, start, and end
# Fetch webpage for that ID
# Parse the webpage for domain hits
# Fetch the stockholm alignment, reporting an error if no alignment found
#

use strict;
use warnings;
use LWP::Simple;
use LWP::UserAgent;
#use Bio::SimpleAlign;
use Getopt::Long;
use XML::Simple;	#Tutorial http://www.techrepublic.com/article/parsing-xml-documents-with-perls-xmlsimple/
use Data::Dumper;

#PFam puts up alignments in stockholm format as plain text webpages, so it should be possible to grab it from online
#and save it without using a complex interface.
#Example of webpage for PFam A: http://pfam.xfam.org/family/PF01102/alignment/full
#Example of webpage for PFam B: http://pfam.xfam.org/pfamb/PB001661/alignment/format?format=stockholm&gaps=default&download=0
#
#
#==================================================================
#	OPTIONS
#==================================================================

my $uniprot = '';
my $start = '';
my $end = '';
my $outputDir = '';
my $logfile = "logfile.txt";

GetOptions (
	'uniprot=s' => \$uniprot,
	'start=i' => \$start,
	'end=i' => \$end,
	'output=s' => \$outputDir,
	'logfile=s' => \$logfile
);

my $ua = LWP::UserAgent->new;
$ua->show_progress(1);
$ua->env_proxy;



my $website = "http://pfam.xfam.org/protein/" . $uniprot . "?output=xml";


my $res = $ua->get( $website );

die $res->status_line, "\n" unless $res->is_success;

my $xml = $res->content;
my $parser = XML::Simple->new;

my $data = $parser->XMLin($xml, KeyAttr => []);#KeyAttr defaults to 'name', 'key', 'id' to fold arrays into hashes unless you give it an empty list and it's driving me nuts
my $domains = $data->{'entry'}->{'matches'}->{'match'};

#Begin parsing the document for domains that fit the criteria
#If the query has multiple matches, they will be stored in an array.  Otherwise the
#single match will be stored in a single hash.
my @hits;
if (ref($domains) eq "ARRAY") {
	foreach my $domain (@{$domains} ) {
		my $location = $domain->{'location'};
		if ($location->{'start'} <= $start and $location->{'end'} >= $end) {
			my %attributes;
			$attributes{'accession'} = $domain->{'accession'};
			$attributes{'type'} = $domain->{'type'};
			$attributes{'uniprot'} = $uniprot;
			push @hits, \%attributes;
			
		}
	}
}

elsif (ref($domains) eq "HASH") {
	my $location = $domains->{'location'};
	if ($location->{'start'} < $start && $location->{'end'} > $end) {
		my %attributes;
		$attributes{'accession'} = $domains->{'accession'};
		$attributes{'type'} = $domains->{'type'};
		$attributes{'uniprot'} = $uniprot;
		push @hits, \%attributes;
		
	}
}
if (!@hits) {
	print "Did not find a domain for $uniprot between $start and $end...\n";
	open LOGFILE, ">>$outputDir/$logfile\n" or die "Cannot open logfile!";
	print LOGFILE "$uniprot,0\n";
	close LOGFILE;
}
foreach my $hit (@hits) {
	my $stockholm;

	#Example of webpage for PFam A: http://pfam.xfam.org/family/PF01102/alignment/full
	#Example of webpage for PFam B: http://pfam.xfam.org/pfamb/PB001661/alignment/format?format=stockholm&gaps=default&download=0
	my $msaWeb;
	if ($hit->{'type'} eq "Pfam-A") {
		$msaWeb = "http://pfam.xfam.org/family/" . $hit->{'accession'} . "/alignment/full";
	} elsif ($hit->{'type'} eq "Pfam-B") {
		$msaWeb = "http://pfam.xfam.org/pfamb/" . $hit->{'accession'} . "/alignment/format?format=stockholm&gaps=default&download=0";
	}
	
	print "$msaWeb\n";


	my $msaPage = $ua->get( $msaWeb );
#	die $msaPage->status_line, "\n" unless $res->is_success;

	$stockholm = $msaPage->content;

	if ($stockholm =~ /^# STOCKHOLM/) {
		open LOGFILE, ">>$outputDir/$logfile\n" or die "Cannot open logfile!";
		print LOGFILE "$hit->{'uniprot'},$hit->{'accession'}\n";
		close LOGFILE;

		my $outputFile;
		if ($outputDir) {
			$outputFile = "$outputDir/$hit->{'accession'}.stk";
			open STOCKHOLM, ">$outputFile";
			print STOCKHOLM $stockholm;
			close STOCKHOLM;
	
		} else {
			print $stockholm;
		}
	} else {
		open LOGFILE, ">>$outputDir/$logfile\n";
		print LOGFILE "$hit->{'uniprot'},$hit->{'accession'},NOTFOUND\n";
		close LOGFILE;
	}
}	



#==============SUBROUTINES=============================
#



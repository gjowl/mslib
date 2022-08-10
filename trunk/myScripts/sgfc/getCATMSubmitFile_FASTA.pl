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


my $seq;
my @seq_array;
while ( $seq = $seq_in->next_seq() ) {
	push(@seq_array, $seq);
}

# Load the TMHMM output file (one line per protein)
# This file has the following tab-separated format
#
# <Fasta ID>	<length>	<# of TM amino acids>	<# of TM AAs in first 60 res>	<# of helices>	<Membrane topology>
#
# gi|441619185|ref|XP_003258060.2|	len=149	ExpAA=44.69	First60=22.12	PredHel=1	Topology=i90-112o
#
# This information will be put into a hash, where the Fasta ID will be the key used to match it up with the sequence data in
# the fasta file.  Later, we will take the rest of the data and parse it to get the individual helices.
#print "$tmFile\n";
open(TMFILE, $tmFile);
my %tmHash;
while (my $line = <TMFILE>) {
	#print $line;
	chomp $line;
	my ($tmId, @tmData) = split (/\s+/, $line);
	@{$tmHash{$tmId}} = @tmData;
}
=pod
foreach my $key ( keys %tmHash ) {
	print "$key\n";
	foreach ( @{$tmHash{$key}} ) {
		print "$_\n";
	}
	print "\n"
}
=cut

my @catmSubmit;
my @mkdirSubmit;
my @tmSequences;

foreach my $seq (@seq_array){


	
#	my ($uniprotId, $tmStartNum, $tmEndNum) = split (/\s+/,$protein);
#	getstore("http://www.uniprot.org/uniprot/$uniprotId.txt","/data00/tmp/$uniprotId.txt");
#	open(FILE,"/data00/tmp/$uniprotId.txt");
#	my @data = <FILE>;
#	close(FILE);
	my ($fastaId) = $seq->display_id;


	my (@tmData) = split(/\s+/, $tmHash{$fastaId});

	#This program can differentiate between files from TMHMM and Phobius if the extensions are set correctly.
	#Phobius Format: SEQENCE ID                     TM SP PREDICTION
	#VarNames:	 $tmId                          0  1  2
	# Check to see if the protein is bitopic--PredHel=1
	my $numberTm;
	my $tmTopology;
	if ($tmExtension eq "phobius") {
		$numberTM = $tmHash{$fastaId}[0];
		$tmTopology = $tmHash{$fastaId}[2];
	} elsif ($tmExtension eq "tmhmm") {
		$numberTM = $tmHash{$fastaId}[3];
		$tmTopology = $tmHash{$fastaId}[4];
	} else {
		die "ERROR: COULD NOT DETERMINE TRANSMEMBRANE PREDICTION \n";
	}
	#print "$numberTM $tmTopology\n";
	if ($numberTM =~/1/) {
	#print "$tmHash{$fastaId}[3] $tmTopology\n";
		# We're interested in the "membrane topology" item of the TMHMM file
		my @tmRange = $tmTopology =~ /[oi][0-9]+-[0-9]+[oi]/g;
		for (my $i = 0; $i < @tmRange; ++$i ) {
			my $fastaTmId = $fastaId;
			if (scalar @tmRange > 1) {
				$fastaTmId = $fastaTmId . "_TM" . ($i+1);
			} 
			#my ($tmStartNum, $tmEndNum) = split (/-/, $tmRange[$i]);
			my ($tmStartNum, $tmEndNum) = $tmRange[$i] =~ /[0-9]+/g;
	
			my($currentId) = "NULL";
			my($uniprotName) = "NULL";
			my($detail);
			my($catmStartNum);
			my($catmEndNum);
			my($ext) = 3; #default to 3 residues on either side of the TM
			my($startExt) = $ext;
			my($endExt) = $ext;
			my($sequence) = $seq->seq;
	
			my($len) = $tmEndNum - $tmStartNum + 1;
			while($tmStartNum <= $startExt) {
				$startExt--;
			}
			while(length($sequence) - $tmEndNum < $endExt) {
				$endExt--;
			}
			$catmStartNum = $tmStartNum - $startExt;
			$catmEndNum = $catmStartNum + $len + $endExt + $startExt - 1;
			my($tmSequence) = substr($sequence,$tmStartNum-1-$startExt,$len + $startExt + $endExt);
	
			my($threadStart) = 35 - (($catmEndNum - $catmStartNum + 1) - 3 - $startExt);
			my($threadEnd) = 33 - $endExt;
	
	
			push(@catmSubmit, "pbsq --qsub \"/data03/CATM/bin/CATM_v21 --backboneCrd /data03/CATM/files/69-gly-residue-helix.crd --logFile ./$fastaTmId/$fastaTmId.log --helixGeoFile /data03/CATM/files/CENTROIDS_0.5_cutoff-4_byHbond --fullSequence $sequence --topFile /data03/CATM/topparsolvhbond/top_all22_prot.inp --parFile /data03/CATM/topparsolvhbond/par_all22_prot.inp --solvFile /data03/CATM/topparsolvhbond/solvpar22.inp --hBondFile /data03/CATM/topparsolvhbond/par_hbond_CA.inp --rotLibFile /data03/CATM/files/EBL_11-2011_CHARMM22.txt --greedyOptimizer=true --greedyCycles=3 --pdbOutputDir ./$fastaTmId --rulesFile /data03/CATM/files/TMRULES_byHbond.txt  --clusterSolutions true --rmsdCutoff=2.0 --startResNum=$catmStartNum --endResNum=$catmEndNum --threadStart=$threadStart --threadEnd=$threadEnd --tmStart=$tmStartNum --tmEnd=$tmEndNum --printTermEnergies true --printAxes true\"\; sleep 0.2\n");
			push (@mkdirSubmit, "mkdir ./$fastaTmId\n");
			my $tmseq = ">$fastaTmId\n";
			$tmseq = $tmseq . substr($sequence, $catmStartNum, ($catmEndNum-$catmStartNum)) . "\n";
			push (@tmSequences, $tmseq);
		}
	
	}
}
open(SUBMITFILE, ">catmSubmit.sh") or die "Did not make submit file \n";
foreach my $submitline (@catmSubmit) {
	print SUBMITFILE $submitline;
}
close (SUBMITFILE);

open(MKDIRFILE, ">catmDirs.sh") or die "Did not make mkdir file \n";
foreach my $mkdirLine (@mkdirSubmit) {
	print MKDIRFILE $mkdirLine;
}
close (MKDIRFILE);

open(TMSEQFILE, ">$fastaTM") or die "Did not make Fasta TM File\n!";
foreach my $fastaTMline (@tmSequences) {
	print TMSEQFILE $fastaTMline;
}
close(TMSEQFILE);

=pod
	my($currentId) = "NULL";
	my($uniprotName) = "NULL";
	my($detail);
	my($catmStartNum);
	my($catmEndNum);
	my($ext) = 3; #default to 3 residues on either side of the TM
	my($startExt) = $ext;
	my($endExt) = $ext;
	my($sequence) = "";

	for(my($i) = 0; $i < @data; $i++) {
		chomp($data[$i]);
		if(substr($data[$i],0,2) eq "ID") {
			my(@toks) = split(/\s+/,$data[$i]); 
			$uniprotName = $toks[1];
		}
		if($tmStartNum == "" || $tmEndNum == "") {	#Only look for TM region if one not supplied
			if(substr($data[$i],0,13) eq "FT   TRANSMEM") {
			### Assume we have only one TM region
				my(@toks) = split(/\s+/,$data[$i]); 
				$tmStartNum = $toks[2];
				$tmEndNum = $toks[3]; 
			}
		}
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
			my($len) = $tmEndNum - $tmStartNum + 1;
			while($tmStartNum <= $startExt) {
				$startExt--;
			}
			while(length($sequence) - $tmEndNum < $endExt) {
				$endExt--;
			}
			$catmStartNum = $tmStartNum - $startExt;
			$catmEndNum = $catmStartNum + $len + $endExt + $startExt - 1;
			my($tmSequence) = substr($sequence,$tmStartNum-1-$startExt,$len + $startExt + $endExt);
	
		}
	}
	
	
	my($threadStart) = 35 - (($catmEndNum - $catmStartNum + 1) - 3 - $startExt);
	my($threadEnd) = 33 - $endExt;
	
	open (SUBMITFILE, ">>./customTM.submit");
	print SUBMITFILE "pbsq --qsub \"/data03/CATM/bin/CATM_v20 --backboneCrd /data03/CATM/files/69-gly-residue-helix.crd --logFile ./$uniprotId/$uniprotId.log --helixGeoFile /data03/CATM/files/CENTROIDS_0.5_cutoff-4_byHbond --fullSequence $sequence --topFile /data03/CATM/topparsolvhbond/top_all22_prot.inp --parFile /data03/CATM/topparsolvhbond/par_all22_prot.inp --solvFile /data03/CATM/topparsolvhbond/solvpar22.inp --hBondFile /data03/CATM/topparsolvhbond/par_hbond_CA.inp --rotLibFile /data03/CATM/files/EBL_11-2011_CHARMM22.txt --greedyOptimizer=true --greedyCycles=3 --pdbOutputDir ./$uniprotId --rulesFile /data03/CATM/files/TMRULES_byHbond.txt --uniprotName $uniprotName --uniprotAccession $uniprotId --clusterSolutions true --rmsdCutoff=2.0 --startResNum=$catmStartNum --endResNum=$catmEndNum --threadStart=$threadStart --threadEnd=$threadEnd --tmStart=$tmStartNum --tmEnd=$tmEndNum --printTermEnergies true --printAxes true\"\; sleep 0.2\n";

}
=cut

#! /usr/bin/perl
#This program takes a UNIPROT ID goes online gets the information for that protein and prints out 
#a CATM submit file for that protein.
#
#This version accepts custom TM start and end sites.
#NOTE: Much of the information in this program is hard coded, and works off of the parameters 
#specified in the CATM HUMAN_UNIPROT_20131120 run.
#
#if(@ARGV != 1 || 3) {
#	print "Usage: getCATMSubmitFile_customTM.pl <protein_TM_list>\n";
#	print "List format: <UniprotId> <tmStartNum> <tmEndNum>\n";
#	exit;
#}
use LWP::Simple;

my $queryFile = shift;

open(QUERY, "$queryFile") or die "Could not open the list!\n";
LINE: while (my $protein = <QUERY>) {
	my ($uniprotId, $tmStartNum, $tmEndNum) = split (/\s+/,$protein);
	getstore("http://www.uniprot.org/uniprot/$uniprotId.txt","/data00/tmp/$uniprotId.txt");
	open(FILE,"/data00/tmp/$uniprotId.txt");
	my @data = <FILE>;
	close(FILE);
	
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

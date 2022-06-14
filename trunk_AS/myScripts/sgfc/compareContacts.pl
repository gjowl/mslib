# This script will accept two lists of contacts and compare the two.
# It will calculate the percentage of matching pairs and plot
# the two maps on top of one another.  The points on the map
# corresponding to shared points will be hilighted.
# This script will be useful for ranking computational models
# by how well they match coevolving positions.


my $program = "/exports/home/scondon/mslib/trunk_AS/bin/identifyInteractingResidues";
my $inputDir = shift;
my $sele1 = "chain a";
my $sele2 = "chain c+d";
my $outputDir = "/tmp/";

my $expContactsFile = shift;
my %expContacts;

open EXPCONTACTS, $expContactsFile or die "Error: did not open experimental contacts file!\n";
while( my $line = <EXPCONTACTS>) {
	chomp $line;
	$expContacts{$line} = "Experimental Contact";
}

for (my $i = 1; $i <= 600; $i++) {
	my $dirIndex = sprintf "%04s", $i;
	my $pdb = "$inputDir/$dirIndex/MC_Final.pdb";
	if (!-e $pdb) {
		next;
	}
	my $command = "$program --pdb $pdb --sele1 $sele1 --sele2 $sele2 --outputDir $outputDir";
	#print "$command\n";
	my @stdout = qx($command);
	my %pdbContacts;
	foreach my $contact (@stdout) {
		chomp($contact);
		#print "$contact\n";	
		$pdbContacts{$contact} = "Structural Contact";
	}
	
	my %sharedContacts;
	my $sharedCounter = 0;
	foreach my $key (keys %expContacts) {
		if (exists $pdbContacts{$key}) {
			#print "FOUND	$key\n";
			$sharedCounter++;
		}
	}
	my $numExpContacts = keys (%expContacts);
	print "$pdb: $sharedCounter/$numExpContacts\n";
}

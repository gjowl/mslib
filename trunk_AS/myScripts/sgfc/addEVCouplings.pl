#!usr/bin/perl
#This script will add up evolutionary couplings from EV-Fold, or a
#subset of them from a list. This will be useful for specifically
#calculating hotspots for inter-chain interactions

my $couplingsFile = shift;
my $rankThreshold = shift;

open COUPLINGS, $couplingsFile or die "ERROR: DID NOT OPEN $couplingsFile!\n";

my %hotspots;
my %resIDs;

#Example line:
#50 1 M 214 D 0 -0.0387938 1 0
#Rank Resi1 Resn1 Resi2 Resn2 MI DI
while (my $line = <COUPLINGS>) {
	chomp $line;
	my ($rank, $resi1, $resn1, $resi2, $resn2, $mi, $di) = split (' ', $line);
	if ($rank <= $rankThreshold) {
		$hotspots{$resi1} += $di;
		$hotspots{$resi2} += $di;
	}
}

foreach my $keys (sort keys %hotspots) {
	print "$keys	$hotspots{$keys}\n";
}


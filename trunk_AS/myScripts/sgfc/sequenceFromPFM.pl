#!usr/bin/perl

#This script will accept a position frequency matrix
#as input and output a sequence consistent with the
#PFM.
#


sub calc_probability {
	my $probRef = shift;
	my @probTable = @{$probRef};

	my $p = rand();
	my $cumulative = 0;
	for (my $i = 0; $i < @probTable; $i++) {
		$cumulative += $probTable[$i];
		if ($p <= $cumulative) {
			return $i;
		}
	}
}


my $pfmFile = shift;	#TSV With sequence labels as first line
open PFM, $pfmFile or die "ERROR: DID NOT OPEN PFM FILE!\n";

my $charLine = <PFM>;
chomp($charLine);

my @chars = split(/\t/, $charLine);
my @pfmTable;
while (my $line = <PFM>) {
	chomp($line);
	my @probs = split(/\t/, $line);

	#Make sure probabilities are normalized
	my @normProbs;
	my $probTotal = 0;
	foreach my $prob (@probs) {
		$probTotal += $prob;
	}
	foreach my $prob (@probs) {
		$normalized = $prob/$probTotal;
		push (@normProbs, $normalized);
	}

	push (@pfmTable, \@normalized);
}

my $sequence = '';

for (my $i = 0; $i < @pfmTable; $i++) {
	my $index = &calc_probability($pfmTable[$i]);
	my $char = $chars[$index];
	$sequence .= $char;
}
print "$sequence\n";


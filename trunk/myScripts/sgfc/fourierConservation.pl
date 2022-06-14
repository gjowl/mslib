# This program will perform Fourier power spectrum analysis on 
# a protein's conservation values at each position to determine
# if there is a helical pattern of conservation, based on the article
# 
#   J.L. Cornette, K.B. Cease, H. Margalit, J.L. Spouge, J.A. Berzofsky, C. DeLisi
#   "Hydrophobicity scales and computational techniques for detecting amphipathic 
#   structures in proteins"
#   J. Mol. Biol., 195 (3) (1987), pp. 659–685, DOI: 10.1016/0022-2836(87)90189-6
# 
# The equation for the fourier transform of a sequence of length n is
# 
#           _                 _ 2     _                 _ 2
#          |  n                |     |  n                |
#   p(w) = |  Sum (Cn*cos(nw)) |  +  |  Sum (Cn*sin(nw)) |  
#          |_ i=1             _|     |_ i=1             _|
# 
# 
# Where	w = helical angle (degrees)
#       Cn = Positional conservation score for position n
# 
# A peak around w ~= 100 indicates a pattern of helical periodicity


#!usr/bin/perl
use Math::Trig;


my @periodicity;
my @sequence;
my $conservationFile = shift;

open CONSERVATION, "$conservationFile" or die "did not open conservation file\n";
while (my $line = <CONSERVATION>) {
	chomp ($line);
	#my ($pos, $res, $Cn, $gaps) = split (/\s+/,$line);
	my ($pos, $res, $Cn, $gaps);
	$Cn =$line;
	if ($gaps =~ /\*/) {	# Columns with too many gaps are not scored
				# and are marked with a '*'
		next;
	} else {
		print "$Cn\n";
		push (@sequence, $Cn);
	}
}
close CONSERVATION;
for (my $w = 0; $w < 180; $w++) {
	my $rad = deg2rad($w);
	my $cos;
	my $sin;
	for (my $res = 0; $res < @sequence; $res++) {
		$cos = $cos + ($sequence[$res] * cos($res * $rad));
		$sin = $sin + ($sequence[$res] * sin($res * $rad));
	}

	my $cosSq = $cos * $cos;
	my $sinSq = $sin * $sin;
	my $p = $cosSq + $sinSq;
	push (@periodicity, $p);
}

for (my $w = 0; $w < @periodicity; $w++) {
	print "$w,$periodicity[$w]\n";
}



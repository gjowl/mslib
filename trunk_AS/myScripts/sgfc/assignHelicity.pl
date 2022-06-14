#!usr/bin/perl
#This script will accept a list of phi, psi, and o-n4 distances
#and determine whether or not these values correspond to a helical
#residue in a protein chain.

#Classification scheme: if all 3 values are within threshold,
#call it a helix. If one or two of the values are outside the
#threshold, call it almost a helix. If all of the values are
#outside the threshold, call it a non-helix.
use POSIX;

my $ramachandranBinFile = '/data05/scondon/Projects/FtsB-FtsL/data/ramachandranBins_helix.tsv';

open BIN, $ramachandranBinFile or die "ERROR: DID NOT OPEN BIN FILE!\n";
my $binHeader = <BIN>;

my %favorableBins;
my %allowedBins;

while (my $line = <BIN>) {
	chomp $line;
	my ($phi, $psi, $class) = split (/\t/, $line);
	
	#Do something with the phi and psi bins
	if ($class eq 'Favorable') {
		push (@{ $favorableBins{$phi} }, $psi);
	} elsif ($class eq 'Allowed') {
		push (@{ $allowedBins{$phi} }, $psi);
	}
}

#foreach my $keys (sort keys %favorableBins) {
#	print "$keys\n";
#	foreach my $psi (@{$favorableBins{$keys}}) {
#		print "$keys $psi\n";
#	}
#}

my $on4min = 2;
my $on4mid = 4.2;
my $on4max = 5;

my $phitable = shift;
my $psitable = shift;
my $on4table = shift;

my @phiframes;
my @psiframes;
my @on4frames;

my @phiheader;
my @psiheader;
my @on4header;

open PHI, $phitable or die "ERROR: DID NOT OPEN $phitable!\n";

my $philine = <PHI>;
chomp $philine;
@phiheader = split (/\t/, $philine);

open PSI, $psitable or die "ERROR: DID NOT OPEN $psitable!\n";
my $psiline = <PSI>;
chomp $psiline;
@psiheader = split (/\t/, $psiline);

open ON4, $on4table or die "ERROR: DID NOT OPEN $on4table!\n";
my $on4line = <ON4>;
chomp $on4line;
@on4header = split ( /\t/, $on4line);


#Check to make sure tables have same position and chain Ids
for (my $i = 0; $i < @phiheader; $i++) {
	if ($phiheader[$i] ne $psiheader[$i] or $phiheader[$i] ne $on4header[$i]) {
		die "ERROR: PHI, PSI, or O-N4 headers do not match!\n";
	}
}
shift @phiheader;
shift @psiheader;
shift @on4header;

my @positionHelicity;
for (my $i = 0; $i < @phiheader; $i++) {
	my %helicity;
	$helicity{'pos'} = $phiheader[$i];
	#@helicity{'timeline'};
	$helicity{2} = 0;
	$helicity{1} = 0;
	$helicity{0} = 0;
	push (@positionHelicity, \%helicity);
}

#Put frame values into a 2d array for each item
while (my $line = <PHI>) {
	chomp $line;
	my @vals = split (/\t/, $line);
	shift @vals;	#get rid of header
	for (my $i = 0; $i < @vals; $i++) {
		push (@{$phiframes[$i]}, $vals[$i]);
	}
}


while (my $line = <PSI>) {
	chomp $line;
	my @vals = split (/\t/, $line);
	shift @vals;	#get rid of header
	for (my $i = 0; $i < @vals; $i++) {
		push (@{$psiframes[$i]}, $vals[$i]);
	}
}


while (my $line = <ON4>) {
	chomp $line;
	my @vals = split (/\t/, $line);
	shift @vals;	#get rid of header
	for (my $i = 0; $i < @vals; $i++) {
		push (@{$on4frames[$i]}, $vals[$i]);
	}
}


#Loop over 2d arrays and compute helicity
my @helicityFrames;
for (my $i = 0; $i < @phiframes; $i++) {
	
	for (my $j = 0; $j < @{$phiframes[$i]}; $j++) {
		
		my $phi_ij = $phiframes[$i][$j];
		my $psi_ij = $psiframes[$i][$j];
		my $on4_ij = $on4frames[$i][$j];


		my $helical = &helicalThreshold($phi_ij, $psi_ij, $on4_ij);
		push (@{ $positionHelicity[$i]->{'timeline'} }, $helical);
		$positionHelicity[$i]->{$helical} += 1;

	}
}


#Print out results
for (my $i = 0; $i < @positionHelicity; $i++) {
	my $posInfo = $positionHelicity[$i];
	print "$posInfo->{'pos'}	$posInfo->{2}	$posInfo->{1}	$posInfo->{0}\n";
}

for (my $i = 0; $i < @positionHelicity; $i++) {
	print "$positionHelicity[$i]->{'pos'}\t";
}
print "\n";
for (my $frame = 0; $frame < @{$phiframes[0]}; $frame++) {
	for (my $i = 0; $i < @positionHelicity; $i++) {
		my $posInfo = $positionHelicity[$i];
		print "$posInfo->{'timeline'}[$frame]\t";
	}
	print "\n";
}









sub helicalThreshold {
	my $phival = shift;
	my $psival = shift;
	my $on4val = shift;

	#Binary values to see if phi psi and on4 are within thresholds
	my $ramaClass = 0;
	my $distClass = 0;

	my $phi_bin = &roundDown($phival, 10);
	my $psi_bin = &roundDown($psival, 10);

	foreach my $favPsi (@{ $favorableBins{$phi_bin} }) {
		if ($favPsi == $psi_bin) {
			$ramaClass = 2;
		}
	}
	foreach my $allowPsi (@{ $allowedBins{$phi_bin} }) {
		if ($allowPsi == $psi_bin) {
			$ramaClass = 1;
		}
	} 
	


	if ($on4val >= $on4min and $on4val <= $on4mid) {
		$distClass = 2;
	} elsif ($on4val >= $on4mid and $on4val <= $on4max) {
		$distClass = 1;
	} else {
		$distClass = 0;
	}

	#CLASS ASSIGNMENT BASED ON RAMA AND DIST:
	#X is not helical
	#h is almost helical
	#H is fully helical
	#
	#             R
	#	   X  h  H  
	#	    0  1  2
	#	X 0 x  x  x	
	#D	h 1 x  h  h
	#	H 2 x  h  H
	#
	my @classes = {$ramaClass, $distClass};
	my $helicalVal = &min($ramaClass, $distClass);

	return $helicalVal;
}

sub roundDown {
	my $number = shift;
	my $roundTo = shift || 1;

	return int(floor($number/$roundTo))*$roundTo;
}

sub min {
    my ($min, @vars) = @_;
    for (@vars) {
        $min = $_ if $_ < $min;
    }
    return $min;
}

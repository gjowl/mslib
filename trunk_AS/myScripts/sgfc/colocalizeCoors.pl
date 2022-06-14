#!usr/bin/perl
#This script will accept two lists of X-Y coordinates and identify which sets of coordinates overlap using a greedy algorithm.
#The centroid of each set of coors will be used to align the two sets and starting from a random spot in one list, the spot
#in the other list with the lowest distance will be selected as the match. The set of matches will then be printed out.
#


my $coorfile_tar = shift;
my $coorfile_src = shift;

open COOR_tar, $coorfile_tar or die "ERROR: DID NOT OPEN COORDINATE FILE _tar!\n";
open COOR_src, $coorfile_src or die "ERROR: DID NOT OPEN COORDINATE FILE _src!\n";

my @coors_tar;
my @coors_src;

my $centroid_tarX;
my $centroid_tarY;
my $centroid_srcX;
my $centroid_srcY;


while (my $line = <COOR_tar>) {
	chomp $line;
	my ($coorId, $coorX, $coorY) = split (/\s+/, $line);
	my %coorInfo;
	$coorInfo{'Id'} = $coorId;
	$coorInfo{'X'} = $coorX;
	$coorInfo{'Y'} = $coorY;
	$centroid_tarX += $coorX;
	$centroid_tarY += $coorY;

	push (@coors_tar, \%coorInfo);
}


while (my $line = <COOR_src>) {
	chomp $line;
	my ($coorId, $coorX, $coorY) = split (/\s+/, $line);
	my %coorInfo;
	$coorInfo{'Id'} = $coorId;
	$coorInfo{'X'} = $coorX;
	$coorInfo{'Y'} = $coorY;
	$centroid_srcX += $coorX;
	$centroid_srcY += $coorY;

	push (@coors_src, \%coorInfo);
}

$centroid_tarX /= @coors_tar;
$centroid_tarY /= @coors_tar;

$centroid_srcX /= @coors_src;
$centroid_srcY /= @coors_src;


#Shift the centroid of coordinates _src to the centroid of coordinates _tar to align
foreach my $coor (@coors_src) {
	$coor->{'X'} -= ($centroid_srcX-$centroid_tarX);
	$coor->{'Y'} -= ($centroid_srcY-$centroid_tarY);
}

#Identify closest match between spots in coors_tar and coors_src
print "COOR_tar	COOR_src	DIST\n";
foreach my $coor (@coors_tar) {
	my $bestId = 'asdf';
	my $bestDist = 99999999999999999999;
	#print "=================================================\n";	
	foreach my $coor_src (@coors_src) {
		my $dist = sqrt(($coor->{'X'} - $coor_src->{'X'})**_src + ($coor->{'Y'} - $coor_src->{'Y'})**_src);
		#print "$coor->{'X'} $coor->{'Y'} $coor_src->{'X'} $coor_src->{'Y'} $dist\n";
		if ($dist < $bestDist) {
			$bestId = $coor_src->{'Id'};
			$bestDist = $dist;
		}
	}

	print "$coor->{'Id'}	$bestId	$bestDist\n";
}

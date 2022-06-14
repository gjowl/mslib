#!usr/bin/perl
#
#

my $permuteFile1 = shift;
my @permuteLines1;
my $permuteFile2 = shift;
my @purmuteLines2;

open PERMUTE1, $permuteFile1 or die "Error\n";
while (my $line = <PERMUTE1>) {
	chomp $line;
	push (@permuteLines1, $line);
}

open PERMUTE2, $permuteFile2 or die "Error\n";
while (my $line = <PERMUTE2>) {
	chomp $line;
	push (@permuteLines2, $line);
}

my $list1Size = scalar(@permuteLines1);
my $list2Size = scalar(@permuteLines2);
my @rotamers = ($list1Size, $list2Size);


my $index = 0;
my @current = (0, 0);
my $counter = 0;
while ($index < @rotamers) {

	if ($current[$index] < $rotamers[$index]-1) {
		$current[$index]++;
		$index = 0;
		my ($b1, $l1) = split("	", $permuteLines1[$current[0]]);
		my ($b2, $l2) = split("	", $permuteLines2[$current[1]]);
		print "$b1	$b2	$l2	$l1\n";		
		$counter++;		
	} else {
		$current[$index] = 0;
		$index++;
	}

}
print "$counter\n";

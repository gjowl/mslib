#!/usr/bin/perl

$program_name    = "mini.pl";
$program_author  = "Alessandro Senes";
$program_version = "1.4.0"; 
$program_date    = "13 November 2007";

use lib "/library/perllibs/";
use Molecule;
use Atomgeometry;
use Getopt::Long;

&GetOptions(
	"crd=s" => \$crdFile,
	"nPatch=s@" => \@nPatch,
	"cPatch=s@" => \@cPatch,
	"cycles=s" => \$cycles,
	"residue=s@" => \@residue,
	"chain=s@" => \@chain,
	"steps=s" => \$steps,
	"atom=s@" => \@atoms, 
	"fixedAtomTypes=s@" => \@fixedAtomTypes,
	"freeChains=s@" => \@freeChains,
	"fixedChains=s@" => \@fixedChains,
	"topfile=s" => \$chrm_toph,
	"parfile=s" => \$chrm_param,
	"noConstrain" => \$no_constrain_flag,
	"constrForce=s" => \$constrForce,
	"charmm_exec=s" => \$charmm_exec,
	"verbose!" => \$verbose_flag,
	"title=s" => \$crdTitle,
	"waterChains=s@" => \@waterChains,
	"dielectric=s" => \$diel,
	);

if (!defined $charmm_exec) {
	$charmm_exec = "charmm";
}

if (!defined $crdTitle) {
	$crdTitle = "Minimized crd";
}

if (!defined $diel) {
	$diel = "80.0";
}

if (!defined $crdFile) {
	$crdFile = shift;
	if (!defined $crdFile) {

		print "Program $program_name v.$program_version ($program_date)\n";
		print "\n";
		print "Minimal usage:\n";
		print "mini.pl --crd <CRD>\n";
		print "\n";
		print "To specify the chain termini (add one --nPatch --cPatch per chain):\n";
		print "mini.pl --crd <CRD> --nPatch <NTER|GLYP|PROP|ACE> --cPatch <CTER|CT1|CT2>\n";
		print "\n";
		print "To use constraint minimization (default is ON , default is 3 cycles)\n";
		print "mini.pl --cord <CRD> --constrForce 100 [--cycles <N>]\n";
		print "\n";
		print "To use constraint minimization (suggested 100, default is 3 cycles)\n";
		print "mini.pl --cord <CRD> --constrForce 100 [--cycles <N>]\n";
		print "\n";
		print "To minimize only certain residues\n";
		print "mini.pl --chain <chainid> --residue <resnum> --chain <chainid> --residue <resnum>\n";
		print "\n";
		print "To minimize only some atoms on one or more residues (multiple atoms must be in quotes)\n";
		print "mini.pl --chain <chainid> --residue <resnum> --atom <NAME> --atom \"<NAME> <NAME> <NAME>...\"\n";
		print "\n";
		print "To held fix an atom type (say the CA or the whole backbone)\n";
		print "mini.pl --fixedAtomTypes <NAME>  --fixedAtomTypes <NAME> --fixedAtomTypes <NAME>\n";
		print "\n";
		print "To minimise only some chains \n";
		print "mini.pl --freeChains <NAME>  --freeChains <NAME> --freeChains <NAME>\n";
		print "\n";
		print "To hold some chains fixed \n";
		print "mini.pl --fixedChains <NAME>  --fixedChains <NAME> --fixedChains <NAME>\n";
		print "\n";
		print "Option --title \"some title\" let you set the title of the minimized crd\n";
		exit;
	}
}
$constrained_flag = 1;
if (defined $no_constrain_flag || (defined $constrForce && $constrForce == 0)) {
	$constrained_flag = 0;
} else {
	if (!defined $constrForce) {
		$constrForce = 100;
	}
	$constrForceHydrogen = $constrForce / 10;
}
	
if (!defined $cycles) {
	if ($constrained_flag) {
		$cycles = 3;
	} else {
		$cycles = 1;
	}
}

if (!defined $steps) {
	$steps = 50;
}

my @atomList;
if (scalar(@chain) > 0) {
	if (scalar(@residue) != scalar(@chain)) {
		die "Number of --chain and --residue options should be identical\n";
	}
	if (scalar(@atom) > 0 && scalar(@atoms) != scalar(@chain)) {
		die "Number of --chain and --atom options should be identical\n";
	}

	for $i (0..$#chain) {
		my @individualAtoms = split(" ",$atoms[$i]);
		for $j (0..$#individualAtoms) {
			push @atomList, {chainid => $chain[$i], resid => $residue[$i], atomname => $individualAtoms[$j]};
		}
	}
}

if (!defined $chrm_toph) {
	$chrm_toph  = "/library/charmmTopPar/top_all22_prot.inp";
}
if (!defined $chrm_param) {
	$chrm_param = "/library/charmmTopPar/par_all22_prot.inp";
}

	
my $tmpName  = $crdFile;
$tmpName =~ s/\.crd$//;;
$charmmScriptFile = $tmpName . "-mini";
$crdOutputFile = $tmpName . "-mini.crd";
if (-e $crdOutputFile) {
	system "mv $crdOutputFile $crdOutputFile.bak";
	print "moved $crdOutputFile to $crdOutputFile.bak\n";
}
$crdOutputFileLowerCase = $crdOutputFile;
$crdOutputFileLowerCase =~ tr/A-Z/a-z/;


###############################################
# Open the crd file
###############################################
$crd = new Molecule;
$crd->ReadCrd($crdFile);
#printf "$crdname has %u atoms\n", $#$crd + 1;

my $prevSegid;
my $prevResnum;
my @fix;
my %seq;
#my %generate;
my @generate;
for $i (0..$#$crd) {

	if ($i == 0 || $$crd[$i]->segid ne $prevSegid || $$crd[$i]->resnum ne $prevResnum) {
		if ($i == 0 || $$crd[$i]->segid ne $prevSegid) {
			push @generate, $$crd[$i]->segid;
		#	if (scalar(@nPatch) < scalar(@generate)) {
		#		if ($$crd[$i]->resid eq "GLY") {
		#			push @nPatch, "GLYP";
		#		} elsif ($$crd[$i]->resid eq "PRO") {
		#			push @nPatch, "PROP";
		#		} else {
		#			push @nPatch, "NTER";
		#		}
		#	}
		#	if (scalar(@CPatch) < scalar(@generate)) {
		#		push @cPatch, "CTER";
		#	}
		}
		push @{$seq{$$crd[$i]->segid}}, $$crd[$i]->resid;
		$prevSegid  = $$crd[$i]->segid;
		$prevResnum = $$crd[$i]->resnum;
		if ($$crd[$i]->resnum ne $#{$seq{$$crd[$i]->segid}} + 1) {
			if ($#{$seq{$$crd[$i]->segid}} < 999) {
				push @{$fix[0]}, sprintf "renam resid x%s sele atom %s %d %s end\n", $#{$seq{$$crd[$i]->segid}} + 1, $$crd[$i]->segid, $#{$seq{$$crd[$i]->segid}} + 1, $$crd[$i]->atomid;
				push @{$fix[1]}, sprintf "renam resid %s sele segid %s .and. resid x%s end\n", $$crd[$i]->resnum, $$crd[$i]->segid, $#{$seq{$$crd[$i]->segid}} + 1;
			} else {
				push @{$fix[0]}, sprintf "renam resid y%s sele atom %s %d %s end\n", $#{$seq{$$crd[$i]->segid}} - 999, $$crd[$i]->segid, $#{$seq{$$crd[$i]->segid}} + 1, $$crd[$i]->atomid;
				push @{$fix[1]}, sprintf "renam resid %s sele segid %s .and. resid y%s end\n", $$crd[$i]->resnum, $$crd[$i]->segid, $#{$seq{$$crd[$i]->segid}} - 999;
			}
		}
	}
	#$generate{$$crd[$i]->segid} = 1;

}


###############################################
# START BUILDING THE CHARMM INPUT FILE
# 
# 
###############################################

###############################################
# Make a title and wrap it to 80 characters 
# ("* " + 78)
###############################################
$time = localtime();
$title = sprintf "Charmm script generated by %s v.%s at %s on %s", $0, $program_version, $time, $ENV{HOST};
$skip = 78;
undef @titleLines;
while (length($title)>78) {
	if (substr($title, $skip, 1) eq " ") {
		push @titleLines, substr($title, 0, $skip);##
		$title = substr($title, $skip + 1, length($title) - $skip - 1);
	}
	$skip--;
	if ($skip == -1) {
		push @titleLines, substr($title, 0, 78);##
		$title = substr($title, 78, length($title) - 78);
	}
}
if (length($title)>0) {
	push @titleLines, $title;
}

for $i (0..$#titleLines) {			
	push @charmmScript, "* $titleLines[$i]\n";
}
push @charmmScript, "*\n";

push @charmmScript, <<"END";

prnl 5

open read unit 1 card name "$chrm_toph"
read rtf card unit 1
close unit 1

open read unit 1 card name "$chrm_param"
read para card unit 1
close unit 1

END

push @charmmScript, <<"END";

! Generate topology of each segment
! =================================

END

my(%waters);
foreach $i (@waterChains) {
	$waters{$i} = 1;
}

#for $segment (sort keys %generate) {
for $i (0..$#generate) {
	push @charmmScript, "read sequence card\n";
	push @charmmScript,  "*\n";
	push @charmmScript, sprintf "%u\n", $#{$seq{$generate[$i]}} + 1;
	undef $tmp;
	for $j (0..$#{$seq{$generate[$i]}}) {
		if (length($tmp) + length($seq{$generate[$i]}[$j]) + 1 > 60) {
			$tmp =~ s/ *$//;
			push @charmmScript, "$tmp\n";
			undef $tmp;
		}	
		$tmp .= "$seq{$generate[$i]}[$j] ";

	}
	$tmp =~ s/ *$//;
	push @charmmScript, "$tmp\n";
	if ($#nPatch >= $i && $#cPatch >= $i) {
		push @charmmScript,  "gene $generate[$i] first $nPatch[$i] last $cPatch[$i] setup\n";
	} else {
		if(exists $waters{$generate[$i]} ) {
			push @charmmScript,  "gene $generate[$i] setup noangle nodihe\n";
		} else {
			push @charmmScript,  "gene $generate[$i] setup\n";
		}
	}
	push @charmmScript,  "\n";
}


push @charmmScript, @{$fix[0]};
push @charmmScript, @{$fix[1]};


push @charmmScript,  <<"END";

! ============================================

open read unit 1 card name "$crdFile"
read coor card resid unit 1
close unit 1

! List any missing non-hydrogen atoms
! ===================================

print coor sele .not. (hydrogen .or. init) end

! Copy coordinates to comparison set
! and save the internal coordinates
!================================
coor copy comp
ic fill preserve
!ic save

ic param 
ic build

END

#if (defined $residue && defined $chain) {
#	push @charmmScript,  <<"END";
#
#define backb select type N .or. type CA -
#.or. type C .or. type O .or. type H -
#.or. type CY .or. type OY .or. type CAY -
#.or. type HN .or. type HY1 -
#.or. type HY2 .or. type HY3 .or. type NT -
#.or. type HT1 .or. type HT2 .or. type HNT -
#.or. type CAT .or. type OTX .or. type OT1 -
#.or. type OT2 .or. type NXT -
#end
#
#defi fixpro select backb .or. .not. ( -
#.byres. atom $chain $residue N) end
#
#
#print coor sele .not. (fixpro) end
#
#! constrain the fix atoms so that no
#! energy terms will be calculated among
#! fixed atoms (this is important later
#! when we use gete instead of inte)
#!================================
#cons fix sele fixpro end
#
#END
#}

if (scalar(@freeChains) > 0) {
	push @charmmScript,  "defi fixchainspro select .not. ( -\n";

		for $i (0..$#freeChains) {
			if ($i == 0) {
				push @charmmScript,  "segid $freeChains[$i] -\n";
			} else {
				push @charmmScript,  ".or. segid $freeChains[$i] -\n";
			}
		}
		
		push @charmmScript,  ") end\n";

		push @charmmScript,  <<"END";

print coor sele .not. (fixchainspro) end

! constrain the fix chains so that no
! energy terms will be calculated among
! fixed chains (this is important later
! when we use gete instead of inte)
!================================
cons fix sele fixchainspro end

END


}

if (scalar(@fixedChains) > 0) {
	push @charmmScript,  "defi fixchainspro select ( -\n";

		for $i (0..$#freeChains) {
			if ($i == 0) {
				push @charmmScript,  "segid $freeChains[$i] -\n";
			} else {
				push @charmmScript,  ".or. segid $freeChains[$i] -\n";
			}
		}
		
		push @charmmScript,  ") end\n";

		push @charmmScript,  <<"END";

print coor sele .not. (fixchainspro) end

! constrain the fix chains so that no
! energy terms will be calculated among
! fixed chains (this is important later
! when we use gete instead of inte)
!================================
cons fix sele fixchainspro end

END


}
if (scalar(@chain) > 0) {
	if (@atomList) {
		# do only selected atoms of selected residues 
		push @charmmScript,  "defi fixpro select .not. ( -\n";

		for $i (0..$#atomList) {
			if ($i == 0) {
				push @charmmScript,  "atom $atomList[$i]{chainid} $atomList[$i]{resid} $atomList[$i]{atomname} -\n";
			} else {
				push @charmmScript,  ".or. atom $atomList[$i]{chainid} $atomList[$i]{resid} $atomList[$i]{atomname} -\n";
			}
		}
		
		push @charmmScript,  ") end\n";

		push @charmmScript,  <<"END";

print coor sele .not. (fixpro) end

! constrain the fix atoms so that no
! energy terms will be calculated among
! fixed atoms (this is important later
! when we use gete instead of inte)
!================================
cons fix sele fixpro end

END
	} else {
		# do whole residues including backbone
		push @charmmScript,  "defi fixpro select .not. ( -\n";
		for $i (0..$#chain) {
			if ($i == 0) {
				push @charmmScript,  ".byres. atom $chain[$i] $residue[$i] N -\n";
			} else {
				push @charmmScript,  ".or. .byres. atom $chain[$i] $residue[$i] N -\n";
			}
		}
		push @charmmScript,  ") end\n";

		push @charmmScript,  <<"END";
print coor sele .not. (fixpro) end

! constrain the fix atoms so that no
! energy terms will be calculated among
! fixed atoms (this is important later
! when we use gete instead of inte)
!================================
cons fix sele fixpro end

END
	}

} elsif (@fixedAtomTypes) {
	push @charmmScript,  "defi fixpro select ( -\n";

	for $i (0..$#fixedAtomTypes) {
		if ($i == 0) {
			push @charmmScript,  "type $fixedAtomTypes[$i] -\n";
		} else {
			push @charmmScript,  ".or. type $fixedAtomTypes[$i] -\n";
		}
	}
	
	push @charmmScript,  ") end\n";

	push @charmmScript,  <<"END";

print coor sele .not. (fixpro) end

! constrain the fix atoms so that no
! energy terms will be calculated among
! fixed atoms (this is important later
! when we use gete instead of inte)
!================================
cons fix sele fixpro end

END
}

push @charmmScript,  <<"END";

! to do three rounds of minimization
!================================
set 2 1
label min_loop
END

if ($constrained_flag) {
	push @charmmScript,  <<"END";
cons harm force $constrForce sele .not. hydrogen end
cons harm force $constrForceHydrogen sele hydrogen end
END
}

push @charmmScript,  <<"END";

update ihbfrq 0 inbfrq 10 -
atom shift vdw vshift -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 -
wmin 1.5 rdiel eps $diel e14fac 1.0 nbxmod 5

mini abnr nstep $steps ihbfrq 0 inbfrq 0

incre 2 by 1
if 2 le $cycles goto min_loop


open write unit 2 card name "$crdOutputFile"
write coor card unit 2
* $crdTitle
*
close unit 2

END

open (FILE, ">$charmmScriptFile\.inp") || die "Cannot open $charmmScriptFile.inp for writing\n";
print FILE @charmmScript;
close FILE;

if ($verbose_flag) {
	print "\n";
	print "######################################################################################\n";
	print " ABOUT TO RUN CHARMM MINIMIZATION WITH THE FOLLOWING PARAMETERS\n";
	print "\n";
	print " - Input $crdFile, output $crdOutputFile\n";
	print " - Charmm executable $charmm_exec, topology file $chrm_toph, parameter file $chrm_param\n";
	if ($constrained_flag) {
		print " - Running $cycles cycle(s) of $steps steps of constrained adopted basis Newton-Raphson minimization, with force constants: non-hydrogen atoms: $constrForce, hydrogen atoms: $constrForceHydrogen\n";
	} else {
		print " - Running $cycles cycle(s) of $steps steps of adopted basis Newton-Raphson minimization\n";
	}
	print " - Chains' terminal patches:";
	for (my $i=0; $i<@nPatch; $i++) {
		print " $nPatch[$i]/$cPatch[$i]"
	}
	print "\n";
	if (scalar(@atoms) > 0) {
		print " - Minimize only residue $chain $residue, atoms";
		for (my $i=0; $i<@atoms; $i++) {
			print " $atoms[$i]";
		}
		print "\n";
	}
	if (scalar(@fixedAtomTypes) > 0) {
		print " - The following atoms types will be kept fixed:";
		for (my $i=0; $i<@fixedAtomTypes; $i++) {
			print " $fixedAtomTypes[$i]";
		}
		print "\n";
	}
	print " - Check $charmmScriptFile.inp for the charmm input file, and $charmmScriptFile.out for the output\n";
	print "\n";
	print "Program $program_name v.$program_version ($program_date)\n";
	print "######################################################################################\n";
	print "\n";
}

system "$charmm_exec < $charmmScriptFile.inp > $charmmScriptFile.out";
if ($crdOutputFile ne $crdOutputFileLowerCase) {
	if (-e $crdOutputFileLowerCase) {
		system "mv $crdOutputFileLowerCase $crdOutputFile";
		print "Minimization run succesfully\n";
	} else {
		print "Minimization did not run succesfully\n";
	}
} else {
	if (-e $crdOutputFile) {
		print "Minimization run succesfully\n";
	} else {
		print "Minimization did not run succesfully\n";
	}
}

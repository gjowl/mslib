# Script that makes a bunch of submit files for condor
# Will read in a config script and make as many copies
# as specified in the argument.  It will copy the contents
# of the config and replace any instances of '~~' with
# a number corresponding to the job.
# usage make_submit <config> <numFiles>
#!usr/bin/perl
use File::Basename;
use Getopt::Long;
use strict;

my $executable;
my $config;
my $outputDir;
my $numFiles;
my $makeOneLiner;

GetOptions (	"executable=s"	=> \$executable,
		"config=s"	=> \$config,
		"outputDir=s"	=> \$outputDir,
		"numFiles=i"	=> \$numFiles,
		"makeOneLiner"	=> \$makeOneLiner)
or die ("perl make_submit.pl --executable <> --config <> --outputDir <> --numFiles <> \n");

		



my ($configName, $configDir, $configType) = fileparse($config, '\..*');
my $numDigits = length ($numFiles);

open (CONFIG,"<", $config) or die "Error: did not open $config!\n";

my $config_contents;
my %config_hash;
my @config_args;
while (my $line = <CONFIG>) {
	chomp $line;
	if ($line ne "") {
		my ($arg, $argval) = split(/\s+/, $line, 2);
		$config_hash{$arg} = $argval;
		push(@config_args, $arg);
	}
}

my $bigsub = $outputDir . "/bigSubmit.sh";
open BIGSUBMIT, ">$bigsub" or die "ERROR: did not write a batch submit file!\n";
for (my $i = 1; $i <= $numFiles; $i++) {
	my $digitFormat = "%0$numDigits" . "s";
	my $fileIndex = sprintf $digitFormat, $i;
	my $subDir = $outputDir . "/$fileIndex";
	if (!-d $subDir) {
		mkdir($subDir, 0755) or die "ERROR: DID NOT MAKE $subDir\n";
	}
	
	
	my $newFileName = $subDir . "/$configName" . ".config";
	#$newFileText =~ s/\Q~~\E/$fileIndex/g;
	#print "$newFileText\n";
	
	my @inputFiles;
	open (NEWSUB, ">", $newFileName) or die "ERROR: did not open new sub file $newFileName for writing!\n";
	foreach my $arg (@config_args) {
		my $argval = $config_hash{$arg};
		$argval =~ s/\Q~~\E/$fileIndex/g;
		if (-f $argval) {
			push(@inputFiles, $argval);
		}
		print NEWSUB "$arg $argval\n";
	}
	close NEWSUB;
	my $inputList = join(", ", @inputFiles);
	my $condorFile = $subDir . "/$configName" . ".condor";
	my $condorOut = $subDir . "/$configName" . ".condorout";
	my $condorErr = $subDir . "/$configName" . ".condorerr";
	my $condorLog = $subDir . "/$configName" . ".condorlog";
	(my $condorText = <<"	END_CONDORTXT") =~ s/^\s//gm;
	executable = $executable
	arguments = "--configfile $newFileName"
	should_transfer_files = Yes
	output = $condorOut
	error = $condorErr
	log = $condorLog
	getenv = True
	queue 1
	END_CONDORTXT
	
	open (CONDORSUB, ">", $condorFile) or die "Error: did not open new condor sub file $condorFile for writing!\n";
	print CONDORSUB $condorText;
	close CONDORSUB;

	print BIGSUBMIT "condor_submit $condorFile\n";
	print BIGSUBMIT "sleep 0.5\n";
}

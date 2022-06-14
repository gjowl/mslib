#!usr/bin/perl
use strict;




#Given undirected network and start node, performs a breadth-first search
#to identify all reachable nodes. Given a parameter for the maximum depth,
#the search will be limited to a path of that length.

#TODO Read graph file
#TODO put graph into adjacency list
#TODO make breadth-first search algorithm with max depth
#TODO print out reachable nodes


#breadthFirstSearch.pl <graphFile> <startNode> <maxDepth>

###############################################################################
#
#			MAIN
#
###############################################################################


my $graphFile = shift;
my $startNode = shift;
my $maxDepth = shift;



#Structure of graph file for the following graph:

#    n3
#     |
#    n1-n2
#        |
#       n4

#n1 n2
#n1 n3
#n2 n4

print "graphFile: $graphFile\n";
print "startNode: $startNode\n";
print "maxDepth: $maxDepth\n";
if (!$maxDepth) {
	$maxDepth = 99999999999999;
}

open GRAPHFILE, $graphFile or die "ERROR: $graphFile not opened!\n";
my @edges;
while (my $line = <GRAPHFILE>) {
	chomp $line;
	push (@edges, $line);
}
close GRAPHFILE;
my $graphRef = &read_graphFile(\@edges);
my %graph = %{$graphRef};
#&print_hashofhashes(\%graph);

my $nthNeighbors = &breadthFirstSearch($startNode, $maxDepth, \%graph);

foreach my $neighbors(keys %{$nthNeighbors}) {
	print "$neighbors\n";
}

#my %paths;
#my %targets;
#foreach my $target(@targetList) {
#	$targets{$target} = 'target';
#}
#foreach my $source (@sources) {
#	my @visited = ($source);
#	my $pathWeight = 0;
#	&depthFirstSearch($source, \@visited, $pathWeight, \%paths, \%graph, \%targets);
#}
#
#&print_paths (\%paths, $k, $pathsFile);
	
#print "$numNodes\n";
#foreach my $source (@sources) {
#	print "$source ";
#}
#print "\n";

#foreach my $target (@targets) {
#	print "$target ";
#}
#print "\n";
#&print_hashofhashes(\%graph);











###############################################################################
#
#			SUBROUTINES
#
###############################################################################



sub read_graphFile{
	
	my %graph;	#Adjacency list with following format:
	#{
	#	1: {2: 0.1, 3: 0.3},
	#	2: {1: 0.1, 4: 0.8},
	#	3: {1: 0.3, 4: 0.2},
	#	4: {2: 0.8, 3: 0.2}
	#}
	my @_edges = @{$_[0]};


	foreach my $edge (@_edges) {
		my ($n1, $n2, $w) = split(/\s+/, $edge);
		$graph{$n1}{$n2} = $w;
		$graph{$n2}{$n1} = $w;
	}
		

	
	return \%graph;
}

sub print_hashofhashes {
	my %hoh = %{$_[0]};
	foreach my $keys (sort keys %hoh) {
		print "$keys : ";
		foreach my $keys2 (sort keys %{ $hoh{$keys} } ) {
			print "$keys2 = $hoh{$keys}{$keys2}\t";
		}
		print "\n";
	}
}

sub breadthFirstSearch {
	my $_startNode = $_[0];
	my $_maxDepth = $_[1];
	my %_graph = %{$_[2]};
	my @queue;
	
	my $nodeDepth = 0;
	my %nodeDepth = ('node' => $_startNode,
			'depth' => $nodeDepth);
	push (@queue, \%nodeDepth);

	my %visited;

	while (scalar(@queue) != 0) {
		my $nextNodeDepth = pop(@queue);
		my $nextNode = $nextNodeDepth->{'node'};
		my $nextDepth = $nextNodeDepth->{'depth'};
		#print "$nextNode	$nextDepth\n";
		$visited{$nextNode} = $nextDepth;

		if($nextNodeDepth->{'depth'} >= $_maxDepth) {
			next;
		}

		my %neighbors = %{$_graph{$nextNode}};
		foreach my $neighborNode (keys %neighbors) {
			#print "$neighborNode\n";
			if (!exists $visited{$neighborNode}) {
				my %neighbor = ('node' => $neighborNode,
						'depth' => $nextDepth +1);
				push (@queue, \%neighbor);
			}
		}
	}

	return \%visited;
}
	


sub depthFirstSearch {
	my $_node = $_[0];
	my @_visited = @{$_[1]};
	my $_pathWeight = $_[2];
	my $_discoveredPaths = $_[3];
	my %_graph = %{$_[4]};
	my %_targets = %{$_[5]};

	foreach my $neighbor (keys %{ $_graph{$_node} }) {
		my $neighborVisited = 0;
		foreach my $visits (@_visited) {
			if ($neighbor == $visits) {
				$neighborVisited = 1;
			}
		}

		if (!$neighborVisited) {
			my @nextVisited = @_visited;
			push (@nextVisited, $neighbor);
			my $updatedWeight = $_pathWeight + $_graph{$_node}{$neighbor};
		
			if (exists $_targets{$neighbor}) {
				
				my $pathString = '';
				for (my $i = 0; $i < @nextVisited-1; $i++) {
					$pathString .= "$nextVisited[$i] -> ";
				}
				$pathString .= "$nextVisited[-1]";
				$_discoveredPaths->{$pathString} = $updatedWeight;
				
			}
			else {
				&depthFirstSearch($neighbor, \@nextVisited, $updatedWeight, $_discoveredPaths, \%_graph, \%_targets);
			}
		}



	}
}

sub print_paths{
	#Structure of nicely formatted output:
	#Sorted in increasing order of cost
	#1 -> 4 -> 9 -> [0.9000]
	#3 -> 5 -> 9 -> [1.0000]
	#1 -> 9 [2.1003]
	my $_paths = $_[0];
	my $_k = $_[1];
	my $ofstream = $_[2];

	open OUT, ">$ofstream" or die "ERROR: DID NOT OPEN $ofstream!\n";
	
	my $rankCounter;
	foreach my $pathString (sort { $_paths->{$a} <=> $_paths->{$b} } (keys %$_paths) ) {
		my $formattedWeight = sprintf("%.4f",$_paths->{$pathString});
		print OUT "$pathString  [$formattedWeight]\n";
		$rankCounter++;
		if ($rankCounter >= $_k) {
			last;
		}
	}
		

}

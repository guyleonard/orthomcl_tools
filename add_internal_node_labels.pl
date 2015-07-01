#!/usr/bin/perl
use strict;
use warnings;

use Bio::Phylo::IO qw(parse);

#my $newick  = '(((trad,necr),didi),(nagr,(chre,phra)));';

my $filename = 'final_tree_4code_namesA-SAR_OpiOut.tree';
open my $newick_file, '<', $filename or die "Could not open file '$filename' $!";

my $newick_tree = '';
while (my $line = <$newick_file>) {
  chomp $line;
  $newick_tree = $line;
}

my $tree    = parse( '-format' => 'newick', '-string' => $newick_tree )->first;
my $number_of_leaves = 64;
my $counter = $number_of_leaves + 1;

$tree->visit(
    sub {
        my $node = shift;
        $node->set_name( $counter++ ) if $node->is_internal;
    }
);

print $tree->to_newick( '-nodelabels' => 1 );
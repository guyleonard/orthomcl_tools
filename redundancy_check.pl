#!usr/bin/perl -w

use warnings;
use strict;

#main

my @fasta_files = glob("*.fasta");

foreach my $file (@fasta_files) {
	&duplicate($file);
}

sub duplicate {
	my $file = shift;
    open FILE, "<", "$file";
    my @array = <FILE>;
    close FILE;
    my @arr = grep /^>/, @array;
    my %h = ();
    map { $h{$_}++ } @arr;
    my @dupes = grep { $h{$_} > 1 } keys %h;
    #@dupes = grep /^[A-Z,a-z,0-9]>/, @array;
    print "$file\n";
    print "@dupes";
    print "\n";
}
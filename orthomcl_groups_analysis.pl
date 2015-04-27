#!/usr/bin/perl
use strict;
use warnings;

use autodie;
use Bio::Index::Fasta;
use Cwd;               # Gets pathname of current working directory
use File::Basename;    # Remove path information and extract 8.3 filename
use Getopt::Std;

# Analyse the output of OrthoMCL
# 1. Assemble all sequences from each Ortholog group into separate files
# 2. Presence Absense Grid
#   Is the genome represented in the ortholog group? 0 or 1
#   Transposed-grid, phylip-like style.
# 3. Count Grid
#	How many representations (genes) are Present in each ortholog group 0...n

our $EMPTY           = q{};
our $WORKING_DIR     = getcwd;
our $VERSION         = '2015-04-27';
our $GOOD_PROTEINS   = $EMPTY;
our $GROUPS_FILE     = $EMPTY;
our $COMPLIANT_FASTA = $EMPTY;
our $OUTPUT_DIR      = $EMPTY;
our $VERBOSE         = $EMPTY;

####
# Use this command to transpose the output....
# perl -F, -lane 'for ( 0 .. $#F ) { $rows[$_] .= $F[$_] }; eof && print map "$_\n", @rows' presence_absense_grid_no_blem.csv > presence_absense_grid_no_blem_transposed.csv

our $ORTHO_DIR = "/home/cs02gl/Dropbox/projects/AGRP/jeremy_test/six_genomes/test_orthomcl_pipeline_6_genomes";

my $good_proteins = "$ORTHO_DIR\/blast_dir\/goodProteins.fasta";

my $groups_file = "$ORTHO_DIR\/groups\/groups.txt";

my $compliant_fasta_dir   = "$ORTHO_DIR\/compliant_fasta";
my @compliant_fasta_files = glob "$compliant_fasta_dir/*.fasta";

# Commandline Options!
my %options = ();
getopts( 'o:g:p:vh', \%options ) or display_help();

if ( $options{h} ) { display_help(); }
if ( $options{v} ) { print "OrthoMCL Groups Analysis $VERSION\n"; }

if ( defined $options { o } ) {

    $OUTPUT_DIR = $options{o};

    if ( defined $options{g} && defined $options{p} ) {

        $GROUPS_FILE   = $options{g};
        $GOOD_PROTEINS = $options{p};

        collate_sequences( $GROUPS_FILE, $GOOD_PROTEINS, $OUTPUT_DIR );
    }

    #&presence_absense_grid( "$groups_file", \@compliant_fasta_files );

    #&count_grid( "$groups_file", \@compliant_fasta_files );

}
else {

    display_help();
}

sub display_help {
    print "Usage for $VERSION:\nMandatory Input:\n";
    print "\t-o output directory\n";
    print "Collate Sequences:\n\t-g orthoMCL groups.txt\n\t-p goodProteins.fasta\n";
    exit 1;
}

sub presence_absense_grid {
    my $groups_file_list = $_[0];
    my @genomes          = @{ $_[1] };

    foreach my $file (@genomes) {
        $file =~ s/$ORTHO_DIR\/compliant_fasta\///ig;
        $file =~ s/\.fasta//ig;
    }

    #print "Genomes:\n@genomes\n";

    open my $presence_absense_file, '>', "$ORTHO_DIR\/presence_absense_grid.csv";
    print $presence_absense_file "Ortho Group,";
    print $presence_absense_file join( ',', @genomes ), "\n";

    open my $groups_file, '<', "$groups_file_list";
    while ( my $line = <$groups_file> ) {
        my %count_hash = map { $_ => 0 } @genomes;
        my @entries     = split( /\:\s{1}/, $line );
        my $ortho_group = $entries[0];
        my @accessions  = split( /\s+/, $entries[1] );

        #print "Group: $ortho_group\nAccessions\n"; #@accessions\n\n";
        foreach my $accession (@accessions) {
            $accession =~ m/(\w{4})\|/ig;
            $accession = $1;

            #print "$accession ";
            $count_hash{$accession} = 1;
        }

        #print "\n";

        print $presence_absense_file "$ortho_group,";

        # Due to the nature of hash's being "random" we need to sort it on the key value
        foreach my $key ( sort ( keys(%count_hash) ) ) {
            print $presence_absense_file "$count_hash{$key},";
        }
        print $presence_absense_file "\n";
    }
}

sub count_grid {
    my $groups_file_list = $_[0];
    my @genomes          = @{ $_[1] };
    foreach my $file (@genomes) {
        $file =~ s/$ORTHO_DIR\/compliant_fasta\///ig;
        $file =~ s/\.fasta//ig;
    }
    open my $presence_absense_file, '>', "$ORTHO_DIR\/count_list.csv";
    print $presence_absense_file "Ortho Group,";
    print $presence_absense_file join( ',', @genomes ), "\n";
    open my $groups_file, '<', "$groups_file_list";
    while ( my $line = <$groups_file> ) {
        my %count_hash = map { $_ => 0 } @genomes;
        my @entries     = split( /\:\s{1}/, $line );
        my $ortho_group = $entries[0];
        my @accessions  = split( /\s+/, $entries[1] );
        foreach my $accession (@accessions) {
            $accession =~ m/(\w{4})\|/ig;
            $accession = $1;
            $count_hash{$accession} += 1;
        }
        print $presence_absense_file "$ortho_group,";

        # Due to the nature of hash's being "random" we need to sort it on the key value
        foreach my $key ( sort ( keys(%count_hash) ) ) {
            print $presence_absense_file "$count_hash{$key},";
        }
        print $presence_absense_file "\n";
    }
}

sub collate_sequences {

    my $groups_file_list = shift;
    my $good_proteins    = shift;
    my $output_dir       = shift;
    my $alignments_dir   = "$output_dir\/alignments";

    # BioPerl to the rescue!
    my $inx = Bio::Index::Fasta->new( -filename => $good_proteins . ".idx", -write_flag => 1 );

    # pass a reference to the critical function to the Bio::Index object
    $inx->id_parser( \&get_id );

    # make the index
    $inx->make_index($good_proteins);

    open my $groups_file, '<', "$groups_file_list";

    #while ( my $line = <$groups_file> ) {
    foreach my $line (<$groups_file>) {

        my @entries   = split( /\:\s{1}/, $line );
        my @accession = split( /\s+/,     $entries[1] );

        print "\n\/----\nRetreiving: Group: $entries[0] with $#accession seqs\n";

        # We are using bioperl and indexing as fastacmd cannot handle non-NCBI fasta headers
        # See - http://blastedbio.blogspot.co.uk/2012/10/my-ids-not-good-enough-for-ncbi-blast.html

        mkdir $alignments_dir unless -d $alignments_dir;

        open my $file_out, '>>', "$alignments_dir\/$entries[0]\.fasta";
        for ( my $x = 0 ; $x <= $#accession ; $x++ ) {

            my $seq       = $inx->fetch("$accession[$x]");
            my $seqstring = $seq->seq;

            print "$accession[$x],";
            print {$file_out} ">$accession[$x]\n$seqstring\n";
        }
        print "\n\\----\n";
        close $file_out;
    }
}

# remove the chevron from the accession
sub get_id {
    my $header = shift;
    $header =~ /^>(.*)/ism;
    return $1;
}

#!/usr/bin/perl
use strict;
use warnings;

use Cwd;               # Gets pathname of current working directory
use File::Basename;    # Remove path information and extract 8.3 filename

#use Bio::SeqIO;
#use IO::String;
#use Bio::SearchIO;
#use Bio::DB::Fasta;
use Bio::Index::Fasta;

# Analyse the output of OrthoMCL
# 1a. Assemble all sequences from each Ortholog group into separate files
# 1b. At this point I can put them into our pipeline to generate alignments, masks, trees.

# 2. Presence Absense Grid
#	Is the genome represented in the ortholog group? 0 or 1

# 3. Count Grid
#	How many representations (genes) are Present in each ortholog group 0...n

####
# Use this command to transpose the output....
# perl -F, -lane 'for ( 0 .. $#F ) { $rows[$_] .= $F[$_] }; eof && print map "$_\n", @rows' presence_absense_grid_no_blem.csv > presence_absense_grid_no_blem_transposed.csv

our $ORTHO_DIR = "/home/cs02gl/Dropbox/projects/gain_loss/orthomcl";

my $good_proteins = "$ORTHO_DIR\/goodProteins.fasta";

my $groups_file = "$ORTHO_DIR\/groups.txt";

my $compliant_fasta_dir   = "$ORTHO_DIR\/compliantFasta_50";
my @compliant_fasta_files = glob "$compliant_fasta_dir/*.fasta";

# Methods
&collate_sequences( "$groups_file", "$good_proteins" );

#&presence_absense_grid( "$groups_file", \@compliant_fasta_files );

#&count_grid( "$groups_file", \@compliant_fasta_files );

sub presence_absense_grid {
    my $groups_file_list = $_[0];
    my @genomes          = @{ $_[1] };

    foreach my $file (@genomes) {
        $file =~ s/$ORTHO_DIR\/compliantFasta_50\///ig;
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
        $file =~ s/$ORTHO_DIR\/compliantFasta_50\///ig;
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

    my $groups_file_list = $_[0];
    my $good_proteins    = $_[1];

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

        #print "\n---\nRetreiving: Group: $entries[0]\n\nEntry: @accession\n---\n";
        print "\n---\nRetreiving: Group: $entries[0] with $#accession seqs\n---\n";

        # WE CAN NOT USE THESE TOOLS DUE TO USELESS NCBI INSISTENCE ON USING THEIR GODDAMN NAMING SCHEME
        # See - http://blastedbio.blogspot.co.uk/2012/10/my-ids-not-good-enough-for-ncbi-blast.html
        #my $results = `blastdbcmd -db $good_proteins -entry $entry -out $ORTHO_DIR\/alignments\/$entries[0]\.fasta`;
        #my $results = `fastacmd -d $good_proteins -s $entry -o $ORTHO_DIR\/alignments\/$entries[0]\.fasta`;

        # BioPERL to the RESCUE!

        open my $file_out, '>>', "$ORTHO_DIR\/alignments\/$entries[0]\.fasta";
        for ( my $x = 0 ; $x <= $#accession ; $x++ ) {

            #print "Getting: $accession[$x]\n";

            my $seq = $inx->fetch("$accession[$x]");

            my $seqstring = $seq->seq;

            #my $desc = $db->header("$accession[$x]");
            print "$accession[$x],";

            print $file_out ">$accession[$x]\n$seqstring\n";
        }
        print "\n----\n";
        close($file_out);
   }
}

sub get_id {
    my $header = shift;
    $header =~ /^>(.*)/;

    #$header =~ /^>(\w{4}\|.*)/;
    #$header =~ /^>.*\bsp\|([A-Z]\d{5}\b)/;
    $1;
}

# Run makeblastdb for each .fas/.fasta file in a directory
#sub make_blast_db {
#    my $good_proteins = shift;
#    print "makeblastdb - Formating Databases...\n";
#    my $results = `makeblastdb -in $good_proteins -dbtype prot`;
#    #my $results = `formatdb -i $good_proteins -p T -o T`;
#    print "makeblastdb - Finished...\n\n";
#}

        # A couple of other methods, both horribly slow compared to indexing...
###
        #my $gb = Bio::SeqIO->new(
        #    -file   => "<$good_proteins",
        #    -format => "fasta"
        #);
        #my $fa = Bio::SeqIO->new(
        #    -file   => ">$ORTHO_DIR\/alignments\/$entries[0]\.fasta",
        #    -format => "fasta",
        #    -flush  => 0
        #);    # go as fast as we can!
        #while ( my $seq = $gb->next_seq ) {
        #    #Sorry! Here would be with problem, if we use this "if (grep {$_=$seq->id} @genes_name;"
        #    $fa->write_seq($seq) if ( grep { $_ eq $seq->id } @accession );
        #}
###
        #my $db = Bio::DB::Fasta->new("$good_proteins");
        #open my $file_out, '>>', "$ORTHO_DIR\/alignments\/$entries[0]\.fasta";
        #for (my $x = 0; $x <= $#accession; $x++) {
        #  print "Getting: $accession[$x]\n";
        #    my $obj        = $db->get_Seq_by_id("$accession[$x]");
        #    my $seqstring  = $obj->seq;
        #    my $desc       = $obj->header;
        #    #my $seqstring = $db->seq("$accession[$x]");
        #    #my $desc = $db->header("$accession[$x]");
        #    print "$desc\n$seqstring\n";
        #}
        #close($file_out);
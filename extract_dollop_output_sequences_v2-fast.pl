#!/usr/bin/perl
use strict;
use warnings;

use Cwd;               # Gets pathname of current working directory
use File::Basename;    # Remove path information and extract 8.3 filename
use Data::Dumper;

our $WORKING_DIR = getcwd;

# This here script will attempt to parse the ancestral state output from a DOLLOP (phylip) analysis
# It's a horrible format and I had attempted to parse it and so had Adam Sardar (who also found it a nightmare)
# https://github.com/adamsardar/perl-libs-custom/blob/master/Supfam/Phylip_Ancestral_States_Parser.pm
# I chose to just continue with the script I had already as I did not want to start importing tree hashes etc
# But borrow a little wisdom from Adam.

# This script is actually going to first export the information about ancestral states
# as a phylip-style file - which is infinitely easier to parse!

our $EMPTY = q{};

# Files
my $dollop_outfile = "/home/cs02gl/Dropbox/projects/AGRP/jeremy_test/six_genomes/test_orthomcl_pipeline_6_genomes/outfile.old";
my $groups_files   = "/home/cs02gl/Dropbox/projects/AGRP/jeremy_test/six_genomes/test_orthomcl_pipeline_6_genomes/group_list.txt";

# User Variables
my $number_of_states = "9881";

# Other Variables - Do not change
#my $concat_line          = "";
my $number_of_char_lines = int( $number_of_states / 40 );
print "$number_of_char_lines\n";
my $count = 0;

&convert_to_phylip_style($dollop_outfile);

#&report_counts("$dollop_outfile\.phy");

my @groups = get_groups("$groups_files");

#&get_group_from_position("$dollop_outfile\.phy");

our $ALIGNMENTS_DIR = "/home/cs02gl/Dropbox/projects/AGRP/jeremy_test/six_genomes/test_orthomcl_pipeline_6_genomes/alignments";

#get_alignments_for_group();

sub get_group_from_position {

    my $phylip_file = shift;
    open my $results_infile, '<', "$phylip_file";

    my $nodes_dir = "$WORKING_DIR\/nodes";

    mkdir $nodes_dir unless -d $nodes_dir;

    foreach my $line (<$results_infile>) {

        # skip blank lines, not needed but speeds up location finding a mite
        next if ( $line =~ m/^[\s|\t]+/ );

        # Lines read in will look similar to below...gap may be space or tab
        # root_1    110100111.111..0

        my @split_line = split( /[\s\t]/, $line );
        my $node       = $split_line[0];
        my $results    = $split_line[1];
        $results =~ s/yes//g;

        # split the results into chars - each char is a position within the
        # groups array...
        my @results_array = split( //, $results );

        my @gain = ();
        my @loss = ();
        my @core = ();

        # There was an "out-by-two" error in the nodes being output
        # i had to +5 to the X position to correct for this
        # 3 from the 'yes' stripped from the output
        # 1 to count itself
        # and 1 to counter the zero-start position of an array e.g. 0..9 = 10

        for ( my $x = 0 ; $x <= $#results_array ; $x++ ) {

            if ( $results_array[$x] eq "." ) {
                push( @core, "$groups[$x+5]" );
            }
            elsif ( $results_array[$x] eq "1" ) {
                push( @gain, "$groups[$x+5]" );
            }
            elsif ( $results_array[$x] eq "0" ) {
                push( @loss, "$groups[$x+5]" );
            }
        }

        # $# need + 1
        open my $general_report, '>>', "$nodes_dir\/general_report\.txt";
        print $general_report "$node\t" . ( $#core + 1 ) . "\t" . ( $#gain + 1 ) . "\t" . ( $#loss + 1 ) . "\n";
        print "$node\t" . ( $#core + 1 ) . "\t" . ( $#gain + 1 ) . "\t" . ( $#loss + 1 ) . "\n";

        my $node_dir = "$nodes_dir\/$node";
        mkdir $node_dir unless -d $node_dir;

        open my $core_report, '>', "$node_dir\/core\.txt";
        print $core_report "@core";
        open my $gain_report, '>', "$node_dir\/gain\.txt";
        print $gain_report "@gain";
        open my $loss_report, '>', "$node_dir\/loss\.txt";
        print $loss_report "@loss";

    }
}

sub get_alignments_for_group {

    my $nodes_dir = "$WORKING_DIR\/nodes";

    # Get a list of everything in dir - also reads in files...
    my @node_subdirs = glob "$nodes_dir/*";

    foreach my $x (@node_subdirs) {

        print "Working on $x\n";

        # Test for directory
        if ( -d $x ) {

            ## Gains
            #
            # Make a directory for alignments to be placed
            mkdir "$x\/gain" unless -d "$x\/gain";
            print "\tMaking gain directory\n\tCopying";

            # Read in file contents to array - it should only be one line.
            open my $gain_in, '<', "$x\/gain\.txt";
            my @gain_line = split( /\s+/, <$gain_in> );
            close($gain_in);

            # Foreach array element, get the corresponding file
            # from the alignments dir...
            foreach my $y (@gain_line) {
                system "cp $ALIGNMENTS_DIR\/$y\.fasta $x\/gain";
                print "+";
            }
            print "\n\n";

            ## Losses
            #
            # Make a directory for alignments to be placed
            mkdir "$x\/loss" unless -d "$x\/loss";
            print "\tMaking loss directory\n\tCopying";

            # Read in file contents to array - it should only be one line.
            open my $loss_in, '<', "$x\/loss\.txt";
            my @loss_line = split( /\s+/, <$loss_in> );
            close($loss_in);

            # Foreach array element, get the corresponding file
            # from the alignments dir...
            foreach my $y (@loss_line) {
                system "cp $ALIGNMENTS_DIR\/$y\.fasta $x\/loss";
                print "-";
            }
            print "\n\n";

            # Default this to off, too many each time to cp
            ## Core
            #
            # Make a directory for alignments to be placed
            #mkdir "$x\/core" unless -d "$x\/core";
            #print "\tMaking core directory\n\tCopying";

            # Read in file contents to array - it should only be one line.
            #open my $core_in, '<', "$x\/core\.txt";
            #my @line = split (/\s+/, <$core_in>);
            #close($core_in);

            # Foreach array element, get the corresponding file
            # from the alignments dir...
            #foreach my $y (@line) {
            #    system "cp $ALIGNMENTS_DIR\/$y\.fasta $x\/core";
            #    print ".";
            #}
            #print "\n\n";
        }
    }

}

sub get_groups {
    my $groups_files = shift;
    open my $group_infile, '<', "$groups_files";
    my @groups = ();
    while (<$group_infile>) {
        my $line = $_;
        chomp($line);
        push( @groups, $line );
    }
    return @groups;
}

sub report_counts {

    my $dollop_phylip = shift;
    print "\tCounting\n";
    open my $dollop_phylip_in, '<', "$dollop_phylip";

    open my $report, '>', "$dollop_phylip\_report\.txt";
    print $report "Node 1\tNode2\tShared\tLoss\tGain\n";

    foreach my $line (<$dollop_phylip_in>) {
        chomp($line);
        next if ( $line =~ m/^\s+/ );

        my @line_array = split( /\t/, $line );

        my @nodes = split( /\_\_/, $line_array[0] );

        # Count number of '.', '0', and '1' in concat_line.
        # '.' = shared, '0' = loss, '1' = gain
        my $dot_count  = () = $line_array[1] =~ m/[\.]/g;
        my $one_count  = () = $line_array[1] =~ m/[1]/g;
        my $zero_count = () = $line_array[1] =~ m/[0]/g;

        print $report "$nodes[0]\t$nodes[1]\t$dot_count\t$zero_count\t$one_count\n";
    }
}

sub convert_to_phylip_style {

    my $dollop_outfile = shift;
    my $concat_line    = $EMPTY;

    print "Parsing: $dollop_outfile\n";

    open my $dollop_outfile_read, '<', "$dollop_outfile";
    open my $dollop_phylip,       '>', "$dollop_outfile\.phy";
    print $dollop_phylip " 0 $number_of_states\n";    #need to add number of nodes...

    my $correct_location = 0;

    foreach my $line (<$dollop_outfile_read>) {
        chomp($line);

        # skip blank lines to speed up location finding
        next if ( $line =~ m/^[\s\t]*$/);

        #if ( $line =~ m/means same as in the node below it on tree/g ) {
        if ( $line =~ m/( . means same as in the node below it on tree)/g ) { 
            $correct_location = 1;
        }

        if ( $correct_location == 1 ) {

            # if the lines starts with root or 2 spaces followed by a digit
            if ( $line =~ m/^(root|\s{2}\d+)\s+/gism ) {

                print ".";
                $count = 0;

                # Split line on white space
                my @node = split( /\s+/, $line );

                # Note the different positions in the array, as we split on space
                if ( $node[0] eq 'root' ) {
                    print $dollop_phylip "$node[0]\_\_$node[1]";
                    # Add the first line of dots, 0s and 1s to the concatenated line
                    $concat_line = $node[3] . $node[4] . $node[5] . $node[6] . $node[7] . $node[8] . $node[9] . $node[10];
                }
                else {
                    print $dollop_phylip "$node[1]\_\_$node[2]";
                    # Add the first line of dots, 0s and 1s to the concatenated line
                    $concat_line = $node[4] . $node[5] . $node[6] . $node[7] . $node[8] . $node[9] . $node[10] . $node[11];
                }
            }
            else {

                $line =~ s/\s//g;    # remove white spaces
                $concat_line = $concat_line . $line;
            }

            if ( $count == $number_of_char_lines ) {
                print $dollop_phylip "\t$concat_line\n";

                # empty the concatened lines for the next set
                $concat_line    = $EMPTY;
            }
            $count++;
        }       
    }
    print "\n";
}
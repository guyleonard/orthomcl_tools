# orthomcl_tools
A couple of tools I made for post orthoMCL/orthAgogue down-stream analysis. The first tool computes a set of CSV files with presence/absence or 'count' information relating to your orthologue groups and taxa. The second tool allows you to take the output of a DOLLOP analysis along with the orthologue groups and plot them to locations on the tree topology as a list and set of alignments for each node. The final two scripts allow you to plot this information in a image created in R using ggtree.

## orthomcl_groups_analysis.pl

Run this script to:
1. Collate all sequences from each ortholog group into separate FASTA files.
2. Generate a Presence Absense Grid
  1. Is the taxa/genome represented in the ortholog group? 0 or 1
  2. Also generates a transposed grid in a phylip-like format which is useful for Dollo Parsimony / ML analyses.
3. Generate Count (tally) Grid
  1. How many representations (genes) are present in each ortholog groups? 0...n

You will need these files:
```
- from orthomcl: goodProteins.fasta
- from orthomcl: groups.txt
- from orthomcl: compliantFasta directory
```

You will need to convert the presence_absence_grid.csv file into a phylip-like format. To do this you must do four things:
 1. transpose the data,
 2. remove the first line (header information),
 3. insert spaces between the 'taxa' names and 0/1s and
 4. add number of taxa and number of 'sites' to the top of the file.

To do this:
```
 1. perl -F, -lane 'for ( 0 .. $#F ) { $rows[$_] .= $F[$_] }; eof && print map "$_\n", @rows' presence_absense_grid.csv > presence_absense_grid_transposed.csv
 2. perl -ni -e 'print unless $. == 1' presence_absense_grid_transposed.csv
 3. cp presence_absense_grid_transposed.csv presence_absense_grid.phy
 4. sed -i 's/\(\w\{4\}\)\(.*\)/\1      \t\2/g' presence_absense_grid.phy
 5. Open in your favourite text editor and add " XX YYYY", where XX = number of taxa, and YYYY = number of ortholog groups
```
Before running extract_dollop_output_sequences_v2-fast.pl you will need to run 'dollop' from the PHYLIP package. This requires the phylip-like file you created in the step above, and possibly a user-specified tree topology.

## extract_dollop_output_sequences_v2-fast.pl

A program to take the "outfile" from a PHYLIP DOLLOP run and parse the output in to a more useful format. It requires two pieces of information: 1) the "outfile", and 2) the number of states being tested. For running with OrthoMCL data you will also need: 1) a list of orthogroup names, and 2) a directory with the *.fasta sequences of each orthogroup
The program outputs four different options:

1. A phylip-like output file, parsed from the state information in the "outfile" from dollop.
2. A tab-separated report listing the node and number of gain, loss and '*core*' states in either of two styles:
  1. 'New-style': This makes the node column one value, which is represented by a internal node number or leaf label. The number is the node number starting from '1' on the first leaf.
  2. 'Old-style': This makes the node column into two values, corresponding to positions across the tree using the node values from the dollop outfile.
3. This outputs a directory, in the 'old-style' format (e.g. root__1 or 3__label), making a directory for each node transition along with three files. One each for loss, gain and core.
4. This copies the *.fasta files from a location for each of the lists in the previous directory, to their appropriate directory.

```
Mandatory Input:
	-i Dollop outfile
	-s Number of states (orthogroups)
	-o Output Directory
Other Options (one or all required):
outfile to phylip-like:
	-c Convert to Phylip-like File (not needed if using -p)
Report Tables:
	-n New-style Report (includes internal tree node numbers)
	-r Old-style report (between-nodes from dollop outfile)
Structured Lists:
	-l Lists of Core/Losses/Gains (requires -g)
	-g Ortholog Group List
	-p phylip-like file (unless -c)
Sequence Collation:
	-f Get .fasta files for groups from directory (requires -d)
	-d The *.fasta directory
	-n Nodes directory (unless -l)
	-x Exclude "core" marked ortholog groups
e.g. Equivalent: program.pl -i input -s number -o output -cr -l -g list.txt or program.pl -i input -s number -o output -rl -g list.txt -p phylip-like.phy

```

## redundancy_check.pl


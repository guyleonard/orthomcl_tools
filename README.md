# orthomcl_tools
A couple of tools I find indispensable for post orthomcl down-stream analysis.

## extract_dollop_output_sequences_v2-fast.pl

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

## orthomcl_groups_analysis.pl
## redundancy_check.pl


library(ggplot2)
library(ggtree)

user_report <- read.csv("outfile.phy_newstyle_report.txt", sep = "\t", header=TRUE, as.is=1, row.names=NULL)
user_tree <- read.tree("tree_4code_internal_labels.tree")
user_p <- ggplot(user_tree, aes(x, y), ladderize=TRUE) + geom_tree() + theme_tree() + geom_tiplab(size=3, align=TRUE, color="purple", x=13) + xlab("") + ylab("") + geom_text(aes(label=Gain, x=branch), size=3, color="springgreen4", vjust=-0.6, subset=.(!isTip)) + geom_text(aes(label=Gain), size=3, color="springgreen4", hjust=0, subset=.(isTip), x=13.5) + geom_text(aes(label=Loss, x=branch), size=3, color="firebrick3", vjust=1.3, subset=.(!isTip)) + geom_text(aes(label=Loss), size=3, color="firebrick3", hjust=0, subset=.(isTip), x=14) + geom_text(aes(label=node), size=2, hjust=-1.5, subset=.(!isTip), color="grey") + scale_x_continuous(expand = c(.1, .2))
user_p <- user_p %<+% user_report
print(user_p)
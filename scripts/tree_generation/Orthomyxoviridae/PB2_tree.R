setwd("C:/Users/nadja/Documents/LaptopAsus/PhD/Chapter_4/")

##################################################
# Dependencies
##################################################

library (ggplot2)
library(ggrepel)
library (phytools)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree", force = TRUE) 
library(BiocManager)
#install("ggtree", force = TRUE)
#creating tree
library(ggtree)
library(tidytree)
library(ape)
library(treeio)

nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")
##################################################

# following files are needed for this workflow:
# contree file from iqtree (or other nex/nwk tree files- with bootstrap values)
# txt file with Orthomyxovirus references with columns accessions and virus name 

## reading in tree
tree <- read.tree("PB2_new_tree.contree") # created with iqtree
## replacing the label descriptions
replace_words <- read.csv("orthomyxo_references.txt", sep = "\t", header = FALSE, col.names = c("Accession", "Virus"))
for(j in seq_along(replace_words$Accession)){
  tree$tip.label <- gsub(replace_words$Accession[j], replace_words$Virus[j], tree$tip.label)
}
tree <- as.phylo(tree)
#outgroup rooting
nodes <- as_tibble(tree)
#tree <- midpoint.root(tree)
tree <- root(tree, node = 317, edgelabel = TRUE)
#### converting the bootstrap value outside ggtree
q <- ggtree(tree) + xlim(NA, 7.5) + expand_limits(y = 100) + coord_cartesian(clip="off")# + geom_text2(aes(label=node), hjust=-.3, size=2)  
#q <- flip(q, 330, 493)
#q <- rotate(q, 463)
#q <- collapse(q, node=487, 'min')
#q <- collapse(q, node=479, 'max')
#q <- collapse(q, node=476, 'max')
#q <- collapse(q, node=473, 'max')
d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 95,] # for significant bootstrap values
d$label <- replace(d$label, d$label > 94,'*') 
q <- q + geom_text(data=d, aes(label=label), hjust = 0.0, vjust = 0.1, size = 3.0) 

show (q)
## annotation location and samples
tipcategories = read.csv("orthomyxo_genomes.txt", 
                         sep = "\t",
                         col.names = c("EVE", "Organism"), 
                         header = FALSE, 
                         stringsAsFactors = FALSE)

dd = as.data.frame(tipcategories)

q <- q %<+% dd +
  geom_tiplab(aes(), hjust = -0.2, size=2.5, alpha=.75) +
  geom_tippoint(aes(color=Organism), size=2, alpha=.75) +
  scale_color_manual(values = c("Anopheles" = "#E31A1C",  "Aedes albopictus" = "#1F78B4", "Aedes aegypti" = "orange"), na.translate=F, name = "", labels = c("Anopheles",  "Aedes aegypti", "Aedes albopictus")) +
  theme(legend.position= c(0.15, 0.5), legend.text = element_text(size = 13), title = element_text(size = 22)) +
  theme(legend.key=element_blank(),legend.background=element_blank()) +
  geom_treescale(x = 0, y = 0, offset = 1, width = 0.5) 
show(q)

### closer look at specific clade for checking
#viewClade(tree_view = NULL, 335, xmax_adjust = 0)
#show (q)
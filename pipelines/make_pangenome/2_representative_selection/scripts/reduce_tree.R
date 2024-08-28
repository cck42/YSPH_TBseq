library(ape)

# reduce_tree.R
args <- commandArgs(trailingOnly = TRUE)

# load your tree
tree <- read.tree(args[1])

# we are removing the refseq and the ougroup from the tree 
# so that parnas only selects reps from the rest of the tree 
# since the refseq and outgroup are automatically included in the tree
tree_modified <- drop.tip(tree, c(args[3], args[4]))

# save the tree to a new file 
write.tree(tree_modified, file=args[2])

library(admixtools)
library(ape)



prefix = "path/to/eigendtrat/files"

my_f2_dir = "path/to/results/folder"

#making
extract_f2(prefix, my_f2_dir, maxmiss=1, overwrite=T)

#reading
f2_blocks = f2_from_precomp(my_f2_dir,remove_na = FALSE)
dim(f2_blocks)
count_snps(f2_blocks)


# compute average of f2 values on all the blocks
f2_mean <- apply(f2_blocks, c(1, 2), mean, na.rm = TRUE)

# change names
rownames(f2_mean) <- rownames(f2_blocks[, , 1])
colnames(f2_mean) <- colnames(f2_blocks[, , 1])


round(f2_mean[1:5, 1:5], 3)

# save f2 matrix for futre anaysis
write.table(f2_mean, file="f2_matrix.txt", quote=F, sep="\t")

# convert to dist class object
f2_dist <- as.dist(f2_mean)

# make neighbour-joining tree
nj_tree <- nj(f2_dist)

# root the tree
nj_rooted <- root(nj_tree, outgroup = "Mbuti", resolve.root = TRUE)
plot(nj_rooted, main = "Rooted NJ tree")

# save the tree
write.tree(nj_rooted, file = "f2_tree.nwk")

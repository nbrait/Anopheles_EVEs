
# Install and load the ape package
install.packages("ape")
library(ape)

# Read the Newick format tree from a file
tree <- read.tree("Anopheles_species.nwk.")

# Plot the cladogram
plot(tree)

# Optional: Add a title or customize the plot
title("Cladogram of Anopheles Species")

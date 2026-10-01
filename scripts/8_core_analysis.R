# The core analysis is applied both across and within taxa, and consists of two 
# distinct stages:
# 1. Pairwise PERMANOVA between GMS clusters on protein shapes (centroid diffs)
# 2. Pairwise betadisper between GMS clusters on protein shapes (dispersion diffs)

# A tile plot illustrating significance of difference between pairwise cluster 
# comparisons is also planned. It will be an nxn tile plot where n is the number of 
# GMS clusters. The lower triangle should contain betadisper results which should be
# coloured with a blue hue. The upper triangle should contain PERMANOVA results
# which should be coloured with a red hue. Intensity of hue should match
# significance. Ideally, tiles matching the same cluster should be greyed out.

# The above tile plot should also be generated for within taxa analyses, and
# tiles that match GMS clusters that are not present in the taxa should be
# greyed out.

# TODO
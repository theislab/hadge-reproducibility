library(SoupX)
# Load data and estimate soup profile
sc = load10X("/Volumes/NAME/deconvolution-project/outs")

# Read the CSV file as a data frame
my_data <- read.csv("/Volumes/NAME/deconvolution-project/gx12_clusters.csv")

# Extract the vector column and convert it to a numeric vector
cluster_labels <- as.numeric(my_data$soupx_groups)

# View the contents of the vector
print(cluster_labels)

sc = setClusters(sc,cluster_labels)
sc = autoEstCont(sc)
out = adjustCounts(sc)

class(out)

library(Matrix)
writeMM(out, "/Volumes/NAME/deconvolution-project/gx12_decontaminated_matrix.mtx")

library(DropletUtils)
DropletUtils:::write10xCounts("/Volumes/NAME/deconvolution-project/gx12_decontaminated_matrix", out)

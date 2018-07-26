library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendsort)

setwd("/Users/stacey/Downloads/Raw and normalized data/csv format")

# Files to work with
norm_nut <- read.csv(file = "Normalized data of DMSO_27UDP + nutlin.csv")
norm_x <- read.csv(file = "Normalized data of DMSO_27UDP no nutlin.csv")

# raw_nut <- read.csv(file = "Raw data of DMSO_27UDP +  nutlin.csv")
# raw_x <- read.csv(file = "Raw data of DMSO_27UDP no nutlin.csv")

# Check out the file

head(norm_nut)
colnames(norm_nut)

# Selecting the columns that we want in our new data frame. In this case the
# first column with the gene name (column 1) and the count columns
# (columns 20:25).
df <- norm_nut[, c(1, 20:25)]
df2 <- norm_x[, c(1, 20:25)]

# Renaming the columns
names <- paste(rep(c("Nut_27", "Nut_DMSO"), each = 3), 1:3, sep ="-")
colnames(df)[2:7] <- names

names2 <- paste(rep(c("X_27", "DMSO"), each = 3), 1:3, sep ="-")
colnames(df2)[2:7] <- names2

# Removing the first two rows
df <- df[-1, ]
df <- df[-1, ]

df2 <- df2[-1, ]
df2 <- df2[-1, ]

# Merge
merged <- merge(df,
                df2,
                by.x = "Probe.Name")

# Changing the type of data from "character" to "numeric"
my_data <- apply(merged[, -1],
                 2,
                 as.numeric)

# Name the rows with the gene name
rownames(my_data) <- merged[, 1]

### Remove Neg and Pos values

pos <- grep("^POS_", row.names(my_data))
my_data <- my_data[-pos, ]

neg <- grep("^NEG_", row.names(my_data))
my_data <- my_data[-neg, ]

# Removing housekeeping genes individually
my_data <- my_data[rownames(my_data) != "GAPDH", ]
my_data <- my_data[rownames(my_data) != "ACTB", ]
my_data <- my_data[rownames(my_data) != "PGK1", ]

#### Filtering for minimum read counts

# The minumum number of reads I want
min_reads <- 5

# The minimum number of samples that must contain the minimum read
# value (min_reads)
min_samples <- 12

# Using apply() to identify the elements in the dataframe that match our
# conditions.
read_filter <- apply(my_data,
                     1,
                     function(samples) length(which(samples >= min_reads)) >= min_samples)

# Subsetting the elements that were identified using our conditions and storing
# it in a new object/dataframe called filtered_data
filtered_data <- my_data[read_filter, ]

# Look at the top of the data and stats.
head(filtered_data)
nrow(filtered_data)

#### Fixing data for heatmap

data <- filtered_data

Nut_27_ave <- rowMeans(data[, 1:3])
Nut_dmso_ave <- rowMeans(data[, 4:6])
X_27_ave <- rowMeans(data[, 7:9])
dmso_ave <- rowMeans(data[, 10:12])

# Adding averages columns to the dataframe
data <- cbind(data, Nut_27_ave, Nut_dmso_ave, X_27_ave, dmso_ave)

# Ordering dataframe by descending value
data <- data[order(-dmso_ave), ]

#### Save table in csv file

write.csv(data, paste(dir, "norm_data_averages.csv", sep = "/"))

#### Heatmap of average values for each treatment

# Adjusting breaks in colorbar
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(data, n = 100)

mat_cluster_rows <- sort_hclust(hclust(dist(data)))
matrix <- data.frame(values = as.numeric(data[, 1]))
dat_colors <- data.frame(
  xmin = mat_breaks[1:(length(mat_breaks) - 1)],
  xmax = mat_breaks[2:length(mat_breaks)],
  ymin = 0,
  ymax = max(density(matrix$values, bw = "SJ")$y),
  fill = viridis(length(mat_breaks) - 1),
  stringsAsFactors = FALSE
)

# Setting directory for saved file
dir <- "/Users/stacey/Downloads/Raw and normalized data/R_data"
dir.create(dir)

pdf(file = paste(dir, "averages_heatmap.pdf", sep = "/"))
heatmap <- pheatmap(
  data[, 13:16],
  color = viridis(length(mat_breaks) - 1),
  breaks = mat_breaks,
  border_color = NA,
  fontsize = 8,
  cluster_cols = FALSE,
  cluster_rows = FALSE
)
print(heatmap)
dev.off()

#### Boxplots of genes in data matrix

library(wesanderson)
library(ggplot2)

print.dir <- paste(dir, "figures", sep = "/")
dir.create(print.dir)

groups <- c("Nut_27", "Nut_DMSO", "X_27", "DMSO")
groups <- factor(groups)
group <- rep(groups, each = 3)

for (index in 1:nrow(data)) {
  gene.name <- rownames(data)[index]

  df <- as.data.frame(data[index, 1:12])
  df$group <- group
  colnames(df) <- c("counts", "group")

  color <- wes_palette(n = 4, name = "Royal1")
  dodge <- position_dodge(width = 0.75)

  plot <- ggplot(df,
                 aes(x = group, y = counts, fill = group))

  boxplot <- plot +
    geom_boxplot() +
    geom_point(position = dodge,
               color = "black",
               size = 1.5)


  adj_boxplot <- boxplot +
    ggtitle(gene.name) +
    labs(x = "Treatment",
         y = "Counts") +
    scale_fill_manual(values = color,
                      name = "Treatment") +
    theme(
      plot.title = element_text(size = 15,
                                face = "bold"),
      axis.title.x = element_text(
        size = 12,
        face = "bold",
        colour = "black",
        vjust = 0
      ),
      axis.text.x = element_text(size = 10,
                                 colour = "black"),
      axis.title.y = element_text(
        size = 12,
        face = "bold",
        colour = "black",
        hjust = 0.5,
        vjust = 1
      ),
      axis.text.y = element_text(
        size = 10,
        colour = "black",
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10, face = "bold"),
      panel.grid.major = element_line(colour = "gray75", linetype = "2222"),
      panel.grid.minor = element_line(colour = "white", linetype = "1111"),
      axis.line = element_line(size = .5, colour = "gray10")
    )

  pdf(paste(print.dir, paste0(gene.name, ".pdf"), sep = "/"))
  print(adj_boxplot)
  dev.off()
}

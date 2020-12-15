library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggdendro)
library(RColorBrewer)

th <- theme_classic() +
	theme(panel.grid.major=element_line(size=0.5, linetype="solid", colour="#eeeeee"),
	panel.grid.minor=element_line(size=0.25, linetype="solid", colour="#eeeeee"))
theme_set(th)

## Read data and prepare data

# read otu table and cast counts to numeric, cor() and others require numeric data
otu <- read.delim("1335.vdb.tab", header=TRUE) %>%
	mutate_if(is.integer, as.numeric)

### Preliminary analysis and data manipulation

# for each row of the otu table
tax <- apply(otu, 1, function(x)
	# select OTUConTaxonomy column
	x[2] %>%
		# split strings and get the result, returned inside a list
		strsplit(";") %>%
		# subset the list
		"[["(1) %>%
		# remove assignment score from the name
		sub(pattern="\\(.*", replacement="")
)
# transpose so OTUs are again the rows
tax <- t(tax)
# bind the new phylum and taxon column to the original table
otu <- cbind(otu, data.frame(Genus=tax[,2], Taxon=tax[,6]))

otumat <- otu %>%
	# the counts are the only numeric columns
	select_if(is.numeric) %>%
	as.matrix

# this same construct will be used later:
# take the numeric count matrix
taxtotals <- otumat %>%
	# for each of its columns (sites)
	apply(2, function(x)
		# for each unique taxon (rows)
		sapply(unique(otu$Taxon), function(t)
			# sum its abundance
			sum(x[otu$Taxon == t])
		)
	# this returns total abundance per taxon per site
	) %>%
	# sum abundance across sites => table with total abundance per taxon
	apply(1, sum)
# order taxons by descending abundance
taxtotals <- taxtotals[order(-taxtotals)]
# then extract the top 19 taxons (removing unclassified)
toptax <- taxtotals[2:20]


## Phylum composition by site

### Color palette

# number of unique phylum, the number of colors needed
ngenus <- length(unique(otu$Genus))
# extract the colors of the Set1 brewer color palette
# this gives a harmless warning since we request too many colors
pal.g <- brewer.pal(ngenus, "Set1")
# interpolate for the remaining elements
# colorRampPalette returns a function which we use directly to extract our colors
cmap.g <- colorRampPalette(pal.g, bias=1, interpolate="spline")(ngenus)

### Plot

# stacked percent barchart
bar.g <- otu %>%
	# get count columns and pivot to long format (similar to melt in reshape2)
	# tidyverse matches() is used to select the columns via a regex
	# site column will contains the sample site, and abundance the abundance values
	pivot_longer(cols=matches("^[a-z]", ignore.case=FALSE), names_to="site", values_to="abundance") %>%
	# x=site, y=abundance and use phylum as group colors for the bar's fill
	ggplot(aes(site, abundance, fill=Genus)) +
		# output percent stacked barchart
		geom_bar(position="fill", stat="identity") +
		# rotate xlab and remove x title
		theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1),
			axis.title=element_blank()) +
		# manually set the colors for each phylum from the palette
		scale_fill_manual(values=cmap.g[order(unique(otu$Genus))])
# save it
ggsave("plots/stacked.genus.by.site.png", bar.g, width=12, height=10)

## Taxonomic composition by site

### Creating the color palettes

# filter the taxonomy table by the names of the top taxa, get phylum (2) and taxon (6) columns
corr <- tax[tax[,6] %in% names(toptax),][,c(2,6)]
# only take unique rows, then convert to data frame and set appropriate column names
corr <- unique(corr) %>%
	as.data.frame %>%
	# "family" seems more appropriate than "genus" here
	rename(Family=V1, Taxon=V2)

# take our correspondence table
corr <- corr %>%
	# convert the family column to an ordered factor according to abundance and importance
	# the levels are from what was obtained above
	mutate(Family=factor(Family,
		levels=c("Proteobacteria", "Actinobacteria", "Firmicutes", "Bacteroidetes"),
		ordered=TRUE)) %>%
	# sort by family then taxon
	arrange(Family, Taxon) %>%
	# now freeze the new taxon order
	mutate(Taxon=factor(Taxon, levels=unique(Taxon), ordered=TRUE))

# add a new column with the colors for the ordering
corr$col <- c(
	# blue, red, green and violet/pink gradients
	# like before, we directly call the function returned with the appropriate number of colors we need
	colorRampPalette(c("lightblue", "darkblue"))(8),
	colorRampPalette(c("#ffff66", "red", "darkred"))(5),
	colorRampPalette(c("lightgreen", "darkgreen"))(4),
	colorRampPalette(c("violet", "blue"))(8)[1:2]
)

### Stacked percent barchart

bar.t <- otu %>%
	# remove columns we don't care about
	select(-Genus, -starts_with("OTU")) %>%
	# extract rows for only top taxa
	filter(Taxon %in% names(toptax)) %>%
	# convert to long format for ggplot as before
	pivot_longer(cols=matches("^[a-t]", ignore.case=FALSE), names_to="site", values_to="abundance") %>%
	# convert to ordered factor as with the correspondence table, using its order
	# counterintuitive and timeconsuming to figure out, but very important to avoid any discrepancies...
	mutate(Taxon=factor(Taxon, levels=levels(corr$Taxon), ordered=TRUE)) %>%
	ggplot(aes(site, abundance, fill=Taxon)) +
		# stacked percent barchart
		geom_bar(position="fill", stat="identity") +
		# set our colors
		scale_fill_manual(values=corr$col) +
		# rotate x legend, remove x axis label
		theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1),
			axis.title=element_blank())
# save the plot
ggsave("plots/stacked.taxon.by.site.png", bar.t, width=12, height=10)

## Phylum correlation map

### Preparing the data

# use the abundance matrix computed at the beginning, cor() needs numerics
gentotals <- otumat %>%
	# for each site (column)
        apply(2, function(x)
		# for each OTU (row),
		# sum abundances of OTUs assigned to each individual phylum
                sapply(unique(otu$Genus), function(g)
                        sum(x[otu$Genus == g])
                )
        )

# calculate a frequency table 
genfreq <- table(otu$Genus)
# filter the most infrequent phylum
gentotals <- gentotals[!rownames(gentotals) %in% names(genfreq)[genfreq < 5],]

### Computing correlation and ordering

gencor <- cor(t(gentotals))

# perform the clustering
hcl <- hclust(as.dist(1 - gencor), method="ward.D2")
# extract order
horder <- rownames(gentotals)[hcl$order]
# convert to data appropriate for ggplot2, which gives us segment coordinates for drawing it directly on top
dhcl <- dendro_data(hcl)
# reorder the matrix
gencor <- gencor[horder,horder]

# lower.tri selects the matrix' lower triangle, and we just remove all the data
gencor[lower.tri(gencor)] <- NA

## Correlation map

# use the ordered correlation matrix
cormap <- gencor %>%
	# convert to dataframe for ggplot
	as.data.frame %>%
	# get the rownames to use for grouping
	rownames_to_column("grp") %>%
	# convert to long format, remove missing values (upper triangle)
	pivot_longer(-grp, names_to="key", values_to="v", values_drop_na=TRUE) %>%
	# make sure the order is preserved like before, otherwise elements may be shifted seemingly at random
	mutate(grp=factor(grp, levels=horder, ordered=TRUE),
		key=factor(key, levels=horder, ordered=TRUE)) %>%
	# pass it directly to ggplot; label will be the rounded correlation value within each cell
	ggplot(aes(key, grp, fill=v, label=round(v,2))) +
			# tilemap background
			geom_tile(color="white") +
			# set a blue to red gradient for [-1;1] correlation values
			scale_fill_gradient2(
				low="blue", high="red", mid="white",
				midpoint=0,
				limit=c(-1,1), space="Lab",
				name="Pearson\nCorrelation") +
			# rotate x axis for easier viewing, remove axis title
			theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1),
				axis.title=element_blank()) +
			# make sure we get a square matrix
			coord_fixed() +
			# add correlation values
			geom_text(color="black", size=3.5) +
			# draw the dendrogram; we offset y coordinates to place it above the tilemap
			geom_segment(data=segment(dhcl), aes(x, y+nrow(gencor)+0.5, xend=xend, yend=yend+nrow(gencor)+0.5), inherit.aes=FALSE)
# save it
ggsave("plots/cormap.pdf", cormap, width=12, height=10)


## Principal component analysis

# compute pca on numeric matrix, scale the variables to unit variance as recommended
pca <- prcomp(otumat, scale=TRUE)
# percent of variance explained by each axis
pca$sdev ^ 2 / sum(pca$sdev ^ 2)
# 40.6% and 28.6% for PC1 and PC2

# render base graphics plot as pdf
pdf("plots/pca.pdf", 12, 10)
# plot the pca
biplot(pca,
	# draw arrows for each eigenvector
	var.axes=TRUE,
	xlabs=rep("*", nrow(otumat)), cex=c(0.9,1.1),
	xlim=c(-0.9,0.2), expand=50,
	col=c("black", "#dd0000"))
dev.off()

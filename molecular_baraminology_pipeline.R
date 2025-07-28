# Load required libraries
library(rentrez)
library(Biostrings)
library(msa)
library(seqinr)
library(stringr)
library(ggplot2)
library(cluster)
library(factoextra)

setwd("C:\\Users\\csmat\\OneDrive\\Desktop\\CreationScience\\Talks\\CRSconferences\\CRS2025\\workshop\\Elephants")

# Define list of Latin species names
# species_list <- c("Homo sapiens", "Mus musculus", "Canis lupus")
# species_list <- c("Lilium davidii", "Notholirion campanulatum", "Fritillaria sinica")
species_list <- read.table(choose.files(filters = Filters[c("txt"),]),sep="\t",check.names=FALSE)
species_list <- species_list$V1

# Search for:
# 1. organelle genome
#    a. mitogenome (clustalw)
#    b. clp genome (MAFFT)
# 2. protein (clustalw)

mito <- FALSE
clp <- TRUE
prot <- FALSE
protein_name <- "TP53" # "NULL"

# Initialize empty list to hold sequences
sequence_list <- list()

# Fetch sequences
for (species in species_list) {
  cat("Fetching:", species, "\n")
  
  species_name <- str_replace_all(species, " ", "_")
  
  # Search term to get the mitochondrial genome
  if (mito == TRUE) {
    search_term <- paste0(species, "[Organism] AND mitochondrion[Title] AND complete genome")
    search_results <- entrez_search(db = "nucleotide", term = search_term, retmax = 1)
  } else if (clp == TRUE) {
    search_term <- paste0(species, "[Organism] AND chloroplast[Title] AND complete genome")
    search_results <- entrez_search(db = "nucleotide", term = search_term, retmax = 1)
  } else if (prot == TRUE) {
    search_term <- paste0(species, "[Organism] AND ", protein_name,"[Gene]")
    search_results <- entrez_search(db = "protein", term = search_term, retmax = 1)
  }
  
  if (length(search_results$ids) == 0) {
    warning(paste("No hit for", species))
    next
  }
  
  # Fetch the sequence in FASTA format
  if (prot == FALSE) {
    fasta_raw <- entrez_fetch(db = "nucleotide", id = search_results$ids[1], rettype = "fasta")
  } else {
    fasta_raw <- entrez_fetch(db = "protein", id = search_results$ids[1], rettype = "fasta")
  }
  
  # Extract just the sequence part
  lines <- unlist(strsplit(fasta_raw, "\n"))
  seq_lines <- lines[!grepl("^>", lines)]
  seq <- paste(seq_lines, collapse = "")
  
  # Add to list
  if (prot == FALSE) {
    sequence_list[[species_name]] <- DNAString(seq)
  } else {
    sequence_list[[species_name]] <- AAString(seq)
  }
}

# Convert to DNAStringSet or AAStringSet
if (prot == FALSE) {
  string_set <- DNAStringSet(sequence_list)
} else {
  string_set <- AAStringSet(sequence_list)
}

# Run alignment using msa
if (mito == TRUE) {
  alignment <- msaClustalW(string_set)
} else if (clp == TRUE) {
  alignment <- msa(string_set, method = "Muscle")
} else if (prot == TRUE) {
  alignment <- msaClustalW(string_set)
}

# Convert to character vector of aligned sequences
aligned_seqs <- as.character(alignment)

# Convert to list of character vectors (as required by seqinr::alignment)
aligned_list <- lapply(aligned_seqs, function(x) unlist(strsplit(x, split = "")))

# Create 'alignment' object for seqinr
seqinr_alignment <- list()
seqinr_alignment$nb <- length(aligned_list)
seqinr_alignment$nam <- names(aligned_list)
seqinr_alignment$seq <- lapply(aligned_list, function(x) paste(x, collapse = ""))
seqinr_alignment$alength <- nchar(seqinr_alignment$seq[[1]])
class(seqinr_alignment) <- "alignment"

idmx <- dist.alignment(seqinr_alignment,"identity")
mx <- 1 - as.matrix(idmx,nrow=n_species,ncol=n_species)


######################## Plots ###########################

message("Drawing plots...")

# Betweenness and withinness plots

### 1. Silhouette plot
jpeg("Silhouette.jpg")
fviz_nbclust(mx, kmeans, method = "silhouette", k.max=20) + theme_classic() # pam|kmeans|hcut|clara
dev.off()

### 2. Elbow plot
jpeg("Elbow.jpg")
fviz_nbclust(mx, kmeans, method = "wss", k.max=20) + theme_classic() # kmeans|pam (kmeans)
dev.off()

### 3. Heatmap
species <- row.names(mx)
myBreaks <- c(seq(0,1,by=0.01))
ceyy = cexx = 0.75

# normalize mx to 0-1
mx_hm <- mx
mx2 <- (mx_hm - min(mx_hm))/(max(mx_hm) - min(mx_hm))
mx_hm <- mx2

# color palette
clr = colorRampPalette(c("yellow","green","darkgreen"))(100) # plant studies
clr = colorRampPalette(c("green","white","red"))(100)
clr = colorRampPalette(c("white","yellow","red"))(100) # mitochondrial studies
clr = colorRampPalette(c("white","yellow","orange","red"))(100)
clr = gray.colors(100) # grayscale

# n is the estimated number of clusters/kinds/baramins
# if you change the estimate, restart from this point!
n <- 5

# clustering methods: ward.D ward.D2 single median average mcquitty complete centroid
clusmeth="ward.D2" 
row.clusters = hclust(dist(mx_hm),method=clusmeth)
col.clusters = hclust(dist(mx_hm),method=clusmeth)
ctk <- cutree(row.clusters,k=n)
clus_clrs <- rainbow(length(as.matrix(unique(ctk))))
clus_clrs_vec <- clus_clrs[ctk[species[sort(row.clusters$order)]]]

### Heatmap code ###
# provide a name for the heatmap
heatmap_name=paste("heatmap_",clusmeth,"_",n,"_250724.jpg",sep="")
jpeg(filename = heatmap_name, height = 3000, width = 3000, units = "px", res=400) # topo.colors(100) 5500, 5000
h <- heatmap(mx_hm, symkey=F, symbreaks=F, scale="none", dendrogram=F, Rowv=T, Colv=T, col = clr, RowSideColors=clus_clrs_vec, breaks = myBreaks, border_color=NA, na.color="white", margin = c(15,15),
             cexRow=cexx,cexCol=ceyy, key=T, trace="none", lmat=rbind( c(4, 3), c(2,1), c(0,0) ), lhei=c(1.8,6.5,1), hclustfun = function(d) hclust(d, method=clusmeth), # dendrogram="none",
             labCol=as.expression(lapply(colnames(mx), function(a) bquote(italic(.(a))))),labRow=as.expression(lapply(rownames(mx), function(a) bquote(italic(.(a))))))
invisible(dev.off())

######################### Stats files #########################

# stats and clusters files
# number of clusters
n <- 7

row.clusters = hclust(dist(mx),method=clusmeth)
ctk <- cutree(row.clusters,k=n)

# clusters from ctk
write.table(ctk, file="kclusters.txt", col.names=F, quote=F, sep="\t")

header = "cluster\tspecies\tmin\tmean\tmax\tSEM\tp-value\tneglog"
write(header, file="kstats.txt", sep="\t", append=T)

cluster_sizes = as.vector(table(ctk))
for (n_cluster in 1:n) {
  csize = cluster_sizes[n_cluster]
  if (csize >= 3) {
    m1 = as.matrix(mx[ctk == n_cluster,ctk == n_cluster])
    
    x = m1[upper.tri(m1)]
    ll = dim(m1)[1]
    
    m2 = as.matrix(cbind(mx[ctk != n_cluster,ctk == n_cluster],t(mx[ctk == n_cluster,ctk != n_cluster])))
    m2b = m2[!duplicated(colnames(m2))]
    
    t = t.test(x,m2b)
    pval = t$p.value
    nglog = -log10(pval)
    min = min(x)
    max = max(x)
    
    mean2 = sprintf("%.3f", mean(x))
    sem2 = sprintf("%.3f", sd(x)/sqrt(csize))
    min2 = sprintf("%.3f", min)
    max2 = sprintf("%.3f", max)
    pval2 = sprintf("%.3f", pval)
    nglog2 = sprintf("%.3f", nglog)
    
    stats = paste(n_cluster, ll, min2, mean2, max2, sem2, pval, nglog2, sep="\t")
    stats2 = gsub("\n","\t",stats)
    write(stats, file="kstats.txt", sep="\t", append = T)
  }
}

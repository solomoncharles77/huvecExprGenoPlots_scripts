# This script takes three arguments
# - Ensemble gene ID
# - rsID of SNP
# - SNP Coord
# This script outputs a hg38 version of the gwas summary statistics file 


# Set global options and load libraries -----------------------------------
options(scipen = 999, "warnPartialMatchDollar"=TRUE)
# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))


# List options
option_list = list(
  make_option(c("-g", "--geneID"), type="character", default=NULL, 
              help="Ensemble gene ID. eg ENSG00000050820", metavar="character"),
  make_option(c("-r", "--rsID"), type="character", default=NULL, 
              help="Rs ID of SNP. eg rs11648176", metavar="character"),
  make_option(c("-p", "--peakType"), type="character", default=NULL, 
              help="Specify what peak to plot (any of: ATACseq, CHIPseq)", metavar="character"),
  make_option(c("-c", "--Coord"), type="character", default = 5e-08, 
              help="pvalue threshold for selecting associations.", metavar="character")
); 


opt_parser <-  OptionParser(option_list=option_list);
opt <-  parse_args(opt_parser);

# Data input error messages

if (is.null(opt$geneID)){
  print_help(opt_parser)
  stop("Please provide gene ID.", call.=FALSE)
}

if (is.null(opt$rsID)){
  print_help(opt_parser)
  stop("Please provide rs ID.", call.=FALSE)
}

if (is.null(opt$peakType)){
  print_help(opt_parser)
  stop("Please provide peak type.", call.=FALSE)
}


# Import interval count --------------------------------------------------
# Fetch rsID details
peakID <- paste0(opt$peakType)
#peakID <- c("ATACseq")

if (peakID == "ATACseq") {
  interval <- data.frame(fread("/scratch/cellfunc/shared/HUVEC_ATACseq/eQTLs/eQTL_deseq2Norm_noMAF/norm_data/general/normalised_phenotype.tsv"))
  # Harmonize names and samples in geno and interval ----------------------------
  mapFile <- read.table("/scratch/cellfunc/shared/HUVEC_ATACseq/eQTLs/eQTL_deseq2Norm_noMAF/norm_data/general/mappings_handeling_repeats.tsv", header = T)
  names(interval) <- mapFile$Genotype[match(names(interval), mapFile$RNA)]
  colnames(interval) <- gsub(".*_S", "E", colnames(interval))
  names(interval)[1] <- c("intervalID")
  
  
} else if (peakID == "CHIPseq") {
  interval <- data.frame(fread("/scratch/cellfunc/shared/HUVEC_ChipSeq/eQTLs/eQTL_deseq2Norm_noMAF/norm_data/general/normalised_phenotype.tsv"))
  # Harmonize names and samples in geno and interval ----------------------------
  mapFile <- read.table("/scratch/cellfunc/shared/HUVEC_ChipSeq/eQTLs/eQTL_deseq2Norm_noMAF/norm_data/general/mappings_handeling_repeats.tsv", header = T)
  names(interval) <- mapFile$Genotype[match(names(interval), mapFile$RNA)]
  colnames(interval) <- gsub(".*_S", "E", colnames(interval))
  names(interval)[1] <- c("intervalID")
  
  
} else {
  cat("Please specify peak type")
}


# Extract interval for cand gene
cand <- paste0(opt$geneID)
#cand <- "Interval_22"
candInterval <- interval[interval$intervalID %in% cand, ]
candName <- candInterval$intervalID

# Fetch rsID details
rs <- paste0(opt$rsID)
#rs <- c("rs643058")
rsFile <- paste0("get_rsID_Coord.txt")
write.table(rs, file = rsFile, quote = F, row.names = F, col.names = F)
#system(paste0("grep -wFf ", rsFile, " /scratch/vasccell/cs806/colocalization/dbSNP/hg38.snp151_All.bed | cut -f4,3,10 > rsLookup.txt"))
rsCoord <- read.table("rsLookup.txt")
colnames(rsCoord) <- c("bp", "rsID", "coord")

# TODO
# I need some if commands around here.
#pos <- as.numeric(rsCoord$bp)
cd <- paste0("chr", rsCoord$coord)

#cd <- rsCoord$coord[2]

a1 <- unlist(strsplit(cd, split = ":"))[3]
a2 <- unlist(strsplit(cd, split = ":"))[4]
hom1 <- paste0(a1, a1)
het <- paste0(a1, a2)
hom2 <- paste0(a2, a2)


# Extract geno for cand SNP
exprLoc <- system(paste0("grep -n ", cd, " huvecGeno/huvec_Geno_snpsloc.txt"), intern = T)
toSkip <- as.numeric(gsub(":.*", "", exprLoc))
toSkip <- toSkip - 1
geno <- data.frame(fread("huvecGeno/huvec_Geno.txt", skip = toSkip, nrows = 1, col.names=scan('huvecGeno/huvecGenoHeader.txt', what ="", sep = "\t", quiet = TRUE)))


# Harmonize geno and expr
colnames(geno) <- gsub(".*_S", "E", colnames(geno))
genoSamps <- colnames(geno)[c(colnames(geno) %in% colnames(candInterval))]
candInterval <- candInterval[, genoSamps]

inteSamps <- colnames(candInterval)[c(colnames(candInterval) %in% colnames(geno))]
geno <- geno[, inteSamps]

candInterval <- data.frame(Gene = cand, candInterval)
geno <- data.frame(SNP = rs, geno)


df <- t(rbind(geno[-1], candInterval[-1]))
colnames(df) = c("snp", "gene_expr")
df <- data.frame(df)
df <- df[complete.cases(df), ]
df$snp = as.factor(df$snp)

dflm = lm(df[,"gene_expr"] ~ as.numeric(df[,"snp"]))

dfPlot <- ggplot(df, aes(snp, gene_expr)) +
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="lightgreen") + 
  labs(x = rs, y = candName)+
  scale_x_discrete(labels=c(hom1,het,hom2))+
  theme_bw()
# theme(axis.text.x = element_text(face="bold", color="#993333", 
#                                  size=14, angle=45))

ggsave(paste0(cand, "_vs_", rs, "_", peakID, "_Plot.png" ), dfPlot,
       scale = 1,
       path = paste0("exprGenoPlot/", peakID),
       width = 80,
       height = 70,
       units = "mm",
       dpi = 300,
       limitsize = TRUE)





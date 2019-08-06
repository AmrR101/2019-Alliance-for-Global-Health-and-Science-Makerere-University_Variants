
## Read in the annotated vcf file
l <- readLines("03-Freebayes/freebayes_annotated.vcf")

## Remove the descriptive portion (begins with #) and keep the variants
d <- l[-grep("^#", l)]
length(d) # check the length, make sure its the same number as reported snps

# extract the header line (this is the final descriptive line, before the variants start)
header <- strsplit(l[tail(grep("^#", l),1)],split="\t")[[1]] 
header

# use the header to determine the number of samples and get sample names
nsamples=length(header)-9
sample_names <- header[10:length(header)]

## parse the vcf table
d2 <- strsplit(d,split="\t")

# extract the SNP depth
DP = as.numeric(sub("DP=","",sapply(strsplit(sapply(d2, "[[", 8L),split=";"), function(x) grep("^DP=", x, value=T))))
summary(DP)

## Filters
## http://bcb.io/2014/05/12/wgs-trio-variant-evaluation/
hist(as.numeric(DP), breaks = 1000, main="Read depth histogram")
abline(v=(avg_depth + 4*sqrt(avg_depth)),col="red")
abline(v=(avg_depth - 4*sqrt(avg_depth)),col="red")
avg_depth <-  mean(DP,trim=0.2) 
keep_filter <- DP <= (avg_depth + 4*sqrt(avg_depth)) & DP >= (avg_depth - 4*sqrt(avg_depth))
table(keep_filter)

qual_cutoff <- avg_depth*2
qual <- sapply(d2, "[[", 6L)
summary(as.numeric(qual))
hist(as.numeric(qual), breaks = 1000, main="Quality histogram")
keep_filter  <- keep_filter & qual>=qual_cutoff 
table(keep_filter)

d2[keep_filter][[1]][9]
#[1] "GT:DP:AD:RO:QR:AO:QA:GL"

chrom <- sapply(d2[keep_filter], "[", 1L)
pos <- sapply(d2[keep_filter], "[[", 2L)
ref <- sapply(d2[keep_filter], "[[", 4L)
alt <- sapply(d2[keep_filter], "[[", 5L)
qual <- sapply(d2[keep_filter], "[[", 6L)
INFO <- sapply(d2[keep_filter], "[[", 8L)

DP = as.numeric(sub("DP=","",sapply(strsplit(sapply(d2[keep_filter], "[[", 8L),split=";"), function(x) grep("^DP=", x, value=T))))

GENO <- t(unname(sapply(d2[keep_filter], function(x) sapply(sapply(x[10:(10+nsamples-1)],strsplit,split=":"),"[[", 1L))))
colnames(GENO) <- sample_names
ord <- order(as.numeric(sub("sample","",colnames(GENO))))
GENO <- GENO[,ord]
class(GENO) <- "numeric" ## change to a numeric matrix
apply(GENO,2,function(x) table(is.na(x)))
GENO2 <- GENO[,-which(sample_names == "sample4")] ## sample 4 clearly bad lets remove

data_df <- data.frame(chrom,as.numeric(pos),ref,alt,as.numeric(qual),DP,INFO, stringsAsFactors = F)

# reference or uncalled
table(apply(GENO2,1, function(x) sum(x %in% c(NA,"0"))))              
boring_snps <- which(apply(GENO2,1, function(x) sum(x %in% c(NA,"0"))) == 14)
data_df <- data_df[-boring_snps,]
GENO2 <- GENO2[-boring_snps,]

table(sapply(strsplit(data_df$alt, split=","), length))
table(nchar(data_df$ref) == nchar(data_df$alt))
table(nchar(data_df$alt) == 1)

data_df$type = "SNP"
data_df$type[(nchar(data_df$ref) == nchar(data_df$alt) & nchar(data_df$ref) > 1  )] = "MNP"
data_df$type[(sapply(strsplit(data_df$alt, split=","), length) == 1 & nchar(data_df$ref) != nchar(data_df$alt))] = "INDEL"
data_df$type[sapply(strsplit(data_df$alt, split=","), length) > 1] = "COMPLEX"
table(data_df$type)

## Lets ignore complex

data_df2 <- data_df[data_df$type != "COMPLEX",]
GENO3 <- GENO2[data_df$type != "COMPLEX",]
dim(data_df2)
dim(GENO3)
table(GENO3,useNA = "always")

## Unieque variants per sample
bar1 <- colSums(sapply(seq.int(1,ncol(GENO3)),function(x) rowSums(GENO3[,-x],na.rm = T) == 0 & GENO3[,x] != 0))
bar1
bar2 <- colSums(sapply(seq.int(1,ncol(GENO3)),function(x) (rowSums(GENO3[,-x],na.rm = T) == 1) == ncol(GENO3) & GENO3[,x] == 0))
bar2
bar_gt <- data.frame(Sample=colnames(GENO3),"Unique SNPS"=bar1+bar2)
bar_gt

## UPGMA Tree
gt_d <- dist(t(GENO3),method = "manhattan")
plot(hclust(gt_d,method="average") , main="UPGMA tree from manhattan distances", xlab="",sub="",)

# Lets pick the ones that are different in all samples from ref
data_df_all <- data_df2[which(apply(GENO3, 1, function(x) all(x == 1))),]
GENO_all <- GENO3[which(apply(GENO3, 1, function(x) all(x == 1))),]

ANN_all <- sub("ANN=","",sapply(strsplit(data_df_all$INFO,split=";"), function(x) grep("^ANN=", x, value=T)))
strsplit(ANN_all, split=",")
strsplit(sapply(strsplit(ANN_all, split=","),"[[", 1L), split="|", fixed = T)

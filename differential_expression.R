library(tidyr)
library(tibble)
library("DEP")
library("dplyr")
library(SummarizedExperiment)
library(ggplot2)
library(collections)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("/Users/user/Documents/proteomics")

AYA<-read.csv("proteomics_gene_names.csv")
#correct fold change always GO over naive (it was mixed in original data)
A <- rowMeans(AYA[,c(23:25)]) # colnames "u87_naive1.1" "u87_naive2.1" "u87_naive3.1",
B <- rowMeans(AYA[,c(20:22)]) # "u87_Go100.1.1" "u87_Go100.2.1" "u87_Go100.3.1"
# fold change: log2 (B-A)/A A= initial value (naive) and B = final value (GO treated), but can just use means 
AYA$correctedfoldchange<- log(B/A, 2)

#Filter to remove the proteins with less than 3 unique peptides
AYA_filtered<-subset(AYA,Peptide.count>=3)
#Filter to use only fold changes >1.5 and p>0.05, need to do separately
# AYA_GO<- subset(AYA_filtered, correctedfoldchange>1.5) NO!
#AYA_GO<-subset(AYA_GO, Anova..p.>0.05) # 107 rows
drops <-c("Max.fold.change")
AYA_GO<-AYA_filtered[, !names(AYA_filtered) %in% drops]
# need to sort by p value

AYA_GO <- AYA_GO[order(AYA_GO$Anova..p., decreasing=TRUE),]

colnames(AYA_GO)[colnames(AYA_GO)=="correctedfoldchange"] <- "Max.fold.change"

# try using non normalised columns.... 
label <- c('u87_Go100.1.1', 'u87_Go100.2.1', 'u87_Go100.3.1', 'u87_naive1.1', 'u87_naive2.1', 'u87_naive3.1')
# conditions:
condition <- c('u87_Go100.1', 'u87_Go100.1', 'u87_Go100.1', 'u87_naive.1', 'u87_naive.1', 'u87_naive.1')
# replicate:
replicate <- c('1', '2', '3','1', '2', '3')

# assume that 1,2 and 3 are different time points, so should be analysed separately?

exp_design <- data.frame(label, condition, replicate)
AYA_GO_unique <- make_unique(AYA_GO, "Gene.name", "X.Accession", delim = ";")
u87_columns <- grep("u87", colnames(AYA_GO_unique))

# DEP functions throw errors about data format, these work for our data
make_se <- function(proteins_unique, columns, expdesign) {
raw <- proteins_unique[, columns]
raw[raw == 0] <- NA
raw <- log2(raw)
expdesign <- mutate(expdesign, condition = make.names(condition)) %>% unite(ID, condition, replicate, remove = FALSE)
rownames(expdesign) <- expdesign$ID
matched <- match(make.names(expdesign$label), make.names(colnames(raw)))
colnames(raw)[matched] <- expdesign$ID
raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
row_data <- proteins_unique[, -columns]
rownames(row_data) <- row_data$name
se <- SummarizedExperiment(assays = as.matrix(raw),colData = expdesign, rowData = row_data)
return(se)
}


filter_missval <- function(se, thr=0, max_repl=3){
thr <- 0
max_repl <-3
 bin_data <- assay(se)
  idx <- is.na(assay(se))
  bin_data[!idx] <- 1
  bin_data[idx] <- 0
keep <- bin_data %>%
    data.frame() %>%
    rownames_to_column() %>%
    gather(ID, value, -rowname) %>%
    left_join(., data.frame(colData(se)), by = "ID") %>%
    group_by(rowname, condition) %>%
    summarize(miss_val = n() - sum(value)) %>%
    filter(miss_val <= 0) %>%
    spread(condition, miss_val)
  se_fltrd <- se[keep$rowname, ]
  return(se_fltrd)
}


data_se <-make_se(AYA_GO_unique, u87_columns[7:12], exp_design)
data_filt <-filter_missval(data_se)
data_norm <- normalize_vsn(data_filt)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_diff <- test_diff(data_imp, type = "control", control = "u87_naive.1")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

plot_volcano(dep, contrast = "u87_Go100.1_vs_u87_naive.1", label_size = 2, add_names = TRUE)
# above not working.... so  

contrast <- "u87_Go100.1_vs_u87_naive.1" 
label_size <- 2
add_names <- TRUE
adjusted <- FALSE
row_data <- rowData(dep, use.names = FALSE)

plot_volcano <- function(contrast, label_size, add_names, adjusted){
 diff <- grep(paste(contrast, "_diff", sep = ""),
               colnames(row_data))
  if(adjusted) {
    p_values <- grep(paste(contrast, "_p.adj", sep = ""),
                     colnames(row_data))
  } else {
    p_values <- grep(paste(contrast, "_p.val", sep = ""),
                     colnames(row_data))
  }
  signif <- grep(paste(contrast, "_significant", sep = ""),
                 colnames(row_data))
# this is where it fails
  df <- data.frame(x = row_data[, diff],
                   y = -log10(row_data[, p_values]),
                   significant = row_data[, signif],
                   name = row_data$name) %>%
    filter(!is.na(significant)) %>%
    arrange(significant)
name1 <- gsub("_vs_.*", "", contrast)
name2 <- gsub(".*_vs_", "", contrast)
p <- ggplot(df, aes(x, y)) +
    geom_vline(xintercept = 0) +
    geom_point(aes(col = significant)) +
    geom_text(data = data.frame(), aes(x = c(Inf, -Inf),
                                       y = c(-Inf, -Inf),
                                       hjust = c(1, 0),
                                       vjust = c(-1, -1),
                                       label = c(name1, name2),
                                       size = 5,
                                       fontface = "bold")) +
    labs(title = contrast,
      x = expression(log[2]~"Fold change")) +
    theme_DEP1() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey"))
  if (add_names) {
    p <- p + ggrepel::geom_text_repel(data = filter(df, significant),
                                      aes(label = name),
                                      size = label_size,
                                      box.padding = unit(0.1, 'lines'),
                                      point.padding = unit(0.1, 'lines'),
                                      segment.size = 0.5)
  }
  if(adjusted) {
    p <- p + labs(y = expression(-log[10]~"Adjusted p-value"))
  } else {
    p <- p + labs(y = expression(-log[10]~"P-value"))
  }
  return(p)
}

data_results <- get_results(dep)

sig <- data_results %>% filter(significant) # 37 significant
prots <- sig$ID

genes <-mapIds(org.Hs.eg.db, as.character(prots), "ENTREZID", "UNIPROT")

gene_id_dict = Dict$new(genes)
 #populate list of ensemble names to pass into enrichGO
 l_entrez = c()
 for (i in prots){
    l_entrez = c(l_entrez, gene_id_dict$get(i))
    }

go_cc <- enrichGO(gene          = l_entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

p <- dotplot(go_cc, orderBy="BgRatio")

go_mf <- enrichGO(gene          = l_entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

go_bp <- enrichGO(gene          = l_entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

p_mf <- dotplot(go_mf, orderBy="GeneRatio")
p_bp <- dotplot(go_bp, orderBy="GeneRatio")
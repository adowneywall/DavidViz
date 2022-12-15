
library(tidyverse)
setwd("Path/To/DAVIDKnowledgebaseFiles")

## BP ##
g_BP <- read_tsv('OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt',col_names = F)
colnames(g_BP) <- c('Gene','Term')
g_BP_tab <- g_BP %>%
  group_by(Term) %>%
  summarize(Term_Type = 'GOTERM_FAT_BP',
            Gene_count = n_distinct(Gene),
            Gene_list = list(Gene))
saveRDS(g_BP_tab,'DAVID_BP_FAT_GOTERMS.RData')

length(unique(unlist(g_BP_tab$Gene_list)))
## MF ##
g_term <- read_tsv('OFFICIAL_GENE_SYMBOL2GOTERM_MF_FAT.txt',col_names = F)
colnames(g_term ) <- c('Gene','Term')
g_term_tab <- g_term  %>%
  group_by(Term) %>%
  summarize(Term_Type = 'GOTERM_FAT_MF',,
            Gene_count = n_distinct(Gene),
            Gene_list = list(Gene))
saveRDS(g_term_tab,'DAVID_MF_FAT_GOTERMS.RData')

## CC ##
g_term <- read_tsv('OFFICIAL_GENE_SYMBOL2GOTERM_CC_FAT.txt',col_names = F)
colnames(g_term ) <- c('Gene','Term')
g_term_tab <- g_term  %>%
  group_by(Term) %>%
  summarize(Term_Type = 'GOTERM_FAT_CC',
            Gene_count = n_distinct(Gene),
            Gene_list = list(Gene))
saveRDS(g_term_tab,'DAVID_CC_FAT_GOTERMS.RData')

## UP_SEQ ##
g_term <- read_tsv('OFFICIAL_GENE_SYMBOL2UP_SEQ_FEATURE.txt',col_names = F)
colnames(g_term ) <- c('Gene','Term')
g_term_tab <- g_term  %>%
  group_by(Term) %>%
  summarize(Term_Type = 'SEQ_UP',
            Gene_count = n_distinct(Gene),
            Gene_list = list(Gene))
saveRDS(g_term_tab,'DAVID_UP_SEQ_GOTERMS.RData')
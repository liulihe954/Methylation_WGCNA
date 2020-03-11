load('gene_symbols_genome.rda')
head(gene_symbols_genome)


library(biomaRt)
biomart="ensembl"
dataset="btaurus_gene_ensembl"
attributes = c("ensembl_gene_id",
               "hgnc_symbol",
               "external_gene_name",
               "entrezgene_id")
database = useMart(biomart)
genome = useDataset(dataset, mart = database)
gene = getBM(attributes,mart = genome)

#searchAttributes(genome,'external_gene_name')

head(gene,300)

checkGeneSymbols('HHAT')
#alias2Symbol('LOC786897',species = '',expand.symbols = F)

table(is.na(gene_symbols_genome))

test = gene_symbols_genome %>% 
  left_join(gene,by = c('ensembl_gene_id' = 'ensembl_gene_id'))

out = c()
for (i in seq_len(dim(test)[1])){
  v1 = test[i,3]
  v2 = test[i,6]
  out[i] =  ifelse(v1 == v2,'1','0')
}

table(out)


fal = which(test$ENTREZID != test$entrezgene_id)

test[fal,c(1,3,6)]
names(test)

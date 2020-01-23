library(tidyverse)
# specify database location
url_template = 'http://software.broadinstitute.org/gsea/msigdb/cards/'
Msig_db_destination = '/ufrc/penagaricano/lihe.liu/Methylation_WGCNA/Msig_db'
# get pathway name index
library(msigdbr)
m_df = msigdbr(species = "Bos taurus")
#head(unique(m_df$gs_name),40)

# obtain name index and paste to all urls
Msig_name_index = unique(m_df$gs_name)
Msig_urls = paste(url_template,Msig_name_index,sep = "")
# prepare R function to retrive description
Get_Descrip = function(URL,
                       selector = "td",
                       pos = 6){
  library(tidyverse)
  library(rvest)
  base_url <- URL
  webpage <- read_html(base_url)
  # Get the artist name
  target_raw <- html_nodes(webpage,selector)
  target_raw <- as.character(html_text( target_raw))
  text = str_replace_all( target_raw,"[\r\n]","")
  final = as.character(text[pos])
  message('try on ', URL)
  return(final)
}

setwd(Msig_db_destination)
load('Msigdb_bta.RData')
# # loop to retrive
# #Msig_urls_t = Msig_urls[1:20]
# All_Descrip = c()
# for (i in seq_along(Msig_urls)){
#   #message(i)
#   All_Descrip[i] = Get_Descrip(Msig_urls[i])
# }
# 
# # massage
# #Msig_name_index_t = Msig_name_index[1:20]
# Universe_Descrip = data.frame(cbind(Msig_name_index,All_Descrip))
# colnames(Universe_Descrip) = c('gs_name','gd_description')
# 
# 
# m_df_all = dplyr::left_join(m_df,
#                             Universe_Descrip,
#                             by = c("gs_name" = "gs_name"))
# 

DB_List = list()
MsigRecords = unique(m_df_all[,c("gs_name","gd_description")]) %>% arrange(gs_name)
MsigID = na.omit(MsigRecords$gs_name)
MsigTerm = na.omit(MsigRecords$gd_description)
for ( p in seq_along(MsigID)){
  tmp = subset(m_df_all, gs_name == MsigID[p])$entrez_gene
  DB_List[[p]] = tmp #
  names(DB_List)[p]  <- paste(MsigID[p],"-")#MeshTerm[p]
  print(p)
}


setwd(Msig_db_destination)
save(m_df_all,
     DB_List,
     file = "Msigdb_bta.RData")

# load('Msigdb_bta.RData')

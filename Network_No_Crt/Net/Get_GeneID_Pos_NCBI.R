
Get_Pos_NCBIweb = function(GeneID,
                           selector = "td",
                           pos=5){
  
  NCBI_template = 'https://www.ncbi.nlm.nih.gov/gene/?term='
  base_url = paste(NCBI_template,GeneID,sep = '')
  webpage <- read_html(base_url)
  target_raw <- html_nodes(webpage,selector)
  target_raw <- as.character(html_text( target_raw))
  text = str_replace_all( target_raw,"[\r\n]","")
  string = as.character(text[pos])
  temp0 =  regexpr(' ', string)
  ID0 = substr(string,1,temp0[1]-1)
  temp = regexpr('\\(', string)
  string2 = substr(string, temp[1], nchar(string))
  temp2 = regexpr('\\.', string2)
  start_pos = as.numeric(substr(string2, 2, temp2[1]-1))
  end_pos = as.numeric(str_replace_all(substr(string2, temp2[1]+2, nchar(string2)),
                                       "\\D+",""))
  message('try on ', GeneID)
  return(list(ID = ID0,
              start = start_pos,
              end = end_pos))
}

GeneID = 'ENSBTAG00000000005'


test_list = c('ENSBTAG00000000005',
              'ENSBTAG00000000008',
              'ENSBTAG00000000009')

Fetch_Pos = Get_Pos_NCBIweb(test_list[1])



GeneBank_template = 'https://www.ncbi.nlm.nih.gov/nuccore/'
#ID = 'NC_037347.1'
Linker1 = '?report=genbank&from='
#Pos1 = 33708356
Linker2 = '&to='
#Pos2 = 33737405

GeneBack_target = paste(GeneBank_template,Fetch_Pos$ID,
                        Linker1,Fetch_Pos$start,Linker2,Fetch_Pos$end,
                        sep = '')

#viewercontent1 > div > div > pre
webpage <- read_html(GeneBack_target)

target_raw =
  webpage %>% html_nodes("body") %>%
  html_nodes("#viewercontent1") %>% 
  magrittr::extract2(1)

  
  magrittr::extract2(1) %>%
  html_nodes("form") %>%
  html_nodes("div") %>% magrittr::extract2(1) %>%
  html_nodes("div") %>% magrittr::extract2(4) %>%
  html_nodes("div") %>%
  html_nodes("div") %>% magrittr::extract2(5) %>%
  html_nodes("div") %>% magrittr::extract2(2) %>%
  html_nodes("div") %>% magrittr::extract2(1) %>%
  html_nodes("div") %>%
  html_nodes("div") %>%
  html_nodes("pre") %>%
  html_nodes("span")%>% magrittr::extract2(14)

/html/body/div/div[1]/form/div[1]/div[4]/div/div[5]/div[2]/div[1]/div/div/pre/span[14]
# webpage %>% 
#   html_nodes(xpath = '//*[@id="NC_037344.1_1"]') %>%
#   xml_attr("value")

target_raw <- html_nodes(webpage,"body")

#target_raw <- html_nodes(webpage,xpath = '/html/body/div/div[1]/form/div[1]/div[4]/div/div[5]/div[2]/div[1]/div/div/pre')
target_raw <- as.character(html_text(target_raw))
text = str_replace_all(target_raw,"[\r\n]","")
length(text)

head(text)

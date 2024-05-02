#----- Human gene editing -----#
#---------1. Analyses ---------#

#library
library('dplyr')
library('tidyr')
library('tidyverse')
library('stringr')
library('tidytext')
library('textdata')
library('xlsx')

#Data
load("Gene_NP.rda")

#Outlet
Outlet <- Gene_NP %>% select(id, outlet) %>%
  count(outlet, sort=TRUE) %>%
  mutate(outlet = reorder(outlet, n))

#write.xlsx(Outlet, file = "Gene_Outlet.xlsx", 
#           sheetName="Combined", append=FALSE)
#write.xlsx(Outlet_Nexis, file = "Gene_Outlet.xlsx",
#           sheetName ="Nexis", append = TRUE)
#write.xlsx(Outlet_ProQuest, file = "Gene_Outlet.xlsx", 
#           sheetName="ProQuest", append=TRUE)

#Year
Year <- Gene_NP%>%select(id, year) %>%
  count(year, sort=TRUE) %>%
  mutate(year = reorder(year, n)) %>%
  filter(n <= 30)

#write.xlsx(Year, file = "Gene_Year.xlsx", 
#           sheetName="Combined", append=FALSE)


#Sentimental analysis

#news data
Gene_NP$article <- tolower(Gene_NP$article)
article <- Gene_NP%>%select(id, article)

article_preprocessed<-article%>%
  mutate(article=str_replace_all(article,"https://t.co/[A-Za-z\\d]+|http://[A-Za-z\\d]+|&amp;|&lt;|&gt;|RT|https", ""))%>% #remove https
  unnest_tokens(word, article, token = "regex", pattern = "([^A-Za-z_\\d#@']|'(?![A-Za-z_\\d#@]))")%>%  #word tokenization
  filter(!word %in% stop_words$word,
         str_detect(word, "[a-z]")) # remove stopwords

#sentiment analysis - nrc
get_sentiments('nrc')

trust <- get_sentiments("nrc") %>% 
  filter(sentiment == "trust")

gene_trust <- article_preprocessed %>%
  inner_join(trust) %>%
  count(word, sort = TRUE)%>%
  ungroup()

gene_trust
sum(gene_trust$n)

#sentiment by document
getCount <- function(data, keyword)
{
  wcount <- str_count(Gene_NP$article, keyword)
  return(data.frame(data, wcount))
}
getCount(Gene_NP, 'trust')

###################################

# Select WSJ, NYT, USA Today
Gene_NP_3 <- Gene_NP%>% filter(str_detect(outlet, "Wall Street|New York Times|USA Today"))
save(Gene_NP_3, file = "C:\\Users\\Ahn Suhwoo\\Desktop\\Gene\\Analyses\\Gene_NP_3.rda")

library(xtable)
print(xtable(Gene_NP_3), type = "html")

write.xlsx(Gene_NP_3, file = "Gene_NP_3.xlsx", 
           sheetName="WSJ_NYT_Today", append=FALSE)

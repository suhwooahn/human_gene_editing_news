#----- Human gene editing -----#
#------ 2. Topic modeling -----#

library(tidyverse)
library(tidytext)
library(stringr)

load("Gene_NP.rda")

# Remove press release
Gene_NP <- Gene_NP[grepl("PR Newswire|ReleaseWire|Marketwired|iCrowdNewswire (English)|PR.com|ews Release Wire|Congressional Press Releases", Gene_NP$outlet),]

# Select legacy media
Gene_NP <- Gene_NP%>% filter(str_detect(outlet, "Wall Street|New York Times|NPR|Business Insider|CNN|Reuters|Los Angeles Times|Forbes|USA Today|NBC|CBS|Bloomberg|ABC|Fox|MSNBC|Politico"))
Gene_NP <- Gene_NP%>% filter(str_detect(outlet, "Wall Street|New York Times|NPR|CNN|Reuters|Los Angeles Times|USA Today|NBC|CBS|ABC|Fox|MSNBC|Politico"))

#------------------------------------#
#------Topic modeling with STM-------#
#------------------------------------#

#https://github.com/dondealban/learning-stm
#https://blog.naver.com/PostView.naver?blogId=statstorm&logNo=222136015862&categoryNo=72&parentCategoryNo=0

Gene_NP$article <- str_replace(Gene_NP$article,"https://t.co/[A-Za-z\\d]+|http://[A-Za-z\\d]+|&amp;|&lt;|&gt;|RT|https", "")
Gene_NP$article <- str_replace(Gene_NP$article,".com", "")

library(stm)

Processed <- textProcessor(Gene_NP$article, metadata=Gene_NP)
Processed <- textProcessor(Gene_NP$title, metadata=Gene_NP)

out <- prepDocuments(Processed$documents, Processed$vocab, Processed$meta)
docs <- out$documents
vocab <- out$vocab
meta <-out$meta

par(mar = c(1, 1, 1, 1))
plotRemoved(Processed$documents, lower.thresh=seq(1,200, by=100))

poliblogPrevFit <- stm(out$documents, out$vocab, K=4, #prevalence=~rating+s(day),
                       max.em.its=75, data=out$meta, init.type="Spectral", 
                       seed=8458159)

par(mfrow = c(1,1))

plot(poliblogPrevFit, type="summary", xlim=c(0,.4))

plot(poliblogPrevFit, type="labels", topics=c(1,2,3))

plot(poliblogPrevFit, type="perspectives", topics=c(1,2))

#------------------------------------#
#------Topic modeling with LDA-------#
#------------------------------------#

# Pre-processing textual data
article_preprocessed <- Gene_NP %>%
  mutate(headline = str_replace(article,"https://t.co/[A-Za-z\\d]+|http://[A-Za-z\\d]+|&amp;|&lt;|&gt;|RT|https", "")) %>% #remove https
  unnest_tokens(word, article, token = "regex", pattern = "([^A-Za-z_\\d#@']|'(?![A-Za-z_\\d#@]))") %>%  #word tokenization
  filter(!word %in% stop_words$word,
         str_detect(word, "[a-z]")) # remove stopwords

## remove customized stop words
new_stopwords<-c("www","http","pm","am")
remove.list <- paste(new_stopwords, collapse = '|')
article_preprocessed<-article_preprocessed %>% filter(!grepl(remove.list, word))

## lemmatization
library(textstem)
article_preprocessed$word=lemmatize_words(article_preprocessed$word, dictionary = lexicon::hash_lemmas)
article_preprocessed

# tf-idf
id_words<-article_preprocessed %>%
  count(id,word,sort=TRUE)

total_words<-id_words %>% group_by(id) %>% summarize(total = sum(n))
id_words<-left_join(id_words, total_words)

id_words<-id_words %>%
  bind_tf_idf(word,id,n)

summary(id_words$tf_idf)
length(unique(id_words$id))

# DTM
## each row represents one document 
## each column represents one term and 
## each value contains the number of appearances of that term in that document
library(tm)

?cast_dtm ### cast()=> turn a tidy one-term-per-row data frame into a matrix; cast_dtm(): convert to a DocumnetTermMatrix            
dtm<-id_words%>%
  count(id, word)%>%
  cast_dtm(id,word,n)

dtm # weighting by term-frequency
Terms<-Terms(dtm)
head(Terms)

#dtm

#To remove sparse terms in the dtm
dtm_removesparseterms <- removeSparseTerms(dtm, 0.99)
dtm_removesparseterms
raw.sum <- apply(dtm_removesparseterms,1,FUN=sum) #sum by raw each raw of the table
dtm_removesparseterms <- dtm_removesparseterms[raw.sum!=0,]
dtm_removesparseterms

#dtm_tfidf <- weightTfIdf(dtm)
#dtm_tfidf
#dtm_tfidf <- removeSparseTerms(dtm_tfidf, 0.95)
#dtm_tfidf

# number of topics
K <- 3

# set random number generator seed
set.seed(9161)

# compute the LDA model, inference via 1000 iterations of Gibbs sampling
library("topicmodels")
topicModel <- LDA(dtm_removesparseterms, K, method="Gibbs", control=list(iter = 500, verbose = 25))

# have a look a some of the results (posterior distributions)
tmResult <- posterior(topicModel)
# format of the resulting object
attributes(tmResult)

# topics are probability distributions over the entire vocabulary
beta <- tmResult$terms   # get beta from results
dim(beta)                # K distributions over nTerms(DTM) terms

# for every document we have a probability distribution of its contained topics
theta <- tmResult$topics 
dim(theta)               # nDocs(DTM) distributions over K topics

doc_topic_prob <- as_data_frame(theta)
doc_topic_prob$id <- row.names(theta)

doc_topic_prob <- doc_topic_prob %>%
  gather(topic, prob, -c(id)) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(rank = rank(-prob, ties.method = "random")) %>%
  ungroup() %>%
  arrange(id, -prob)

top_topic <- doc_topic_prob %>%
  filter(rank == 1)

top_topic %>%
  group_by(topic) %>%
  summarize(freq = length(id)) %>%
  ungroup() %>%
  mutate(perc = freq/sum(freq))

# Check the top 10 most likely tersm with the term probabilities beta of the inferred topics
terms(topicModel, 20)

### Term-Topic matrix ###
beta
trans_beta <- t(beta)

library('xlsx')
write.xlsx(trans_beta, file = "gene_topic_terms.xlsx")

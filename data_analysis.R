library(mixOmics)
library(MetabFUN)
library(readxl)
library(stringr)

if (dir.exists('C:/Users/dbraas/Dropbox/projects')==T) dropbox <- 'C:/Users/dbraas/Dropbox/projects/'
if (dir.exists('C:/Users/Danie/Dropbox/projects')==T) dropbox <- 'C:/Users/Danie/Dropbox/projects/'
if (dir.exists('C:/Users/Daniel/Dropbox/projects')==T) dropbox <- 'C:/Users/Daniel/Dropbox/projects/'
if (dir.exists('D:/Dropbox/projects')==T) dropbox <- 'D:/Dropbox/projects/'
Abbrev <- read.csv(paste0(dropbox,'Metabolomics/Abbreviations.csv'), header = T)

wd <- paste0(dropbox,'Christina Charles/Maven/')
setwd(wd)

info <- read_excel("Metabolomics List- Braas 5-9-17.xlsx", sheet = 1)
info$Sample <- paste('AS',str_pad(1:96, width=3, pad='0'),sep='.')

neg_data <- data.frame(read.table("negpeakstringent110717.csv", header=T, sep=','))
neg_data$Polarity <- 'neg'

#removal of isotopologues
neg_data <- arrange(neg_data, round(medRt,1), medMz)

Rt_diff <- as.numeric()
Mz_diff <- as.numeric()
for (i in 1:nrow(neg_data)){
  if (i == 1) {
    Rt_diff <- append(Rt_diff, round(neg_data$medRt[1],1))
    Mz_diff <- append(Mz_diff, neg_data$medMz[1])
  } else {
    Rt_diff <- append(Rt_diff, round(neg_data$medRt[i]-neg_data$medRt[i-1]),1)
    Mz_diff <- append(Mz_diff, abs(neg_data$medMz[i]-neg_data$medMz[i-1]))
  }
}
neg_data$Rt_diff <- Rt_diff
neg_data$Mz_diff <- Mz_diff

delete <- as.numeric()
for (i in 2:nrow(neg_data)){
  if ((neg_data$Rt_diff[i] <= 0.2 & 
       1.005 > neg_data$Mz_diff[i] &
       neg_data$Mz_diff[i] > 1.002) == T)
    delete <- append(delete, i)
}

neg_data <- neg_data[-delete,]

neg_data <- arrange(neg_data, medMz, round(medRt,1))

Rt_diff <- as.numeric()
Mz_diff <- as.numeric()
for (i in 1:nrow(neg_data)){
  if (i == 1) {
    Rt_diff <- append(Rt_diff, round(neg_data$medRt[1]),1)
    Mz_diff <- append(Mz_diff, neg_data$medMz[1])
  } else {
    Rt_diff <- append(Rt_diff, round(neg_data$medRt[i]-neg_data$medRt[i-1],1))
    Mz_diff <- append(Mz_diff, abs(neg_data$medMz[i]-neg_data$medMz[i-1]))
  }
}
neg_data$Rt_diff <- Rt_diff
neg_data$Mz_diff <- Mz_diff

delete <- as.numeric()
for (i in 2:nrow(neg_data)){
  if ((neg_data$Rt_diff[i] <= 0.2 & 
      1.005 > neg_data$Mz_diff[i] &
      neg_data$Mz_diff[i] > 1.002) == T)
    delete <- append(delete, i)
}

neg_data <- neg_data[-delete,]

neg_data <- neg_data %>% 
  select(groupId, medMz, medRt, maxQuality, Polarity, AS.001:AS.096) %>% 
  gather(Sample, Area, -groupId, -medMz, -medRt, -maxQuality, -Polarity) %>% 
  left_join(., info, by='Sample')

neg <- select(neg_data, medMz, medRt, Sample, Area) %>% 
  unite(Mz_Rt, c(medMz, medRt), sep='_') %>% 
  spread(Mz_Rt, Area)

neg_mat <- as.matrix(neg[,2:length(neg)])
rownames(neg_mat) <- neg$Sample

neg_pca <- prcomp(neg_mat, scale=T, center=T)
neg_PC <- neg_pca$x %>% 
  as.data.frame() %>% 
  mutate(Sample=rownames(.)) %>% 
  left_join(., info, by='Sample')

ggplot(neg_PC, aes(PC1, PC2, fill=Condition))+
  geom_point(shape=21, size=3)

neg_PLSDA <- plsda(X = neg_mat,
                   Y = info$Condition,
                   scale=T)

a <- neg_PLSDA$variates$X %>% 
  data.frame() %>% 
  mutate(Sample = rownames(.)) %>% 
  left_join(., info, by='Sample') %>% 
  ggplot(., aes(comp.1, comp.2, fill=Condition, label=Sample))+
  geom_point(shape=21, size=5)+
  geom_text(vjust=-1, size=3)+
  scale_fill_brewer('Condition', palette = 'Set1')+
  theme_bw()+
  theme(text=element_text(face = 'bold'))+
  labs(title = 'PLSDA model from negative ions only',
       x=paste0('Comp 1 (', round(neg_PLSDA$explained_variance$X[1]*100,1),'%)'),
       y=paste0('Comp 2 (', round(neg_PLSDA$explained_variance$X[2]*100,1),'%)'))
  
neg_CCP <- cor(neg_PLSDA$X, neg_PLSDA$variates$X, use='pairwise') %>% 
  data.frame() %>% 
  mutate(Ion = rownames(.),
         Corr = sqrt(comp.1^2+comp.2^2))  

b <- filter(neg_CCP, Corr >= 0.5) %>% 
  ggplot(., aes(comp.1,comp.2, label=Ion))+
  geom_point()+
  geom_text(vjust=-1, size=3, color='navy')+
  geom_hline(yintercept = 0, lty=3)+
  geom_vline(xintercept = 0, lty=3)+
  geom_point(aes(x=0,y=0), color='darkgreen')+
  xlim(-1,1)+
  ylim(-1,1)+
  theme_bw()+
  theme(text=element_text(face = 'bold'))+
  labs(title = 'PLSDA model from negative ions only',
       x=paste0('Comp 1 (', round(neg_PLSDA$explained_variance$X[1]*100,1),'%)'),
       y=paste0('Comp 2 (', round(neg_PLSDA$explained_variance$X[2]*100,1),'%)'))

gridExtra::grid.arrange(a,b, nrow=1)  


# doing the analysis with positive mode ions ------------------------------

pos_data <- read.table("stringentpos110817.csv", header=T, sep=',')
pos_data$Polarity <- 'pos'#removal of isotopologues
pos_data <- arrange(pos_data, round(medRt,1), medMz)

Rt_diff <- as.numeric()
Mz_diff <- as.numeric()
for (i in 1:nrow(pos_data)){
  if (i == 1) {
    Rt_diff <- append(Rt_diff, round(pos_data$medRt[1],1))
    Mz_diff <- append(Mz_diff, pos_data$medMz[1])
  } else {
    Rt_diff <- append(Rt_diff, round(pos_data$medRt[i]-pos_data$medRt[i-1]),1)
    Mz_diff <- append(Mz_diff, abs(pos_data$medMz[i]-pos_data$medMz[i-1]))
  }
}
pos_data$Rt_diff <- Rt_diff
pos_data$Mz_diff <- Mz_diff

delete <- as.numeric()
for (i in 2:nrow(pos_data)){
  if ((pos_data$Rt_diff[i] <= 0.2 & 
       1.005 > pos_data$Mz_diff[i] &
       pos_data$Mz_diff[i] > 1.002) == T)
    delete <- append(delete, i)
}

pos_data <- pos_data[-delete,]

pos_data <- arrange(pos_data, medMz, round(medRt,1))

Rt_diff <- as.numeric()
Mz_diff <- as.numeric()
for (i in 1:nrow(pos_data)){
  if (i == 1) {
    Rt_diff <- append(Rt_diff, round(pos_data$medRt[1]),1)
    Mz_diff <- append(Mz_diff, pos_data$medMz[1])
  } else {
    Rt_diff <- append(Rt_diff, round(pos_data$medRt[i]-pos_data$medRt[i-1],1))
    Mz_diff <- append(Mz_diff, abs(pos_data$medMz[i]-pos_data$medMz[i-1]))
  }
}
pos_data$Rt_diff <- Rt_diff
pos_data$Mz_diff <- Mz_diff

delete <- as.numeric()
for (i in 2:nrow(pos_data)){
  if ((pos_data$Rt_diff[i] <= 0.2 & 
       1.005 > pos_data$Mz_diff[i] &
       pos_data$Mz_diff[i] > 1.002) == T)
    delete <- append(delete, i)
}

pos_data <- pos_data[-delete,]

pos_data <- pos_data %>% 
  select(groupId, medMz, medRt, maxQuality, Polarity, AS.001:AS.096) %>% 
  gather(Sample, Area, -groupId, -medMz, -medRt, -maxQuality, -Polarity) %>% 
  left_join(., info, by='Sample')

pos <- select(pos_data, medMz, medRt, Sample, Area) %>% 
  unite(Mz_Rt, c(medMz, medRt), sep='_') %>% 
  spread(Mz_Rt, Area)

pos_mat <- as.matrix(pos[,2:length(pos)])
rownames(pos_mat) <- pos$Sample

pos_pca <- prcomp(pos_mat, scale=T, center=T)
pos_PC <- pos_pca$x %>% 
  as.data.frame() %>% 
  mutate(Sample=rownames(.)) %>% 
  left_join(., info, by='Sample')

ggplot(pos_PC, aes(PC1, PC2, fill=Condition))+
  geom_point(shape=21, size=3)

pos_PLSDA <- plsda(X = pos_mat,
                   Y = info$Condition,
                   scale=T)

c <- pos_PLSDA$variates$X %>% 
  data.frame() %>% 
  mutate(Sample = rownames(.)) %>% 
  left_join(., info, by='Sample') %>% 
  ggplot(., aes(comp.1, comp.2, fill=Condition, label=Sample))+
  geom_point(shape=21, size=5)+
  geom_text(vjust=-1, size=3)+
  scale_fill_brewer('Condition', palette = 'Set1')+
  theme_bw()+
  theme(text=element_text(face = 'bold'))+
  labs(title = 'PLSDA model from positive ions only',
       x=paste0('Comp 1 (', round(pos_PLSDA$explained_variance$X[1]*100,1),'%)'),
       y=paste0('Comp 2 (', round(pos_PLSDA$explained_variance$X[2]*100,1),'%)'))

pos_CCP <- cor(pos_PLSDA$X, pos_PLSDA$variates$X, use='pairwise') %>% 
  data.frame() %>% 
  mutate(Ion = rownames(.),
         Corr = sqrt(comp.1^2+comp.2^2))  

d <- filter(pos_CCP, Corr >= 0.5) %>% 
  ggplot(., aes(comp.1,comp.2, label=Ion))+
  geom_point()+
  geom_text(vjust=-1, size=3, color='navy')+
  geom_hline(yintercept = 0, lty=3)+
  geom_vline(xintercept = 0, lty=3)+
  geom_point(aes(x=0,y=0), color='darkgreen')+
  xlim(-1,1)+
  ylim(-1,1)+
  theme_bw()+
  theme(text=element_text(face = 'bold'))+
  labs(title = 'PLSDA model from positive ions only',
       x=paste0('Comp 1 (', round(pos_PLSDA$explained_variance$X[1]*100,1),'%)'),
       y=paste0('Comp 2 (', round(pos_PLSDA$explained_variance$X[2]*100,1),'%)'))

gridExtra::grid.arrange(c,d, nrow=1)  

# data analysis with both positive and negative ions ----------------------

all_data <- bind_rows(neg_data, pos_data)
all_data$Polarity <- ifelse(all_data$Polarity=='neg', '-','+')

all <- select(all_data, Polarity, medMz, medRt, Sample, Area) %>% 
  unite(Pol_Mz_Rt, c(Polarity,medMz, medRt), sep='_') %>% 
  spread(Pol_Mz_Rt, Area)

all_mat <- as.matrix(all[,2:length(all)])
rownames(all_mat) <- all$Sample

all_pca <- prcomp(all_mat, scale=T, center=T)
all_PC <- all_pca$x %>% 
  as.data.frame() %>% 
  mutate(Sample=rownames(.)) %>% 
  left_join(., info, by='Sample')

ggplot(all_PC, aes(PC1, PC2, fill=Condition))+
  geom_point(shape=21, size=3)

all_PLSDA <- plsda(X = all_mat,
                   Y = info$Condition,
                   scale=T)

e <- all_PLSDA$variates$X %>% 
  data.frame() %>% 
  mutate(Sample = rownames(.)) %>% 
  left_join(., info, by='Sample') %>% 
  ggplot(., aes(comp.1, comp.2, fill=Condition, label=Sample))+
  geom_point(shape=21, size=5)+
  geom_text(vjust=-1, size=3)+
  scale_fill_brewer('Condition', palette = 'Set1')+
  theme_bw()+
  theme(text=element_text(face = 'bold'))+
  labs(title = 'PLSDA model from all ions',
       x=paste0('Comp 1 (', round(neg_PLSDA$explained_variance$X[1]*100,1),'%)'),
       y=paste0('Comp 2 (', round(neg_PLSDA$explained_variance$X[2]*100,1),'%)'))

all_CCP <- cor(all_PLSDA$X, all_PLSDA$variates$X, use='pairwise') %>% 
  data.frame() %>% 
  mutate(Ion = rownames(.),
         Corr = sqrt(comp.1^2+comp.2^2))  

f <- filter(all_CCP, Corr >= 0.5) %>% 
  ggplot(., aes(comp.1,comp.2, label=Ion))+
  geom_point()+
  geom_text(vjust=-1, size=3, color='navy')+
  geom_hline(yintercept = 0, lty=3)+
  geom_vline(xintercept = 0, lty=3)+
  geom_point(aes(x=0,y=0), color='darkgreen')+
  xlim(-1,1)+
  ylim(-1,1)+
  theme_bw()+
  theme(text=element_text(face = 'bold'))+
  labs(title = 'PLSDA model from all ions',
       x=paste0('Comp 1 (', round(neg_PLSDA$explained_variance$X[1]*100,1),'%)'),
       y=paste0('Comp 2 (', round(neg_PLSDA$explained_variance$X[2]*100,1),'%)'))

gridExtra::grid.arrange(e,f, nrow=1) 

pdf('PLS-DA models of all ions without isotopologues.pdf', width=14, height=10)
gridExtra::grid.arrange(a,b, nrow=1)
gridExtra::grid.arrange(c,d, nrow=1)
gridExtra::grid.arrange(e,f, nrow=1)
dev.off()
# making graphs with the ions of interest ---------------------------------

top30 <- arrange(all_CCP, desc(abs(Corr))) %>% 
  .[1:30,]

names(all_data)[15] <- 'Patient'
pdf('Top 30 predictors from PLS-DA model without isotopologues.pdf', width=14, height=10)
all_data %>% 
  unite(Pol_Mz_Rt, c(Polarity,medMz, medRt), sep='_') %>% 
  filter(Pol_Mz_Rt %in% top30$Ion) %>% 
  mutate(Pol_Mz_Rt = gsub('_',' ', Pol_Mz_Rt)) %>% 
  ggplot(., aes(Condition, Area, fill=Condition))+
  geom_boxplot()+
  facet_wrap(~Pol_Mz_Rt, scales='free')+
  theme_bw()+
  scale_fill_manual('Condition',values=colors)+
  labs(x='',y='Response (AU)', title='Top 30 predictors from PLS-DA model')+
  theme(text = element_text(face='bold'),
        axis.text.x = element_blank())
dev.off()  

  
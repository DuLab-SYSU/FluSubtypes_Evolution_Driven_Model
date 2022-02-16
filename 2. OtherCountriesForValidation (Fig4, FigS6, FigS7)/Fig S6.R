remove(list = ls())
library(ggplot2)
library(ggpubr)

GermanySeqSummary <- read.csv('GermanySeqSummary.csv')
GermanySeqSummary$month <- as.Date(GermanySeqSummary$month,format = '%Y-%m-%d')


AustriaSeqSummary <- read.csv('AustriaSeqSummary.csv')
AustriaSeqSummary$month <- as.Date(AustriaSeqSummary$month, format = '%Y-%m-%d')


NorwaySeqSummary <- read.csv('NorwaySeqSummary.csv')
NorwaySeqSummary$month <- as.Date(NorwaySeqSummary$month, format = '%Y-%m-%d')

starttime <- as.Date('2014-10-01',format = '%Y-%m-%d')
endtime <- as.Date('2019-09-01',format = '%Y-%m-%d')

fGermany <- GermanySeqSummary %>%
  filter(month>=starttime) %>%
  filter(month<=endtime) %>%
  ggplot(aes(x=month,y=sum, fill=Subtype))+
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(x='Time(monthly)', y='The number of sequence',title='Germany')+
  scale_fill_manual(values = c('LightBlue','LightCoral'))+
  theme_bw() 

fAustria <- AustriaSeqSummary %>%
  filter(month>=starttime) %>%
  filter(month<=endtime) %>%
  ggplot(aes(x=month,y=sum, fill=Subtype))+
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(x='Time(monthly)', y='The number of sequence',title='Austria')+
  scale_fill_manual(values = c('LightBlue','LightCoral'))+
  theme_bw() 

fNorway <- NorwaySeqSummary %>%
  filter(month>=starttime) %>%
  filter(month<=endtime) %>%
  ggplot(aes(x=month,y=sum, fill=Subtype))+
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(x='Time(monthly)', y='The number of sequence',title='Norway')+
  scale_fill_manual(values = c('LightBlue','LightCoral'))+
  theme_bw() 

ggarrange(fGermany, fAustria, fNorway, ncol=1, labels = c('A','B','C'))

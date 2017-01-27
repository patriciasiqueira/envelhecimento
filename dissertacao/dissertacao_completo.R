# Dados do Sul/Sudoeste de MG sobre envelhecimento
# AA usando variáveis originais e distância de Mahalanobis
# métodos: Ward e k-médias

# fase preliminar ---------------------------------------------------------
# diretório de trabalho
setwd("/home/patricia/Dropbox/nupis/projetos/2015-larissa")
# dados Sul/Sudoeste de Minas
load("sul_dem_nomes_antigos.rda")
# variáveis fectot no lugar de tft, razdep no lugar de rd etc.

# pacotes
library(MVN)          # normalidade multivariada
library(mvShapiroTest)# normalidade multivariada
library(corrplot)     # correlações
library(ecodist)      # Ward - D . Mahalanobis
library(ggplot2)      # gráficos
library(ggdendro)     # dendrogramas
library(dplyr)        # manipulação de dados
library(NbClust)      # Número de grupos ideais
library(cluster)      # pam
library(Rcpp)         # ACP  
library(FactoMineR)   # ACP

# limpar dados e ajeitar --------------------------------------------------
# tirar pop
sul <- sul[,-13]
# nomes dos municípios como rótulos
sul <- as.data.frame(sul)
row.names(sul) <- sul$nome.mun
# resumo estatístico
summary(sul[,-c(1:4)])
# CV
Xb <- apply(sul[,-c(1:4)],2,mean)
S <- apply(sul[,-c(1:4)],2,sd)
CV <- S/Xb*100
CV


# olhar separadamente para as variáveis
names(sul[,-c(1:4)])
boxplot(sul$espvida)
boxplot(sul$fectot)
sort(sul$fectot)
sul$nome.mun[which.min(sul$fectot)]
sul$nome.mun[which.max(sul$fectot)]
boxplot(sul$mort1)
boxplot(sul$mort5)
boxplot(sul$razdep)
boxplot(sul$sobre40)
boxplot(sul$sobre60)
boxplot(sul$t_env)
sort(sul$t_env)
sul$nome.mun[which.min(sul$t_env)]
sul$nome.mun[which.max(sul$t_env)]
# identificar outliers
# distância de Mahalanobis
result <- mvOutlier(sul[,-c(1:4)], qqplot=TRUE,method="quan")
result
# normalidade
uniPlot(sul[,-c(1:4)], type = "qqplot")  # QQPlots
uniPlot(sul[,-c(1:4)], type = "histogram") # hisogramas
uniNorm(sul[,-c(1:4)], type = "SW", desc = TRUE) # Shapiro-Wilk univariado
# teste multivariado de Royston 
roystonTest(sul[,-c(1:4)])
# Shapiro Wilk multivariado
X <- as.matrix(sul[,-c(1:4)])
mvShapiro.Test(X)
graphics.off()
resultado <- mardiaTest(X, qqplot = T)
resultado
# correlações
graphics.off()
cor(sul[,-c(1:4)])
sul.r <- cor(sul[,-c(1:4)])
round(sul.r,2)
corrplot.mixed(sul.r,tl.col = "black",
               number.digits = 4, number.font = 0.1,
               tl.cex = 1, is.corr = F)
# retirar as variáveis correlacionadas com outras (mort5 e sobre40)
# e que, de alguma forma, calculam a mesma coisa
names(sul)
sul <- sul[,-c(8,10)]
# correlações após retirar as variáveis
graphics.off()
cor(sul[,-c(1:4)])
sul.r <- cor(sul[,-c(1:4)])
corrplot.mixed(sul.r,tl.col = "black",
               number.digits = 4, number.font = 0.1,
               tl.cex = 1, is.corr = F)

# variáveis originais - Mahalanobis ---------------------------------------

# matriz de distâncias - Mahalanobis
dis <- (distance(sul[,-c(1:4)],method="mahalanobis"))^0.5   #Mahalanobis

# usar distância de Mahalanobis
# Ward
wr <- hclust(dis, method = "ward.D2")
ggdendrogram(wr)
# dendrograma com opções - Ward
hc <- dis %>% hclust(method = "ward.D2") 
ggdendrogram(hc)
ddata <- hc %>% as.dendrogram() %>% dendro_data()
ggdendrogram(hc) + geom_text(size= 2, aes(x = x, y = y, label = label, 
                                          angle = -90, vjust=0, hjust = 0), data= label(ddata)) +
  scale_y_continuous(expand = c(0.6, 0)) +
  scale_x_continuous(expand = c(0.1, 0)) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) 
# com distâncias
ggdendrogram(hc) + geom_text(size= 2, aes(x = x, y = y, label = label, 
                                          angle = -90, vjust=0, hjust = 0), data= label(ddata)) +
  scale_y_continuous(expand = c(0.6, 0)) +
  scale_x_continuous(expand = c(0.1, 0)) +
  theme(axis.text.x = element_blank()) 


# gráficos - dois primeiros CPs - Mahalanobis
X_cp <- princomp(sul[,-c(1:4)], cor=T)
summary(X_cp) # os dois primeiros contabilizam 76% da var.
X_cp$loadings
escores <- X_cp$scores[,1:2]
escores <- as.data.frame(escores)
# usando ggplot2
graf_CP <- ggplot(data=escores,aes(x=escores[,1],y=escores[,2],label=rownames(escores)))
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour="black",alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw()+ylab("componente principal 2") +
  xlab("componente principal 1")

# gráficos de 2 grupos
# ward
lab2_wr <- cutree(wr, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_wr,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# avaliação dos grupos obtidos - k=2 ------------------------
# ward
g2_wr <- cutree(wr, k=2); table(g2_wr)
# Municípios dentro de cada grupo
sapply(unique(g2_wr),function(g)row.names(sul[,-c(1:4)])[g2_wr == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g2_wr),mean)
aggregate(sul[,-c(1:4)],list(g2_wr),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g2_wr),mean)
S <- aggregate(sul[,-c(1:4)],list(g2_wr),sd)
CV <- S/Xb*100
CV
# boxplots - ward
# espvida
sul2_wr <- sul
sul2_wr$grupo <- lab2_wr
ggplot(data = sul2_wr, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_wr, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_wr, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_wr, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_wr, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_wr, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))

# gráficos - Mahalanobis - k=3 --------------------------------------------
# ward
lab3_wr <- cutree(wr, k=3)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_wr,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")

# avaliação dos grupos obtidos - Mahalanobis - k=3 ------------------------
# ward
g3_wr <- cutree(wr, k=3); table(g3_wr)
# Municípios dentro de cada grupo
sapply(unique(g3_wr),function(g)row.names(sul[,-c(1:4)])[g3_wr == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g3_wr),mean)
aggregate(sul[,-c(1:4)],list(g3_wr),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g3_wr),mean)
S <- aggregate(sul[,-c(1:4)],list(g3_wr),sd)
CV <- S/Xb*100
CV
# boxplots
# espvida
sul3_wr <- sul
sul3_wr$grupo <- lab3_wr
ggplot(data = sul3_wr, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_wr, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_wr, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_wr, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_wr, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_wr, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))

# avaliação dos grupos obtidos - Mahalanobis - k=4 ------------------------
# ward
g4_wr <- cutree(wr, k=4); table(g4_wr)
# Municípios dentro de cada grupo
sapply(unique(g4_wr),function(g)row.names(sul[,-c(1:4)])[g4_wr == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g4_wr),mean)
aggregate(sul[,-c(1:4)],list(g4_wr),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g4_wr),mean)
S <- aggregate(sul[,-c(1:4)],list(g4_wr),sd)
CV <- S/Xb*100
CV
# boxplots
lab4_wr <- cutree(wr, k=4)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab4_wr,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# boxplots
# espvida
sul4_wr <- sul
sul4_wr$grupo <- lab4_wr
ggplot(data = sul4_wr, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul4_wr, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul4_wr, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul4_wr, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul4_wr, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul4_wr, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))

# ward, k = 5
# ward
lab5_wr <- cutree(wr, k=5)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab5_wr,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
g5_wr <- cutree(wr, k=5); table(g5_wr)
# Municípios dentro de cada grupo
sapply(unique(g5_wr),function(g)row.names(sul[,-c(1:4)])[g5_wr == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g5_wr),mean)
aggregate(sul[,-c(1:4)],list(g5_wr),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g5_wr),mean)
S <- aggregate(sul[,-c(1:4)],list(g5_wr),sd)
CV <- S/Xb*100
CV
# boxplots
# espvida
sul5_wr <- sul
sul5_wr$grupo <- lab5_wr
ggplot(data = sul5_wr, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul5_wr, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul5_wr, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul5_wr, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul5_wr, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul5_wr, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))

# k-médias - k=2, k=3, k=4, k=5 ----------------------------------
set.seed(1)
autoplot(kmeans(sul[,-c(1:4)], 2, nstart=25), data = sul, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")
set.seed(1)
autoplot(kmeans(sul[,-c(1:4)], 3, nstart=25), data = sul, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")
set.seed(1)
autoplot(kmeans(sul[,-c(1:4)], 4, nstart=25), data = sul, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")
# kmeans com ggplot
set.seed(1)
k_2 <- kmeans(sul[,-c(1:4)], 2, nstart=25)
labk2 <- k_2$cluster
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=labk2,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# avaliar grupos obtidos
table(labk2)
# Municípios dentro de cada grupo
sapply(unique(labk2),function(g)row.names(sul[,-c(1:4)])[labk2 == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(labk2),mean)
aggregate(sul[,-c(1:4)],list(labk2),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(labk2),mean)
S <- aggregate(sul[,-c(1:4)],list(labk2),sd)
CV <- S/Xb*100
CV
# boxplots
# k-médias, k = 2
# espvida
sul2_km <- sul
sul2_km$grupo <- labk2
ggplot(data = sul2_km, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_km, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_km, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_km, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_km, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_km, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# gráfico das SQ
wssplot <- function(data, nc=10, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}
wssplot(sul[,-c(1:4)])

# k = 3
set.seed(1)
k_3 <- kmeans(sul[,-c(1:4)], 3, nstart=25)
labk3 <- k_3$cluster
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=labk3,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# avaliar grupos obtidos
table(labk3)
# Municípios dentro de cada grupo
sapply(unique(labk3),function(g)row.names(sul[,-c(1:4)])[labk3 == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(labk3),mean)
aggregate(sul[,-c(1:4)],list(labk3),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(labk3),mean)
S <- aggregate(sul[,-c(1:4)],list(labk3),sd)
CV <- S/Xb*100
CV
# boxplots
# k-médias, k = 3
# espvida
sul3_km <- sul
sul3_km$grupo <- labk3
ggplot(data = sul3_km, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_km, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_km, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_km, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_km, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_km, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))

# k = 4
set.seed(1)
k_4 <- kmeans(sul[,-c(1:4)], 4, nstart=25)
labk4 <- k_4$cluster
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=labk4,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# avaliar grupos obtidos
table(labk4)
# Municípios dentro de cada grupo
sapply(unique(labk4),function(g)row.names(sul[,-c(1:4)])[labk4 == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(labk4),mean)
aggregate(sul[,-c(1:4)],list(labk4),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(labk4),mean)
S <- aggregate(sul[,-c(1:4)],list(labk4),sd)
CV <- S/Xb*100
CV
# boxplots
# k-médias, k = 4
# espvida
sul4_km <- sul
sul4_km$grupo <- labk4
ggplot(data = sul4_km, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul4_km, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul4_km, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul4_km, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul4_km, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul4_km, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))

# k = 5
set.seed(1)
k_5 <- kmeans(sul[,-c(1:4)], 5, nstart=25)
labk5 <- k_5$cluster
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=labk5,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# avaliar grupos obtidos
table(labk5)
# Municípios dentro de cada grupo
sapply(unique(labk5),function(g)row.names(sul[,-c(1:4)])[labk5 == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(labk5),mean)
aggregate(sul[,-c(1:4)],list(labk5),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(labk5),mean)
S <- aggregate(sul[,-c(1:4)],list(labk5),sd)
CV <- S/Xb*100
CV
# boxplots
# k-médias, k = 5
# espvida
sul5_km <- sul
sul5_km$grupo <- labk5
ggplot(data = sul5_km, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul5_km, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul5_km, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul5_km, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul5_km, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul5_km, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))






# Número de grupos --------------------------------------------------------
set.seed(1)
ng <- NbClust(sul[,-c(1:4)], distance = NULL, 
              diss = dis, method="ward.D2")
ng$Best.nc
table(ng$Best.nc)
barplot(table(ng$Best.nc[1,]),
        xlab="Número de grupos", ylab="Número de critérios")
# gap
ng_gap <- NbClust(sul[,-c(1:4)], distance = NULL, 
                  diss = dis, method="ward.D2", index = "gap")
ng_gap$Best.nc
# ch
ng_ch <- NbClust(sul[,-c(1:4)], distance = NULL, 
                  diss = dis, method="ward.D2", index = "ch")
ng_ch$Best.nc

# k-médias
set.seed(1)  # para que os resultados sempre sejam os mesmos
ng_km <- NbClust(sul[,-c(1:4)], min.nc=2, max.nc=15, method="kmeans")
ng_km$Best.nc
table(ng_km$Best.nc)
barplot(table(ng_km$Best.nc[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos por 26 critérios")

# gap
ng_gap <- NbClust(sul[,-c(1:4)], method="kmeans", index = "gap")
ng_gap$Best.nc
# ch
ng_ch <- NbClust(sul[,-c(1:4)],method="kmeans", index = "ch")
ng_ch$Best.nc

# pam
op <- par(mfrow=c(1,2))
# 2 grupos
pam2 <- pam(sul[,-c(1:4)], 2)
plot(pam2)
# 3 grupos
pam3 <- pam(sul[,-c(1:4)], 3)
plot(pam3)
# 4 grupos
pam4 <- pam(sul[,-c(1:4)], 4)
plot(pam4)
# 5 grupos
pam5 <- pam(sul[,-c(1:4)], 5)
plot(pam5)

# para obter o silhouette plot dos agrupamentos de Ward 
# e k-médias lado a lado
# Segundo Everitt (cluster analysis), quando s(i) está próximo 
# de 1, a heterogeneidade do grupo é muito menor do que sua
# separação e i é considerado bem classificado
# quando s(i) está próximo de -1, o grupo está mal classificado
# quando próximo de 0 não está claro se o objeto deveria estar 
# no grupo ou no próximo
dissE <- daisy(sul[,-c(1:4)]) # distância euclidiana
op <- par(mfrow=c(1,2))
plot(silhouette(cutree(wr,2),dis)) # ward, k = 2
plot(silhouette(labk2, dissE))  # k-médias, k = 2
plot(silhouette(cutree(wr,3),dis)) # ward, k = 3
plot(silhouette(labk3, dissE))  # k-médias, k = 3
plot(silhouette(cutree(wr,4),dis)) # ward, k = 4
plot(silhouette(labk4, dissE))  # k-médias, k = 4
plot(silhouette(cutree(wr,5),dis)) # ward, k = 5
plot(silhouette(labk5, dissE))  # k-médias, k = 5


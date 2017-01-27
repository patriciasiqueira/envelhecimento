# Dados do Sul/Sudoeste de MG sobre envelhecimento
# AA usando variáveis originais e escores dos componentes principais

# fase preliminar ---------------------------------------------------------

# diretório de trabalho --------------------------------------------------
setwd("/home/patricia/Dropbox/nupis/projetos/2015-larissa")

# dados Sul/Sudoeste de Minas -------------------------------------------------------------------
load("sul_dem.rda")

# pacotes -----------------------------------------------------------------
library(corrplot)     # Correlações
library(ClustOfVar)   # Agrupamento de Variáveis
library(Rcpp)         # ACP  
library(FactoMineR)   # ACP
library(stringi)      
library(ggplot2)      # Gráficos
library(ggdendro)     # Dendrogramas
library(NbClust)      # Número de grupos ideais
library(ggfortify)    # Agrupamentos
library(ecodist)      # Ward - D . Mahalanobis
library(dplyr)        # manipulação de dados
library(cluster)      # pam
library(MVN)          # normalidade multivariada
library(mvShapiroTest)# normalidade multivariada

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

# olhar separadamente para as variávies
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
corrplot(sul.r,tl.col = "black",type="lower")
corrplot.mixed(sul.r, tl.col = "black")

# ordem diferente
corrplot(sul.r,tl.col = "black", method = "number",
         order = "AOE")
corrplot(sul.r,tl.col = "black", method = "number",
         order = "AOE", number.digits = 3, number.font = 0.1,
         tl.cex = 1, is.corr = F)
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

# realizar todos os agrupamentos hierárquicos usando
# distância de Mahalanobis
# single linkage
vp <- hclust(dis, method = "single")
ggdendrogram(vp)
# complete linkage
vd <- hclust(dis, method = "complete")
ggdendrogram(vd)
# distância média
dm <- hclust(dis, method = "average")
ggdendrogram(dm)
# centroide
ct <- hclust(dis, method = "centroid")
ggdendrogram(ct)
# Ward
wr <- hclust(dis, method = "ward.D2")
ggdendrogram(wr)

# dendrograma com opções - centroide
hc <- dis %>% hclust(method = "centroid") 
ddata <- hc %>% as.dendrogram() %>% dendro_data()
ggdendrogram(hc) + geom_text(size= 2, aes(x = x, y = y, label = label, 
                                          angle = -90, vjust=0, hjust = 0), data= label(ddata)) +
  scale_y_continuous(expand = c(0.6, 0)) +
  scale_x_continuous(expand = c(0.1, 0)) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) 
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


# gráficos - dois primeiros CPs - Mahalanobis - k=2 -----------------------
X_cp <- princomp(sul[,-c(1:4)], cor=T)
summary(X_cp) # os dois primeiros contabilizam 76% da var.
X_cp$loadings
# R básico
# plot(X_cp$scores[,1:2], type = "n")
# text(X_cp$scores[,1:2], cex = 1, label = row.names(sul))
escores <- X_cp$scores[,1:2]
# usando ggplot2
graf_CP <- ggplot(data=escores,aes(x=escores[,1],y=escores[,2],label=rownames(escores)))
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour="black",alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw()+ylab("componente principal 2") +
  xlab("componente principal 1")

# gráficos de 2 grupos
# vizinho mais próximo
lab2_vp <- cutree(vp, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_vp,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# vizinho mais distante
lab2_vd <- cutree(vd, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_vd,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# distância média
lab2_dm <- cutree(dm, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_dm,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# centroide
lab2_ct <- cutree(ct, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_ct,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# ward
lab2_wr <- cutree(wr, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_wr,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")

# avaliação dos grupos obtidos - Mahalanobis - k=2 ------------------------
# ligação simples
g2_vp <- cutree(vp, k=2); table(g2_vp)
# Municípios dentro de cada grupo
sapply(unique(g2_vp),function(g)row.names(sul[,-c(1:4)])[g2_vp == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g2_vp),mean)
aggregate(sul[,-c(1:4)],list(g2_vp),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g2_vp),mean)
S <- aggregate(sul[,-c(1:4)],list(g2_vp),sd)
CV <- S/Xb*100
CV
# ligação completa
g2_vd <- cutree(vd, k=2); table(g2_vd)
# Municípios dentro de cada grupo
sapply(unique(g2_vd),function(g)row.names(sul[,-c(1:4)])[g2_vd == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g2_vd),mean)
aggregate(sul[,-c(1:4)],list(g2_vd),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g2_vd),mean)
S <- aggregate(sul[,-c(1:4)],list(g2_vd),sd)
CV <- S/Xb*100
CV
# distância média
g2_dm <- cutree(dm, k=2); table(g2_dm)
# Municípios dentro de cada grupo
sapply(unique(g2_dm),function(g)row.names(sul[,-c(1:4)])[g2_dm == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g2_dm),mean)
aggregate(sul[,-c(1:4)],list(g2_dm),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g2_dm),mean)
S <- aggregate(sul[,-c(1:4)],list(g2_dm),sd)
CV <- S/Xb*100
CV
# centroide
g2_ct <- cutree(ct, k=2); table(g2_ct)
# Municípios dentro de cada grupo
sapply(unique(g2_ct),function(g)row.names(sul[,-c(1:4)])[g2_ct == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g2_ct),mean)
aggregate(sul[,-c(1:4)],list(g2_ct),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g2_ct),mean)
S <- aggregate(sul[,-c(1:4)],list(g2_ct),sd)
CV <- S/Xb*100
CV
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

# boxplots
# vizinho mais próximo
# espvida
sul2_vp <- sul
sul2_vp$grupo <- lab2_vp
ggplot(data = sul2_vp, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_vp, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_vp, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_vp, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_vp, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_vp, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# vizinho mais distante
# espvida
sul2_vd <- sul
sul2_vd$grupo <- lab2_vd
ggplot(data = sul2_vd, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_vd, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_vd, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_vd, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_vd, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_vd, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# distância média
# espvida
sul2_dm <- sul
sul2_dm$grupo <- lab2_dm
ggplot(data = sul2_dm, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_vd, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_vd, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_vd, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_vd, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_vd, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# centroide
# espvida
sul2_ct <- sul
sul2_ct$grupo <- lab2_ct
ggplot(data = sul2_ct, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_ct, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_ct, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_ct, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_ct, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_ct, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# ward
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
# vizinho mais próximo
lab3_vp <- cutree(vp, k=3)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_vp,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# vizinho mais distante
lab3_vd <- cutree(vd, k=3)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_vd,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# distância média
lab3_dm <- cutree(dm, k=3)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_vp,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# centroide
lab3_ct <- cutree(ct, k=3)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_ct,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# ward
lab3_wr <- cutree(wr, k=3)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_wr,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")

# avaliação dos grupos obtidos - Mahalanobis - k=3 ------------------------
# ligação simples
g3_vp <- cutree(vp, k=3); table(g3_vp)
# Municípios dentro de cada grupo
sapply(unique(g3_vp),function(g)row.names(sul[,-c(1:4)])[g3_vp == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g3_vp),mean)
aggregate(sul[,-c(1:4)],list(g3_vp),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g3_vp),mean)
S <- aggregate(sul[,-c(1:4)],list(g3_vp),sd)
CV <- S/Xb*100
CV
# boxplots
# vizinho mais próximo
# espvida
sul3_vp <- sul
sul3_vp$grupo <- lab3_vp
ggplot(data = sul3_vp, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_vp, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_vp, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_vp, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_vp, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_vp, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# vizinho mais distante
# espvida
sul3_vd <- sul
sul3_vd$grupo <- lab3_vd
ggplot(data = sul3_vd, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_vd, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_vd, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_vd, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_vd, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_vd, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# distância média
# espvida
sul3_dm <- sul
sul3_dm$grupo <- lab3_dm
ggplot(data = sul3_dm, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_vd, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_vd, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_vd, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_vd, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_vd, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# centroide
# espvida
sul3_ct <- sul
sul3_ct$grupo <- lab3_ct
ggplot(data = sul3_ct, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_ct, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_ct, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_ct, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_ct, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_ct, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# ward
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

# gráficos de 4 grupos
# vizinho mais próximo
lab <- cutree(vp, k=4)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# vizinho mais distante
lab <- cutree(vd, k=4)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# distância média
lab <- cutree(dm, k=4)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# centroide
lab <- cutree(ct, k=4)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# ward
lab <- cutree(wr, k=4)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")

# ward, k = 4
# ward
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
# ligação simples
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
ggplot(data = sul2_vp, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))

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
# ligação simples
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
# ligação simples
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
# ligação simples
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


# variáveis padronizadas - distância euclidiana ---------------------------
# padronizar as variáveis - subtraindo a média e dividindo pelo
# desvio padrão
sul.p <- scale(sul[,-c(1:4)],center=T,scale=T)
sul.p <- as.data.frame(sul.p)
# sul.p é o data frame com as variáveis padronizadas
# transformar nomes dos municípios em row.names
row.names(sul.p) <- sul$nome.mun
names(sul.p)


# gráficos - 2 CPs --------------------------------------------------------
X_cp_p <- princomp(sul.p, cor=T)
summary(X_cp_p) # os dois primeiros contabilizam 76% da var.
escores_p <- X_cp_p$scores[,1:2]
# usando ggplot2
graf_CP_e <- ggplot(data=escores_p,
                  aes(x=escores_p[,1],y=escores_p[,2],
                      label=rownames(escores_p)))
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour="black",alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw()+ylab("componente principal 2") +
  xlab("componente principal 1")

# realizar todos os agrupamentos hierárquicos usando
# dados padronizados e distância euclidiana
dis_e <- dist(sul.p)
# single linkage
vp_e <- hclust(dis_e, method = "single")
ggdendrogram(vp_e)
# complete linkage
vd_e <- hclust(dis_e, method = "complete")
ggdendrogram(vd_e)
# distância média
dm_e <- hclust(dis_e, method = "average")
ggdendrogram(dm_e)
# centroide
ct_e <- hclust(dis_e, method = "centroid")
ggdendrogram(ct_e)
# Ward
wr_e <- hclust(dis_e, method = "ward.D2")
ggdendrogram(wr_e)


# gráficos - 2 CPs - dados padronizados - euclidiana - k=2 ----------------
# vizinho mais próximo
lab_e <- cutree(vp_e, k=2)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# vizinho mais distante
lab_e <- cutree(vd_e, k=2)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# distância média
lab_e <- cutree(dm_e, k=2)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# centroide
lab_e <- cutree(ct_e, k=2)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# ward
lab_e <- cutree(wr_e, k=2)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")

# gráfico - dois primeiros CPs
# k = 2
X_cp <- princomp(sul[,-c(1:4)], cor=T)
summary(X_cp) # os dois primeiros contabilizam 76% da var.
X_cp$loadings
# R básico
# plot(X_cp$scores[,1:2], type = "n")
# text(X_cp$scores[,1:2], cex = 1, label = row.names(sul))
escores <- X_cp$scores[,1:2]
# usando ggplot2
graf_CP <- ggplot(data=escores,aes(x=escores[,1],y=escores[,2],label=rownames(escores)))
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour="black",alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw()+ylab("componente principal 2") +
  xlab("componente principal 1")

# gráficos de 2 grupos
# vizinho mais próximo
lab2_vp <- cutree(vp_e, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_vp,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# vizinho mais distante
lab2_vd <- cutree(vd_e, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_vd,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# distância média
lab2_dm <- cutree(dm_e, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_dm,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# centroide
lab2_ct <- cutree(ct_e, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_ct,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# ward
lab2_wr <- cutree(wr_e, k=2)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab2_wr,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")

# avaliação dos grupos obtidos - euclidiana - k=2 -------------------------
# ligação simples
g2_vp <- cutree(vp_e, k=2); table(g2_vp)
# Municípios dentro de cada grupo
sapply(unique(g2_vp),function(g)row.names(sul.p)[g2_vp == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g2_vp),mean)
aggregate(sul[,-c(1:4)],list(g2_vp),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g2_vp),mean)
S <- aggregate(sul[,-c(1:4)],list(g2_vp),sd)
CV <- S/Xb*100
CV
# boxplots
# vizinho mais próximo
# espvida
sul2_vp <- sul
sul2_vp$grupo <- lab2_vp
ggplot(data = sul2_vp, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_vp, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_vp, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_vp, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_vp, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_vp, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# vizinho mais distante
# espvida
sul2_vd <- sul
sul2_vd$grupo <- lab2_vd
ggplot(data = sul2_vd, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_vd, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_vd, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_vd, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_vd, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_vd, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# distância média
# espvida
sul2_dm <- sul
sul2_dm$grupo <- lab2_dm
ggplot(data = sul2_dm, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_vd, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_vd, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_vd, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_vd, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_vd, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# centroide
# espvida
sul2_ct <- sul
sul2_ct$grupo <- lab2_ct
ggplot(data = sul2_ct, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul2_ct, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul2_ct, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul2_ct, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul2_ct, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul2_ct, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# ward
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

# gráficos de 3 grupos
# vizinho mais próximo
lab3_vp <- cutree(vp_e, k=3)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_vp,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# vizinho mais distante
lab3_vd <- cutree(vd_e, k=3)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_vd,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# distância média
lab3_dm <- cutree(dm_e, k=3)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_dm,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# centroide
lab3_ct <- cutree(ct_e, k=3)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_ct,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# ward
lab3_wr <- cutree(wr_e, k=3)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab3_wr,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")

# avaliação dos grupos obtidos - euclidiana - k=3 -------------------------
# ligação simples
g3_vp <- cutree(vp_e, k=3); table(g3_vp)
# Municípios dentro de cada grupo
sapply(unique(g3_vp),function(g)row.names(sul.p)[g3_vp == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(g3_vp),mean)
aggregate(sul[,-c(1:4)],list(g3_vp),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(g3_vp),mean)
S <- aggregate(sul[,-c(1:4)],list(g3_vp),sd)
CV <- S/Xb*100
CV
# boxplots
# vizinho mais próximo
# espvida
sul3_vp <- sul
sul3_vp$grupo <- lab3_vp
ggplot(data = sul3_vp, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_vp, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_vp, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_vp, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_vp, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_vp, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# vizinho mais distante
# espvida
sul3_vd <- sul
sul3_vd$grupo <- lab3_vd
ggplot(data = sul3_vd, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_vd, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_vd, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_vd, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_vd, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_vd, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# distância média
# espvida
sul3_dm <- sul
sul3_dm$grupo <- lab3_dm
ggplot(data = sul3_dm, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_vd, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_vd, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_vd, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_vd, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_vd, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# centroide
# espvida
sul3_ct <- sul
sul3_ct$grupo <- lab3_ct
ggplot(data = sul3_ct, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))
# fectot
ggplot(data = sul3_ct, aes(as.factor(grupo), fectot)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$fectot),max(sul$fectot)))
# mort1
ggplot(data = sul3_ct, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
ggplot(data = sul3_ct, aes(as.factor(grupo), razdep)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$razdep),max(sul$razdep)))
# sobre60
ggplot(data = sul3_ct, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
ggplot(data = sul3_ct, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))
# ward
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

# gráficos de 4 grupos
# vizinho mais próximo
lab_e <- cutree(vp_e, k=4)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# vizinho mais distante
lab_e <- cutree(vd_e, k=4)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# distância média
lab_e <- cutree(dm_e, k=4)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# centroide
lab_e <- cutree(ct_e, k=4)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# ward
lab_e <- cutree(wr_e, k=4)
graf_CP_e +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab_e,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")

# kmeans
set.seed(1)
autoplot(kmeans(sul.p, 2, nstart=25), data = sul.p, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")
set.seed(1)
autoplot(kmeans(sul.p, 3, nstart=25), data = sul.p, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")
set.seed(1)
autoplot(kmeans(sul.p, 4, nstart=25), data = sul.p, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")
# k-médias com ggplot
# k = 2, dados padronizados
set.seed(1)
k_2_p <- kmeans(sul.p, 2, nstart=25)
labk2 <- k_2_p$cluster
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=labk2,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# boxplots
# k-médias, padronizados, k = 2
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

# k-médias, padronizada, k = 3
set.seed(1)
k_3_p <- kmeans(sul.p, 3, nstart=25)
labk3 <- k_3_p$cluster
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=labk3,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")
# boxplots
# k-médias, padronizada, k = 3
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





set.seed(1)
k_4_p <- kmeans(sul.p, 4, nstart=25)
lab <- k_4_p$cluster
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")








# avaliar grupos obtidos de Ward
grupos.2 <- cutree(wr, k=2); table(grupos.2)
# Municípios dentro de cada grupo
sapply(unique(grupos.2),function(g)row.names(sul[,-c(1:4)])[grupos.2 == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(grupos.2),mean)
aggregate(sul[,-c(1:4)],list(grupos.2),median)
# CV
Xb <- aggregate(sul[,-c(1:4)],list(grupos.2),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.2),sd)
CV <- S/Xb*100
CV

# avaliar os 3 grupos de Ward
grupos.3 <- cutree(wr, k=3); table(grupos.3)
# Municípios dentro de cada grupo
sapply(unique(grupos.3),function(g)row.names(sul[,-c(1:4)])[grupos.3 == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(grupos.3),mean)
aggregate(sul[,-c(1:4)],list(grupos.3),median)

# Os grupos 2 e 3 estão muito parecidos!!!

# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.3),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.3),sd)
CV <- S/Xb*100
CV

# http://strata.uga.edu/software/pdf/clusterTutorial.pdf
library(cluster)     
# padrão
mydata.agnes <- agnes(sul[,-c(1:4)])
plot(mydata.agnes, which.plots = 2, main= "Dendrogram of my data")
# manhattan, Ward
mydata.agnes.ALT <- agnes(sul[,-c(1:4)], metric="manhattan", method="ward")
plot(mydata.agnes.ALT, which.plots = 2, main= "Dendrogram of my data")
# euclidiana, Ward
mydata.agnes.ALT <- agnes(sul[,-c(1:4)], metric="euclidian", stand = F, method="ward")
plot(mydata.agnes.ALT, which.plots = 2, main= "Dendrogram of my data")
# mahalanobis, ward
dis.m <- (distance(sul[,-c(1:4)],method="mahalanobis"))^0.5   #Mahalanobis
agnes.mw <- agnes(dis.m, diss = F, method="ward", stand = F)
plot(agnes.mw, which.plots = 4, main= "Dendrogram of my data")
# pacote vegan
library(vegan)
# mahalanobis, ward - resultado diferente
dist.vegan.m <- vegdist(sul[,-c(1:4)], method="mahalanobis")
mydata.bray.agnes <- agnes(dist.vegan.m, method = "ward", stand = F)
plot(mydata.bray.agnes, which.plots = 2, main= "Dendrogram of my data")


# variáveis padronizadas - euclidiana -------------------------------------

# usando dados padronizados e distância euclidiana
dis <- dist(sul.p)

# métodos hierárquicos

# ligação simples
vp <- hclust(dis, method = "single")
library(ggdendro)
ggdendrogram(vp)

# escolha do número de grupos
# Número de grupos pelo critérios
set.seed(1)
ng <- NbClust(sul.p, method="single")

ng$Best.nc

# Retirar os métodos Dindex e Hubert para fazer barplot porque são gráficos 
table(ng$Best.nc[1,c(-23,-25)])
barplot(table(ng$Best.nc[1,c(-23,-25)]),
        xlab="Número de grupos", ylab="Número de critérios")


# gamma
ng_gamma<- NbClust(sul.p, method="single",index = "gamma")
ng_gamma$Best.nc

# gplus
ng_gplus<- NbClust(sul.p, method="single", index = "gplus")
ng_gplus$Best.nc

# gap
ng_gap <- NbClust(sul.p, method="single", index = "gap")
ng_gap$Best.nc

# 5 critérios: 2 grupos, maioria dos métodos: 5 grupos

# ch:5 ; duda:2; cindex:4; gamma:9; beale:2 ; ccc:2 ; PtBiserial:12 ;
# gplus:2; DB:12; gap:2


# obter o número de objetos dentro de cada grupo
grupos.2 <- cutree(vp, k=2); table(grupos.2)
grupos.5 <- cutree(vp, k=5); table(grupos.5)
# Municípios dentro de cada grupo
sapply(unique(grupos.2),function(g)row.names(sul.p)[grupos.2 == g])
sapply(unique(grupos.5),function(g)row.names(sul.p)[grupos.5 == g])

#Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(grupos.2),median)
aggregate(sul[,-c(1:4)],list(grupos.5),median)

# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.2),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.2),sd)
CV <- S/Xb*100
CV

# 5 grupos
Xb <- aggregate(sul[,-c(1:4)],list(grupos.5),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.5),sd)
CV <- S/Xb*100
CV


# Ligação completa
vd <- hclust(dis, method = "complete")
ggdendrogram(vd)

# escolha do número de grupos
# Número de grupos pelo critérios
set.seed(1)
ng <- NbClust(sul.p, method="complete")

ng$Best.nc

# Retirar os métodos Dindex e Hubert para fazer barplot porque são gráficos 
table(ng$Best.nc[1,c(-23,-25)])
barplot(table(ng$Best.nc[1,c(-23,-25)]),
        xlab="Número de grupos", ylab="Número de critérios")


# gamma
ng_gamma<- NbClust(sul.p, method="complete",index = "gamma")
ng_gamma$Best.nc

# gplus
ng_gplus<- NbClust(sul.p, method="complete", index = "gplus")
ng_gplus$Best.nc

# gap
ng_gap <- NbClust(sul.p, method="complete", index = "gap")
ng_gap$Best.nc

# 4 critérios: 4 grupos , maioria dos métodos 3 métodos

# ch:2; duda:4; cindex:3; gamma:15; beale:4; ccc:2; PtBiserial:4;
# gplus:15; DB:4; gap:2


# obter o número de objetos dentro de cada grupo
grupos.4 <- cutree(vp, k=4); table(grupos.4)

# Municípios dentro de cada grupo
sapply(unique(grupos.4),function(g)row.names(sul.p)[grupos.4 == g])


#Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(grupos.4),median)

# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.4),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.4),sd)
CV <- S/Xb*100
CV

# distancia media
dm <- hclust(dis, method = "average")
library(ggdendro)
ggdendrogram(dm)

# escolha do número de grupos
# Número de grupos pelo critérios
set.seed(1)
ng <- NbClust(sul.p, method="average")

ng$Best.nc

# Retirar os métodos Dindex e Hubert para fazer barplot porque são gráficos 
table(ng$Best.nc[1,c(-23,-25)])
barplot(table(ng$Best.nc[1,c(-23,-25)]),
        xlab="Número de grupos", ylab="Número de critérios")


# gamma
ng_gamma<- NbClust(sul.p, method="average",index = "gamma")
ng_gamma$Best.nc

# gplus
ng_gplus<- NbClust(sul.p, method="average", index = "gplus")
ng_gplus$Best.nc

# gap
ng_gap <- NbClust(sul.p, method="average", index = "gap")
ng_gap$Best.nc

# 7 critérios:2 grupos , maioria dos métodos 2 métodos

# ch:5; duda:2; cindex:2; gamma:2; beale:2; ccc:15; PtBiserial:7;
# gplus:2; DB:2; gap:2


# obter o número de objetos dentro de cada grupo
grupos.2 <- cutree(dm, k=2); table(grupos.2)

# Municípios dentro de cada grupo
sapply(unique(grupos.2),function(g)row.names(sul.p)[grupos.2 == g])


#Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(grupos.2),median)

# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.2),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.2),sd)
CV <- S/Xb*100
CV


# Centróide
ct <- hclust(dis, method = "centroid")
ggdendrogram(ct)

# escolha do número de grupos
# Número de grupos pelo critérios
set.seed(1)
ng <- NbClust(sul.p, method="centroid")

ng$Best.nc

# Retirar os métodos Dindex e Hubert para fazer barplot porque são gráficos 
table(ng$Best.nc[1,c(-23,-25)])
barplot(table(ng$Best.nc[1,c(-23,-25)]),
        xlab="Número de grupos", ylab="Número de critérios")


# gamma
ng_gamma<- NbClust(sul.p, method="centroid",index = "gamma")
ng_gamma$Best.nc

# gplus
ng_gplus<- NbClust(sul.p, method="centroid", index = "gplus")
ng_gplus$Best.nc

# gap
ng_gap <- NbClust(sul.p, method="centroid", index = "gap")
ng_gap$Best.nc

# 7 critérios: 2 grupos , maioria dos métodos  métodos

# ch:2; duda:2; cindex:4; gamma:2; beale:2; ccc:2; PtBiserial:12;
# gplus:2; DB:4; gap:2


# obter o número de objetos dentro de cada grupo
grupos.2 <- cutree(ct, k=2); table(grupos.2)

# Municípios dentro de cada grupo
sapply(unique(grupos.2),function(g)row.names(sul.p)[grupos.2 == g])


#Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(grupos.2),median)

# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.2),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.2),sd)
CV <- S/Xb*100
CV

# v) Ward
wr <- hclust(dis, method = "ward.D2")
ggdendrogram(wr)

# escolha do número de grupos
# Número de grupos pelo critérios
set.seed(1)
ng <- NbClust(sul.p, method="ward.D2")

ng$Best.nc

# Retirar os métodos Dindex e Hubert para fazer barplot porque são gráficos 
table(ng$Best.nc[1,c(-23,-25)])
barplot(table(ng$Best.nc[1,c(-23,-25)]),
        xlab="Número de grupos", ylab="Número de critérios")


# gamma
ng_gamma<- NbClust(sul.p, method="ward.D2",index = "gamma")
ng_gamma$Best.nc

# gplus
ng_gplus<- NbClust(sul.p, method="ward.D2", index = "gplus")
ng_gplus$Best.nc

# gap
ng_gap <- NbClust(sul.p, method="ward.D2", index = "gap")
ng_gap$Best.nc

#  critérios:  grupos , maioria dos métodos 3 grupos

# ch:2; duda:6; cindex:3; gamma:14; beale:14; ccc:15; PtBiserial:6;
# gplus:15; DB:13; gap:2


# obter o número de objetos dentro de cada grupo
grupos.3 <- cutree(wr, k=3); table(grupos.3)

# Municípios dentro de cada grupo
sapply(unique(grupos.3),function(g)row.names(sul.p)[grupos.3 == g])


#Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(grupos.3),median)


# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.3),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.3),sd)
CV <- S/Xb*100
CV

# plotar os grupos a partir das variáveis originais
acp <- princomp(sul.p, cor = T)
plot(acp$scores[,1:2],type="n")
lab <- cutree(wr, k=3)
text(acp$scores[,1:2],labels = lab, cex = 1)
# outra opção
ggplot(acp$scores[,1:2], aes(acp$scores[,1],acp$scores[,2])) +  
  geom_text(aes(label = names(lab), color=lab), check_overlap = T, size = 3) + theme_bw() +
  scale_color_distiller(palette='Dark2') + 
  theme(legend.position="none")

library(RColorBrewer)
display.brewer.all()

# Método não hierárquico

# variáveis padronizadas - euclidiana - k-médias
set.seed(1)  # para que os resultados sempre sejam os mesmos
ng_km <- NbClust(sul.p, min.nc=2, max.nc=15, method="kmeans")
ng_km$Best.nc
table(ng_km$Best.nc[1,c(-23,-25)])
barplot(table(ng_km$Best.nc[1,c(-23,-25)]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos por 26 critérios")

# gamma
ng_gamma<- NbClust(sul.p, method="kmeans",index = "gamma")
ng_gamma$Best.nc

# gplus
ng_gplus<- NbClust(sul.p, method="kmeans", index = "gplus")
ng_gplus$Best.nc

# gap
ng_gap <- NbClust(sul.p, method="kmeans", index = "gap")
ng_gap$Best.nc

# 4 critérios: 2 grupos , maioria dos métodos 2 grupos

# ch:2; duda:2; cindex:5; gamma:11; beale:2; ccc:12; PtBiserial:4;
# gplus:12; DB:14; gap:2


# centers: número de grupos
# nstart: número de configurações iniciais diferentes
# k = 3, pois Ward resultou em 3 grupos ou k= 2 pelo NbClust?

km <- kmeans(sul.p, centers = 2, nstart = 25)   
km
km$centers  #centróides de cada grupo
km$cluster
km$withinss
km$size
# gráficos de dispersão
autoplot(kmeans(sul.p, 2, nstart=25), data = sul.p, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw() +ylab("componente principal 2")+
  xlab("componente principal 1")

set.seed(1)
autoplot(kmeans(sul.p, 3, nstart=25), data = sul.p, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")

# avaliar medidas estatísticas com o número k escolhido
set.seed(1)
grupos.3.K <- kmeans(sul.p, 3, nstart=25)$cluster
table(grupos.3.K)
sort(grupos.3.K)
aggregate(sul[,-c(1:4)],list(grupos.3.K),median)
# CV































op <- par(mfrow=c(1,2))
# 2 grupos
pam2 <- pam(sul[,-c(1:4)], 2)
plot(pam2)



# padronizar as variáveis - subtraindo a média e dividindo pelo
# desvio padrão
sul.p <- scale(sul[,-c(1:4)],center=T,scale=T)
sul.p <- as.data.frame(sul.p)
# sul.p é o data frame com as variáveis padronizadas

# transformar nomes dos municípios em row.names
row.names(sul.p) <- sul$nome.mun


# para obter o silhouette plot do agrupamento 
plot(silhouette(cutree(dm,2),dis))
plot(silhouette(cutree(dm,3),dis))

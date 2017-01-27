# Dissertação - Larissa ---------------------------------------------------
# Dados do Sul/Sudoeste de MG sobre envelhecimento
# AA usando variáveis originais e escores dos componentes principais

# diretório de trabalho --------------------------------------------------
setwd("/home/patricia/Dropbox/nupis/projetos/2015-larissa")

# pacotes -----------------------------------------------------------------

library(corrplot)  # Correlações
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

# library(rgl)          # Escores em 3 dimensões
# library(stargazer)    # Tabelas LaTEX

# dados Sul/Sudoeste de Minas -------------------------------------------------------------------
load("sul_dem.rda")

# agrupamento todas as variáveis ------------------------------------------

# padronizar as variáveis
sul.p <- scale(sul[,-c(1:4)],center=T,scale=T)
sul.p <- as.data.frame(sul.p)
# sul.p é o data frame com as variáveis padronizadas
# correlações
sul.p.r <- cor(sul.p)
round(sul.p.r,2)
corrplot(sul.p.r,tl.col = "black",type="lower")
corrplot.mixed(sul.p.r)

# transformar nomes dos municípios em row.names
row.names(sul.p) <- sul$nome.mun

# o mesmo para os dados originais (sul)
sul <- as.data.frame(sul)
row.names(sul) <- sul$nome.mun

# resumo estatístico
summary(sul[,-c(1:4)])
cor(sul[,-c(1:4)])

# correlações
sul.r <- cor(sul[,-c(1:4)])
round(sul.r,2)
corrplot(sul.r,tl.col = "black",type="lower")
corrplot.mixed(sul.r)

# testes de hipóteses
cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
res1 <- cor.mtest(sul[,-c(1:4)],0.95)
## specialized the insignificant value according to the significant level
corrplot(sul.r, p.mat = res1[[1]], sig.level=0.05)
## add all p-values
corrplot(sul.r, p.mat = res1[[1]], insig = "p-value", sig.level=-1)


# agrupamento de Variáveis - padronizadas
X.quanti <- sul.p
tree <- hclustvar(X.quanti)
plot(tree,main="",cex=1,ylab="altura")     # dendrograma todas as variáveis

# agrupamento de Variáveis - originais
X.quanti <- sul[,-c(1:4)]
str(X.quanti)
tree <- hclustvar(X.quanti)
plot(tree,main="",cex=1,ylab="altura")     # dendrograma todas as variáveis

sul.v <- sul.p


# variáveis originais -----------------------------------------------------
# agrupamento usando variáveis originais e distância de Mahalanobis

# matriz de distâncias - Mahalanobis
dis <- (distance(sul[,-c(1:4)],method="mahalanobis"))^0.5   #Mahalanobis

# distância média
dm <- hclust(dis, method = "average")
ggdendrogram(dm)
set.seed(1)
ng <- NbClust(sul[,-c(1:4)], min.nc=2, max.nc=20, method="average", 
              diss = dis, distance = NULL)
table(ng$Best.n[1,])
graphics.off()
barplot(table(ng$Best.n[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos pelos critérios")
# um dos grupos com observações demais

# centroide
ct <- hclust(dis, method = "centroid")
ggdendrogram(ct)
set.seed(1)
ng <- NbClust(sul[,-c(1:4)], min.nc=2, max.nc=20, method="centroid", 
              diss = dis, distance = NULL)
table(ng$Best.n[1,])
barplot(table(ng$Best.n[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos pelos critérios")
# um dos grupos com observações demais

# Ward
wr <- hclust(dis, method = "ward.D2")
ggdendrogram(wr)
set.seed(1)
ng <- NbClust(sul[,-c(1:4)], min.nc=2, max.nc=20, method="ward.D2", 
              diss = dis, distance = NULL)
table(ng$Best.n[1,])
graphics.off()
barplot(table(ng$Best.n[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos pelos critérios")
# resultados dos diferentes métodos de números de grupos
ng$All.index
ng$Best.nc
ng$All.CriticalValues
ng$Best.partition

# só GAP - deu erro
ng <- NbClust(sul[,-c(1:4)], min.nc=2, max.nc=10, method="ward.D2", 
              diss = dis, distance = NULL, index = "gap")

# avaliar os 2 grupos
grupos.2 <- cutree(wr, k=2); table(grupos.2)
# Munic?pios dentro de cada grupo
sapply(unique(grupos.2),function(g)row.names(sul[,-c(1:4)])[grupos.2 == g])
# Avaliar medidas estatisticas - vari?veis originais - 
aggregate(sul[,-c(1:4)],list(grupos.2),mean)
aggregate(sul[,-c(1:4)],list(grupos.2),median)
aggregate(sul[,-c(1:4)],list(grupos.2),sd)
aggregate(sul[,-c(1:4)],list(grupos.2),max)
aggregate(sul[,-c(1:4)],list(grupos.2),min)
# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.2),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.2),sd)
CV <- S/Xb*100
CV

# http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
# pacote cluster
op <- par(mfrow=c(1,2))
# 2 grupos
pam2 <- pam(sul[,-c(1:4)], 2)
plot(pam2)
pamk.best <- pamk(sul[,-c(1:4)])
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(sul[,-c(1:4)], pamk.best$nc))  # best: 2 grupos
plot(pam(sul[,-c(1:4)], 2))  # best: 2 grupos
plot(pam(sul[,-c(1:4)], 3))  # com 3
plot(pam(sul[,-c(1:4)], 4))  # com 4

# http://strata.uga.edu/software/pdf/clusterTutorial.pdf
library(cluster)     
mydata.agnes <- agnes(sul[,-c(1:4)])
plot(mydata.agnes, which.plots = 2, main= "Dendrogram of my data")
# manhattan, Ward
mydata.agnes.ALT <- agnes(sul[,-c(1:4)], metric="manhattan", method="ward")
plot(mydata.agnes.ALT, which.plots = 2, main= "Dendrogram of my data")
# euclidiana, Ward
mydata.agnes.ALT <- agnes(sul[,-c(1:4)], metric="euclidian", method="ward")
plot(mydata.agnes.ALT, which.plots = 2, main= "Dendrogram of my data")
# mahalanobis, ward
dis.m <- (distance(sul[,-c(1:4)],method="mahalanobis"))^0.5   #Mahalanobis
agnes.mw <- agnes(dis.m, diss = F, method="ward")
plot(agnes.mw, which.plots = 4, main= "Dendrogram of my data")


# Variáveis padronizadas --------------------------------------------------
# usando dados padronizados e distância euclidiana
dis <- dist(sul.p)

# distância média
dm <- hclust(dis, method = "average")
ggdendrogram(dm)
set.seed(1)
ng <- NbClust(sul.p, min.nc=2, max.nc=20, method="average")
table(ng$Best.n[1,])
barplot(table(ng$Best.n[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos pelos critérios")
# um dos grupos com observações demais
# avaliar os 2 grupos
grupos.2 <- cutree(dm, k=2); table(grupos.2)
# Munic?pios dentro de cada grupo
sapply(unique(grupos.2),function(g)row.names(sul[,-c(1:4)])[grupos.2 == g])
# Avaliar medidas estatisticas - vari?veis originais - 
aggregate(sul[,-c(1:4)],list(grupos.2),mean)
aggregate(sul[,-c(1:4)],list(grupos.2),median)
aggregate(sul[,-c(1:4)],list(grupos.2),sd)
aggregate(sul[,-c(1:4)],list(grupos.2),max)
aggregate(sul[,-c(1:4)],list(grupos.2),min)
# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.2),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.2),sd)
CV <- S/Xb*100
CV

# centroide
ct <- hclust(dis, method = "centroid")
ggdendrogram(ct)
set.seed(1)
ng <- NbClust(sul.p, min.nc=2, max.nc=20, method="centroid")
table(ng$Best.n[1,])
barplot(table(ng$Best.n[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos pelos critérios")
# um dos grupos com observações demais
# avaliar os 2 grupos
grupos.2 <- cutree(ct, k=2); table(grupos.2)
# Munic?pios dentro de cada grupo
sapply(unique(grupos.2),function(g)row.names(sul[,-c(1:4)])[grupos.2 == g])
# Avaliar medidas estatisticas - vari?veis originais - 
aggregate(sul[,-c(1:4)],list(grupos.2),mean)
aggregate(sul[,-c(1:4)],list(grupos.2),median)
aggregate(sul[,-c(1:4)],list(grupos.2),sd)
aggregate(sul[,-c(1:4)],list(grupos.2),max)
aggregate(sul[,-c(1:4)],list(grupos.2),min)
# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.2),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.2),sd)
CV <- S/Xb*100
CV

# v) Ward
wr <- hclust(dis, method = "ward.D2")
ggdendrogram(wr)
set.seed(1)
ng <- NbClust(sul.p, min.nc=2, max.nc=10, method="ward.D2")
table(ng$Best.n[1,])
barplot(table(ng$Best.n[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos pelos critérios")
# só gap
ng <- NbClust(sul.p, min.nc=2, max.nc=10, method="ward.D2", index = "gap")
ng$Best.nc
ng$All.CriticalValues
ng$Best.partition

# avaliar os 3 grupos
grupos.3 <- cutree(wr, k=3); table(grupos.3)
# Munic?pios dentro de cada grupo
sapply(unique(grupos.3),function(g)row.names(sul[,-c(1:4)])[grupos.3 == g])
# Avaliar medidas estatisticas - vari?veis originais - 
aggregate(sul[,-c(1:4)],list(grupos.3),mean)
aggregate(sul[,-c(1:4)],list(grupos.3),median)
aggregate(sul[,-c(1:4)],list(grupos.3),sd)
aggregate(sul[,-c(1:4)],list(grupos.3),max)
aggregate(sul[,-c(1:4)],list(grupos.3),min)
# avaliar o CV de cada grupo (quanto menor, mais homogeneo)
Xb <- aggregate(sul[,-c(1:4)],list(grupos.3),mean)
S <- aggregate(sul[,-c(1:4)],list(grupos.3),sd)
CV <- S/Xb*100
CV


# pam
op <- par(mfrow=c(1,2))
# 2 grupos
pam2 <- pam(sul.p, 2)
plot(pam2)
pamk.best <- pamk(sul.p)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(sul.p, pamk.best$nc))  # best: 2 grupos
plot(pam(sul.p, 3))  # com 3
plot(pam(sul.p, 4))  # com 4

library(fpc)
asw <- numeric(20)
for (k in 2:20)
  asw[[k]] <- pam(sul.p, k)$silinfo$avg.width
k.best <- which.max(asw)
cat("silhouette-optimal number of clusters:", k.best, "\n")

# diferente padronização --------------------------------------------------
# usando dados padronizados pela amplitude e distância euclidiana
sul.v <- sul
names(sul.v)
summary(sul.v)
sul.v$espvida <- with(sul.v, max(sul.v$espvida)-sul.v$espvida)
sul.v$fectot <- with(sul.v, max(sul.v$fectot)-sul.v$fectot)
sul.v$mort1 <- with(sul.v, max(sul.v$mort1)-sul.v$mort1)
sul.v$mort5 <- with(sul.v, max(sul.v$mort5)-sul.v$mort5)
sul.v$razdep <- with(sul.v, max(sul.v$razdep)-sul.v$razdep)
sul.v$sobre40 <- with(sul.v, max(sul.v$sobre40)-sul.v$sobre40)
sul.v$sobre60 <- with(sul.v, max(sul.v$sobre60)-sul.v$sobre60)
sul.v$t_env <- with(sul.v, max(sul.v$t_env)-sul.v$t_env)
sul.v$pop <- scale(sul.p$pop,center=T,scale=T)

sul.v <- as.data.frame(sul.v)
row.names(sul.v) <- sul.v$nome.mun
dis <- dist(sul.v[,-c(1:4)])

# distância média
dm <- hclust(dis, method = "average")
ggdendrogram(dm)
set.seed(1)
ng <- NbClust(sul.p, min.nc=2, max.nc=20, method="average")
table(ng$Best.n[1,])
barplot(table(ng$Best.n[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos pelos critérios")
# um dos grupos com observações demais

# centroide
ct <- hclust(dis, method = "centroid")
ggdendrogram(ct)
set.seed(1)
ng <- NbClust(sul.v, min.nc=2, max.nc=20, method="centroid")
table(ng$Best.n[1,])
barplot(table(ng$Best.n[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos pelos critérios")
# um dos grupos com observações demais

# v) Ward
wr <- hclust(dis, method = "ward.D2")
library(ggdendro)
ggdendrogram(wr)
set.seed(1)
ng <- NbClust(sul.v[,-c(1:4)], min.nc=2, max.nc=10, method="ward.D2")
table(ng$Best.n[1,])
barplot(table(ng$Best.n[1,]),
        xlab="Número de grupos", ylab="Número de critérios",
        main="Número de grupos pelos critérios")
# só gap
ng <- NbClust(sul.p, min.nc=2, max.nc=10, method="ward.D2", 
              diss = dis, distance = NULL, index = "gap")
ng$Best.nc
ng$All.CriticalValues
ng$Best.partition

# pam
op <- par(mfrow=c(1,2))
# 2 grupos
pam2 <- pam(sul.v[,-c(1:4)], 2)
plot(pam2)
pamk.best <- pamk(sul.v[,-c(1:4)])
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(sul.v[,-c(1:4)], pamk.best$nc))  # best: 2 grupos
plot(pam(sul.v[,-c(1:4)], 3))  # com 3
plot(pam(sul.v[,-c(1:4)], 4))  # com 4

library(fpc)
asw <- numeric(20)
for (k in 2:20)
  asw[[k]] <- pam(sul.v[,-c(1:4)], k)$silinfo$avg.width
k.best <- which.max(asw)
cat("silhouette-optimal number of clusters:", k.best, "\n")

# pacote clusterSim
# install.packages("clusterSim")
library(clusterSim)
###   VER DEPOIS

# Análise de Componentes Principais - ACP ---------------------------------------
# ACP com dados originais e usando matriz de correlações R
acp1 <- prcomp(sul[,-c(1:4)],scale=T)
summary(acp1)  # 2 CPs
acp1$rotation[,1:3]      # escores CPs
# scree plot
lamb <- acp1$sdev^2
lambp <- lamb/sum(lamb)
screep <- data.frame(lamb,lambp) 
ggplot(data=screep,aes(x=1:9,y=lambp))+geom_line()+theme_bw()+
  ylab("% da variação explicada")+xlab("componentes principais")
# padrão é usar a matriz R e não a S
acp11 <- PCA(sul[,-c(1:4)],graph=T)
plot(acp11,title=" ",choix="var",ylab = "componente principal 1",
     xlab = "componente principal 2")
summary(acp11)
acp11$eig        #autovalores
# correlações entre CPs e X
acp11$var$cor[,1:2]
# escores
acp11$ind$coord
# Criar data frame com escores
scores <- as.data.frame(acp11$ind$coord[,1:2])
# Diagrama de dispersão dos municípios de acordo com os escores do CP1 e do CP2
ggplot(data=scores,aes(x=Dim.1,y=Dim.2,label=rownames(scores)))+
  geom_hline(yintercept=0,colour="gray65")+
  geom_vline(xintercept=0,colour="gray65")+
  geom_text(colour="black",alpha=0.8,size=4)+
  ggtitle(" ") + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")

# ACP com dados padronizados e usando matriz de covariâncias S 
# parece o mesmo resultado que usar os dados originais e utilizar R
acp2 <- prcomp(sul.p,scale=F)
summary(acp2)  # 2 CPs
acp2$rotation[,1:2]      # escores CPs
# scree plot
lamb <- acp2$sdev^2
lambp <- lamb/sum(lamb)
screep <- data.frame(lamb,lambp) 
ggplot(data=screep,aes(x=1:9,y=lambp))+geom_line()+theme_bw()+
  ylab("% da variação explicada")+xlab("componentes principais")
# padrão é usar a matriz R e não a S
acp21 <- PCA(sul.p,scale.unit = F, graph=T)
plot(acp21,title=" ",choix="var",ylab = "componente principal 1",xlab = "componente principal 2")
summary(acp21)
acp21$eig        #autovalores
# correlações entre CPs e X
acp21$var$cor[,1:2]
# escores
acp21$ind$coord
# Criar data frame com escores
scores <- as.data.frame(acp21$ind$coord[,1:2])
# Diagrama de dispersão dos municípios de acordo com os escores do CP1 e do CP2
ggplot(data=scores,aes(x=Dim.1,y=Dim.2,label=rownames(scores)))+
  geom_hline(yintercept=0,colour="gray65")+
  geom_vline(xintercept=0,colour="gray65")+
  geom_text(colour="black",alpha=0.8,size=4)+
  ggtitle(" ") + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")

# ACP com dados padronizados e usando matriz R
acp3 <- prcomp(sul.p,scale=T)
summary(acp3)
acp3$rotation[,1:2]      #Escores CPs
#Scree plot
lamb <- acp3$sdev^2
lambp <- lamb/sum(lamb)
screep <- data.frame(lamb,lambp) 
ggplot(data=screep,aes(x=1:9,y=lambp))+geom_line()+theme_bw()+
  ylab("% da variação explicada")+xlab("componentes principais")
acp31 <- PCA(sul.p,graph=T)
plot(acp31,title=" ",choix="var",ylab = "componente principal 1",xlab = "componente principal 2")
summary(acp31)
acp31$eig        #autovalores
#correlações entre CPs e X
acp31$var$cor[,1:2]
# Escores
acp31$ind$coord
# Criar data frame com escores
scores <- as.data.frame(acp$ind$coord[,1:2])
# Diagrama de dispersão dos municípios de acordo com os escores do CP1 e do CP2
ggplot(data=scores,aes(x=Dim.1,y=Dim.2,label=rownames(scores)))+
  geom_hline(yintercept=0,colour="gray65")+
  geom_vline(xintercept=0,colour="gray65")+
  geom_text(colour="black",alpha=0.8,size=4)+
  ggtitle(" ") + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")

# data frame com os escores dos 2 CPs
scores2 <- scores

# Análise de agrupamento - Método de Ward - Dist. Euclidiana ---------------------------------------
dis <- dist(scores2)     #Distância eucliadiana
wr <- hclust(dis, method = "ward.D2")

grupos.2 <- cutree(wr, k=2)  #2 grupos
grupos.3 <- cutree(wr, k=3)  #3 grupos
#Recuperar as observações
#Todas de uma vez
sapply(unique(grupos.2),function(g)row.names(sul.p)[grupos.2 == g])
sapply(unique(grupos.3),function(g)row.names(sul.p)[grupos.3 == g])


#Dendrograma
hc <- dis %>% hclust(method = "ward.D2") 
ddata <- hc %>% as.dendrogram() %>% dendro_data()
ggdendrogram(hc) + geom_text(size= 2.5, aes(x = x, y = y, label = label, 
                                            angle = -90, vjust=0, hjust = 0), data= label(ddata)) +
  scale_y_continuous(expand = c(0.6, 0)) +
  scale_x_continuous(expand = c(0.1, 0)) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) 


#NbClust - Número de grupos ideal segundo vários critérios
set.seed(1)
ng <- NbClust(scores2, min.nc=2, max.nc=15, method="ward.D2")

#Criar vetor para plotar em GGplot2
freqg <- as.numeric(table(ng$Best.n[1,]))
ngrupo <- as.numeric(names(table(ng$Best.n[1,])))
ngruposW <- data.frame(ngrupo,freqg)
vetor <- c()
for(i in 1:nrow(ngruposW)){
  vetor <- c(vetor,rep(ngruposW[i,1],ngruposW[i,2]))
}
vetor <- data.frame(vetor)
ggplot(data=vetor,aes(x=vetor))+geom_bar()+theme_bw()+
  ylab("número de critérios")+xlab("número de grupos")

#Avaliar medidas estatisticas - variáveis originais - interessante
aggregate(sul.v,list(grupos.3),mean)
aggregate(sul.v,list(grupos.3),median)
aggregate(sul.v,list(grupos.3),sd)

#Obter gráfico escores colorido método de Ward
set.seed(1)
teste <-  kmeans(scores2, 3, nstart=25)

teste$cluster <- grupos.3
names

autoplot(teste, data = scores2, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")


# Análise de agrupamento - K-médias ---------------------------------------

#número de grupos definido pelo método hierárquico aglomerativo de Ward

#agrupamentos
set.seed(1)
autoplot(kmeans(scores2, 2, nstart=25), data = scores2, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw() +ylab("componente principal 2")+
  xlab("componente principal 1")
set.seed(1)
autoplot(kmeans(scores2, 3, nstart=25), data = scores2, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")


#Avaliar medidas estatísticas
set.seed(2)
grupos.3.K <- kmeans(scores2, 3, nstart=25)$cluster
aggregate(sul.v,list(grupos.3.K),mean)
aggregate(sul.v,list(grupos.3.K),sd)


# Data frame agrupamentos -------------------------------------------------

#Ward
x1 <- as.data.frame(grupos.2)
y1 <- as.data.frame(grupos.3)

#K-médias
set.seed(8)
grupos.2K <- kmeans(scores2, 2, nstart=25)$cluster
set.seed(8)
grupos.3K <- kmeans(scores2, 3, nstart=25)$cluster
x2 <- as.data.frame(grupos.2K)
y2 <- as.data.frame(grupos.3K)

x <- data.frame(nome.mun=row.names(x1), dois.W=as.factor(x1$grupos.2), 
                tres.W=as.factor(y1$grupos.3),dois.K=as.factor(x2$grupos.2K), 
                tres.K=as.factor(y2$grupos.3K))

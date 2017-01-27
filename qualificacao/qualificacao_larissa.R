# Dissertação - Larissa ---------------------------------------------------
# Dados do Sul/Sudoeste de MG sobre envelhecimento
# AA usando variáveis originais e escores dos componentes principais

# diretório de trabalho --------------------------------------------------
setwd("/home/patricia/Dropbox/nupis/projetos/2015-aa-larissa")

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
library(rgl)          # Escores em 3 dimensões
library(stargazer)    # Tabelas LaTEX

# dados Sul Sudoeste de Minas -------------------------------------------------------------------
load("sul_dem.rda")

# agrupamento todas as variáveis ------------------------------------------

sul.p <- scale(sul[,-c(1:4)],center=T,scale=T)
sul.p <- as.data.frame(sul.p)

# transformar nomes dos municípios em row.names
row.names(sul.p) <- sul$nome.mun

# resumo estatístico
summary(sul.p)
cor(sul.p)

# correlações
sul.r <- cor(sul.p)
round(sul.r,2)
corrplot(sul.r,tl.col = "black",type="lower")

# agrupamento de Variáveis 
X.quanti <- sul.p
tree <- hclustvar(X.quanti)
plot(tree,main="",cex=1,ylab="altura")     # dendrograma todas as variáveis

sul.v <- sul.p

# Variáveis selecionadas --------------------------------------------------
# sul.v <- sul.v[,c(1,6,7,8,9)]

# variáveis padronizadas
#sul.p <- scale(sul.v,center=T,scale=T)
#sul.p <- as.data.frame(sul.p)
#names(sul.p)
#head(sul.p)

#transformar nomes dos municípios em row.names
#row.names(sul.p) <- sul$nome.mun
#head(sul.p)
#names(sul.p)

# análise de Componentes Principais - ACP ---------------------------------------

acp1 <- prcomp(sul.p,scale=T)
summary(acp1)
acp1$rotation[,1:2]      # escores CPs

# scree plot
lamb <- acp1$sdev^2
lambp <- lamb/sum(lamb)
screep <- data.frame(lamb,lambp) 
ggplot(data=screep,aes(x=1:9,y=lambp))+geom_line()+theme_bw()+
  ylab("% da variação explicada")+xlab("componentes principais")

# padrão é usar a matriz R e não a S
acp <- PCA(sul.p,graph=T)
plot(acp,title=" ",choix="var",ylab = "componente principal 1",
     xlab = "componente principal 2")
summary(acp)
acp$eig        # autovalores

# correlações entre CPs e X
acp$var$cor[,1:2]

# escores
acp$ind$coord

# criar data frame com escores
scores <- as.data.frame(acp$ind$coord)

#dDiagrama de dispersão dos municípios de acordo com os escores do CP1 e do CP2
ggplot(data=scores,aes(x=Dim.1,y=Dim.2,label=rownames(scores)))+
  geom_hline(yintercept=0,colour="gray65")+
  geom_vline(xintercept=0,colour="gray65")+
  geom_text(colour="black",alpha=0.8,size=4)+
  ggtitle(" ") + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")

# data frame com os escores dos 3 CPs
scores2 <- scores[,1:2]

# análise de agrupamento - Método de Ward - Dist. Euclidiana ---------------------------------------

dis <- dist(scores2)     # distância eucliadiana
wr <- hclust(dis, method = "ward.D2")
ggdendrogram(wr)

grupos.2 <- cutree(wr, k=2)  # 2 grupos
grupos.3 <- cutree(wr, k=3)  # 3 grupos
# Recuperar as observações
# Todas de uma vez
sapply(unique(grupos.2),function(g)row.names(sul.p)[grupos.2 == g])
sapply(unique(grupos.3),function(g)row.names(sul.p)[grupos.3 == g])

# dendrograma
hc <- dis %>% hclust(method = "ward.D2") 
ddata <- hc %>% as.dendrogram() %>% dendro_data()
ggdendrogram(hc) + geom_text(size= 2.5, aes(x = x, y = y, label = label, 
                                            angle = -90, vjust=0, hjust = 0), data= label(ddata)) +
  scale_y_continuous(expand = c(0.6, 0)) +
  scale_x_continuous(expand = c(0.1, 0)) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) 

# NbClust - Número de grupos ideal segundo vários critérios
set.seed(1)
ng <- NbClust(scores2, min.nc=2, max.nc=15, method="ward.D2")

# criar vetor para plotar em GGplot2
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
# a maioria dos critérios sugeriu ?

# avaliar medidas estatisticas - variáveis originais - interessante
aggregate(sul.v,list(grupos.3),mean)
aggregate(sul.v,list(grupos.3),median)
aggregate(sul.v,list(grupos.3),sd)


#Obter gráfico escores colorido método de Ward
set.seed(8)
teste <-  kmeans(scores3, 3, nstart=25)

teste$cluster <- grupos.3
names

autoplot(teste, data = scores3, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")


# Análise de agrupamento - K-médias ---------------------------------------

#número de grupos definido pelo método hierárquico aglomerativo de Ward

#agrupamentos
set.seed(8)
autoplot(kmeans(scores3, 2, nstart=25), data = scores3, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw() +ylab("componente principal 2")+
  xlab("componente principal 1")
set.seed(8)
autoplot(kmeans(scores3, 3, nstart=25), data = scores3, label = TRUE, shape = FALSE,
         label.size = 4) + theme_bw()+ylab("componente principal 2")+
  xlab("componente principal 1")


#Avaliar medidas estatísticas
set.seed(8)
grupos.3.K <- kmeans(scores3, 3, nstart=25)$cluster
aggregate(sul.v,list(grupos.3.K),mean)
aggregate(sul.v,list(grupos.3.K),sd)


# Data frame agrupamentos -------------------------------------------------

#Ward
x1 <- as.data.frame(grupos.2)
y1 <- as.data.frame(grupos.3)

#K-médias
set.seed(8)
grupos.2K <- kmeans(scores3, 2, nstart=25)$cluster
set.seed(8)
grupos.3K <- kmeans(scores3, 3, nstart=25)$cluster
x2 <- as.data.frame(grupos.2K)
y2 <- as.data.frame(grupos.3K)

x <- data.frame(nome.mun=row.names(x1), dois.W=as.factor(x1$grupos.2), 
                tres.W=as.factor(y1$grupos.3),dois.K=as.factor(x2$grupos.2K), 
                tres.K=as.factor(y2$grupos.3K))


# Mapas -------------------------------------------------------------------

dd <- filter(atlas, meso==3110)
dd$nome.mun <- as.factor(dd$nome.mun)
x$nome.mun <- as.factor(x$nome.mun)
xx <- inner_join(dd, x, by="nome.mun")
names(xx)

#mgm <-readOGR(dsn="/home/lincoln/Dropbox/dados/malhas",layer="31mu2500gsr")
mgm <-readOGR(dsn="C:/Users/Asus/Dropbox/Unifal/NUPIS/Mapas/31mu2500gsr",layer="31mu2500gsr")
names(mgm)
mgm1<-subset(mgm, MESORREGIÃ.==3110)
mgm1@data<-rename(mgm1@data, codmun7=GEOCODIG_M)
mgm1$codmun7<-as.factor(mgm1$codmun7)
mga <- inner_join(mgm1@data, xx, by="codmun7")  
mgf<-fortify(mgm1, region="codmun7") # demora um pouco
mga$id<-mga$codmun7
mgf <-inner_join(mgf, mga, by = "id")
mapa <- ggplot(mgf, aes(long, lat, group=group, fill = tres.W)) + 
  geom_polygon() + coord_equal() + 
  theme_void() + scale_fill_brewer(palette = "RdYlGn")
mapa


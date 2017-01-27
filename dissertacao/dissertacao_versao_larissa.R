# Dados do Sul/Sudoeste de MG sobre envelhecimento
# AA usando variáveis originais e distância de Mahalanobis
# métodos: Ward e k-médias

# fase preliminar ---------------------------------------------------------
# diretório de trabalho
setwd("/home/patricia/Dropbox/nupis/projetos/2015-larissa")
# dados Sul/Sudoeste de Minas
load("sul_dem.rda")

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
# olhar separadamente para as variávies
names(sul[,-c(1:4)])
boxplot(sul$espvida)
boxplot(sul$tft)
sort(sul$tft)
sul$nome.mun[which.min(sul$tft)]
sul$nome.mun[which.max(sul$tft)]
boxplot(sul$mort1)
boxplot(sul$mort5)
boxplot(sul$rd)
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
escores=as.data.frame(escores)
# usando ggplot2
graf_CP <- ggplot(data=escores,aes(x=escores[,1],y=escores[,2],label=rownames(escores)))
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour="black",alpha=0.8,size=4) +
  ggtitle(" ") + theme_few()+ylab("componente principal 2") +
  xlab("componente principal 1")


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

#CPS

lab4_wr <- cutree(wr, k=4)
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=lab4_wr,alpha=0.8,size=4) +
  ggtitle(" ") + theme_bw() + ylab("componente principal 2") +
  xlab("componente principal 1")


require(gridExtra)


# k-médias - k=2, k=3, k=4, k=5 ----------------------------------
# k = 4
set.seed(1)
k_4 <- kmeans(sul[,-c(1:4)], 4, nstart=25)
labk4 <- k_4$cluster
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour=labk4,alpha=0.8,size=4) +
  ggtitle(" ") + theme_few() + ylab("componente principal 2") +
  xlab("componente principal 1")

# avaliar grupos obtidos
table(labk4)
# Municípios dentro de cada grupo
sapply(unique(labk4),function(g)row.names(sul[,-c(1:4)])[labk4 == g])
# Avaliar medidas estatisticas - variáveis originais - 
aggregate(sul[,-c(1:4)],list(labk4),mean)
aggregate(sul[,-c(1:4)],list(labk4),median)
aggregate(sul[,-c(1:4)],list(labk4),min)
aggregate(sul[,-c(1:4)],list(labk4),max)


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
#a <- ggplot(data = sul4_km, aes(as.factor(grupo), espvida)) + 
 # geom_boxplot(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge", outlier.colour = "red", 
  #             outlier.shape = 16, outlier.size = 4, notch = FALSE, notchwidth = 0.5) + 
  #scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida))) 
#a + geom_text(label=sul4_km$nome.mun) 

a<- ggplot(data = sul4_km, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida)))

# fectot
b <- ggplot(data = sul4_km, aes(as.factor(grupo), tft)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$tft),max(sul$tft)))
# mort1
c <- ggplot(data = sul4_km, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1)))
# razdep
d <- ggplot(data = sul4_km, aes(as.factor(grupo), rd)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$rd),max(sul$rd)))
# sobre60
e <- ggplot(data = sul4_km, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60)))
# t_env
f <- ggplot(data = sul4_km, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env)))

grid.arrange(a,b,c,d,e,f,ncol=3)



########## OUTLIERS

sul4_km <- sul
sul4_km$grupo <- labk4
a <- ggplot(data = sul4_km, aes(as.factor(grupo), espvida)) + 
  geom_boxplot(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge", outlier.colour = "red", 
               outlier.shape = 16, outlier.size = 4, notch = FALSE, notchwidth = 0.5) + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida))) 
a + geom_text(label=sul4_km$nome.mun) 

# fectot
b <- ggplot(data = sul4_km, aes(as.factor(grupo), tft)) + 
  geom_boxplot(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge", outlier.colour = "red", 
               outlier.shape = 16, outlier.size = 4, notch = FALSE, notchwidth = 0.5) + 
  scale_y_continuous(limits = c(min(sul$tft),max(sul$tft))) 
b + geom_text(label=sul4_km$nome.mun)


# mort1
c <- ggplot(data = sul4_km, aes(as.factor(grupo), mort1)) + 
  geom_boxplot(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge", outlier.colour = "red", 
               outlier.shape = 16, outlier.size = 4, notch = FALSE, notchwidth = 0.5) + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1))) 
c + geom_text(label=sul4_km$nome.mun)


# razdep

d <- ggplot(data = sul4_km, aes(as.factor(grupo), rd)) + 
  geom_boxplot(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge", outlier.colour = "red", 
               outlier.shape = 16, outlier.size = 4, notch = FALSE, notchwidth = 0.5) + 
  scale_y_continuous(limits = c(min(sul$rd),max(sul$rd))) 
d + geom_text(label=sul4_km$nome.mun)

# sobre60
e <- ggplot(data = sul4_km, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge", outlier.colour = "red", 
               outlier.shape = 16, outlier.size = 4, notch = FALSE, notchwidth = 0.5) + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60))) 
e + geom_text(label=sul4_km$nome.mun)


# t_env

f <- ggplot(data = sul4_km, aes(as.factor(grupo), t_env)) + 
  geom_boxplot(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge", outlier.colour = "red", 
               outlier.shape = 16, outlier.size = 4, notch = FALSE, notchwidth = 0.5) + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env))) 
f + geom_text(label=sul4_km$nome.mun)


grid.arrange(a,b,c,d,e,f,ncol=3)


############## população e renda da população ocupada ###############

dados_n <- read.csv("teste7.csv",h=T,sep=";")
save(dados_n, file = "sul_dem_n.Rda")
load("sul_dem_n.Rda")
# nomes como rótulos
row.names(dados_n) <- dados_n$nome_mun
# tirar nome_mun
dados_n <- dados_n[,-1]
save(dados_n, file = "sul_dem_n.Rda")
load("sul_dem_n.Rda")
# resumo estatístico
aggregate(dados_n,list(labk4),median)
aggregate(dados_n,list(labk4),min)
aggregate(dados_n,list(labk4),max)
# CV
Xb <- aggregate(dados_n,list(labk4),mean)
S <- aggregate(dados_n,list(labk4),sd)
CV <- S/Xb*100
CV

# boxplots
# k-médias, k = 4
# espvida
dados_n$grupo <- labk4
a <- ggplot(data = dados_n, aes(as.factor(grupo), pop)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(dados_n$pop),max(dados_n$pop)))
# fectot
b <- ggplot(data = dados_n, aes(as.factor(grupo), renocup)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(dados_n$renocup),max(dados_n$renocup)))
grid.arrange(a,b,ncol=2)
# mort1
c <- ggplot(data = dados_n, aes(as.factor(grupo), e_anosestudo)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(dados_n$e_anosestudo),max(dados_n$e_anosestudo)))
grid.arrange(a,b,c,ncol=3)





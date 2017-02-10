# Dados do Sul/Sudoeste de MG sobre envelhecimento
# AA usando variáveis originais e distância de Mahalanobis
# métodos: Ward e k-médias

# fase preliminar ---------------------------------------------------------
# carregar dados Sul/Sudoeste de Minas
load("/home/patricia/Dropbox/nupis/projetos/2015-larissa/sul_dem.rda")

# pacotes utilizados
library(corrplot)     # correlações
library(ecodist)      # Ward - distância de Mahalanobis
library(ggplot2)      # gráficos
library(ggthemes)     # gráficos
library(ggdendro)     # dendrogramas
library(dplyr)        # manipulação de dados
library(gridExtra)    # apresentação de gráficos combinados

# limpar dados 
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
# correlações
graphics.off()
sul.r <- cor(sul[,-c(1:4)])
corrplot.mixed(sul.r,tl.col = "black",
               number.digits = 4, number.font = 0.1,
               tl.cex = 1, is.corr = F)
# retirar as variáveis correlacionadas com outras (mort5 e sobre40)
# e que, de alguma forma, calculam a mesma coisa
sul <- sul[,-c(8,10)]
# correlações após retirar as variáveis
graphics.off()
sul.r <- cor(sul[,-c(1:4)])
corrplot.mixed(sul.r,tl.col = "black",
               number.digits = 4, number.font = 0.1,
               tl.cex = 1, is.corr = F)

# corrplot com valores em preto
corrplot(sul.r, add=T, type="lower", method="number",
         col="black", diag=F, tl.pos="n", cl.pos="n", 
         number.digits = 4, number.font = 0.5)

# variáveis originais - Mahalanobis ---------------------------------------
# matriz de distâncias - Mahalanobis
dis <- (distance(sul[,-c(1:4)],method="mahalanobis"))^0.5
# agrupamento hierárquico de Ward
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
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept = 13, color = 'red') +
  geom_hline(yintercept = 11.5, color = 'gray') +
  geom_hline(yintercept = 8.8, color = 'gray') +
  annotate("text", label = "4 grupos", x = 135, y = 13.5, size = 3, colour = "red") +
  annotate("text", label = "5 grupos", x = 135, y = 12, size = 3, colour = "gray") +
  annotate("text", label = "7 grupos", x = 135, y = 9.3, size = 3, colour = "gray")

# gráficos - dois primeiros CPs - Mahalanobis
X_cp <- princomp(sul[,-c(1:4)], cor=T)
summary(X_cp) # os dois primeiros contabilizam 76% da variância total
escores <- X_cp$scores[,1:2]
escores <- as.data.frame(escores)
graf_CP <- ggplot(data=escores,aes(x=escores[,1],y=escores[,2],label=rownames(escores)))
graf_CP +  geom_hline(yintercept=0,colour="gray65") +
  geom_vline(xintercept=0,colour="gray65") +
  geom_text(colour="black",alpha=0.8,size=4) +
  ggtitle(" ") + theme_few() + ylab("componente principal 2") +
  xlab("componente principal 1")

# k-médias - k = 4
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
# municípios dentro de cada grupo
sapply(unique(labk4),function(g)row.names(sul[,-c(1:4)])[labk4 == g])
# avaliar medidas estatisticas - variáveis originais - 
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
a <- ggplot(data = sul4_km, aes(as.factor(grupo), espvida)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$espvida),max(sul$espvida))) +
  xlab("grupo") + ylab("esperança de vida")

# tft
b <- ggplot(data = sul4_km, aes(as.factor(grupo), tft)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$tft),max(sul$tft))) +
  xlab("grupo") + ylab("taxa de fecundidade total")

# mort1
c <- ggplot(data = sul4_km, aes(as.factor(grupo), mort1)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$mort1),max(sul$mort1))) +
  xlab("grupo") + ylab("mortalidade até 1 ano")

# rd
d <- ggplot(data = sul4_km, aes(as.factor(grupo), rd)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$rd),max(sul$rd))) +
  xlab("grupo") + ylab("razão de dependência")

# sobre60
e <- ggplot(data = sul4_km, aes(as.factor(grupo), sobre60)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$sobre60),max(sul$sobre60))) +
  xlab("grupo") + ylab("sobrevivência até 60 anos")

# t_env
f <- ggplot(data = sul4_km, aes(as.factor(grupo), t_env)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(sul$t_env),max(sul$t_env))) +
  xlab("grupo") + ylab("taxa de envelhecimento")

# combinando os gráficos
grid.arrange(a,b,c,d,e,f,ncol=3)

# outras variáveis não usadas no agrupamento:
# pop, renocup e e_anosestudo
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
# pop
dados_n$grupo <- labk4
a <- ggplot(data = dados_n, aes(as.factor(grupo), pop)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(dados_n$pop),max(dados_n$pop))) +
  xlab("grupo") + ylab("população")

# renocup
b <- ggplot(data = dados_n, aes(as.factor(grupo), renocup)) + 
  geom_boxplot() + 
  scale_y_continuous(limits = c(min(dados_n$renocup),max(dados_n$renocup))) +
  xlab("grupo") + ylab("renda dos ocupados")

grid.arrange(a,b,ncol=2)

# mapas
library(maptools) 
library(rgdal)
library(dplyr)
library(rgeos)

# criar data frame com nome do município e a qual grupo pertence
z <- as.data.frame(labk4)
x <- data.frame(nome.mun=row.names(z), grupo=as.factor(z$labk4))
x
# incluir as outras informações dos municípios
sul$nome.mun <- as.factor(sul$nome.mun)
x$nome.mun <- as.factor(x$nome.mun)
xx <- inner_join(sul, x, by="nome.mun")

# faça o download dos arquivos necessários para criar o mapa em:
# ftp://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2015/UFs/MG/mg_municipios.zip
# descompacte os arquivos em uma pasta e ajuste o endereço abaixo:
mgm <- readOGR(dsn="/home/patricia/Dropbox/nupis/projetos/2015-larissa/dissertacao/mapas",
               layer="31MUE250GC_SIR") 
mgm@data <- rename(mgm@data, codmun7=CD_GEOCMU)
mgm$codmun7 <- as.factor(mgm$codmun7)
mga <- right_join(mgm@data, xx, by="codmun7")  
mgf <- fortify(mgm, region="codmun7")

mga$id <- mga$codmun7
mgf <- inner_join(mgf, mga, by = "id")
ggplot(mgf, aes(long, lat, group=group, fill=grupo)) + 
  geom_polygon(colour='black') + coord_equal() + theme_void() + scale_fill_brewer(name="",palette = "Spectral", labels=c("Grupo 1", "Grupo 2", "Grupo 3", "Grupo 4"))


library(FactoMineR)   # ACP
library(ggplot2)      # Gr?ficos
library(maptools)     # Mapas
library(rgdal)        # Mapas
library(dplyr)        # Mapas
library(rgeos)

load("/home/patricia/Dropbox/nupis/projetos/2015-larissa/sul_dem.rda")

# tirar pop
sul <- sul[,-13]
# nomes dos munic??pios como rótulos
sul <- as.data.frame(sul)
row.names(sul) <- sul$nome.mun

# retirar as variáveis correlacionadas com outras (mort5 e sobre40)
# e que, de alguma forma, calculam a mesma coisa
names(sul)
sul <- sul[,-c(8,10)]


# k-médias - k=2, k=3, k=4, k=5 ----------------------------------

set.seed(1)
grupos.1 <- kmeans(sul[,-c(1:4)], 1, nstart=25)
labk1 <- grupos.1$cluster

set.seed(1)
grupos.2 <- kmeans(sul[,-c(1:4)], 2, nstart=25)
labk2 <- grupos.2$cluster

set.seed(1)
grupos.3 <- kmeans(sul[,-c(1:4)], 3, nstart=25)
labk3 <- grupos.3$cluster

set.seed(1)
grupos.4 <- kmeans(sul[,-c(1:4)], 4, nstart=25)
labk4 <- grupos.4$cluster

set.seed(1)
grupos.5 <- kmeans(sul[,-c(1:4)], 5, nstart=25)
labk5 <- grupos.5$cluster

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


library(dplyr)
library(ggplot2)
library(rgeos) # se der erro ao instalar veja acima
library(rgdal) # se der erro ao instalar veja acima
library(maptools)


load("C:/Users/Larissa/Documents/mapa/atlas.rda")
mgm <-readOGR(dsn="C:/Users/Larissa/Documents/mapa",layer="31mu2500gsr")

mgm@data = rename(mgm@data, codmun7=GEOCODIG_M)

mgm1 = subset(mgm, MESORREGIã==3110)

mg = filter(atlas, nome.uf=='minas gerais')
mg = droplevels(mg)
mg$codmun7  = as.factor(mg$codmun7)
mgm$codmun7 = as.factor(mgm$codmun7)
mga = left_join(mgm@data, mg, by="codmun7")
mgf = fortify(mgm, region="codmun7")
mga$id = mga$codmun7
mmg = left_join(mgf,mga)

mgf1    = fortify(mgm1, region = "codmun7")
mgf1$id = as.factor(mgf1$id)
ssm     = filter(atlas, meso == 3110)
ssm$id  = ssm$codmun7

mssm    = inner_join(mgf1, ssm, by = "id")

cnomes = aggregate(cbind(long, lat) ~ nome.mun, data= mssm, FUN=mean)
ggplot(mssm, aes(long, lat)) +
  geom_polygon(aes(group=group), colour='lightgray', fill = NA) +
  coord_equal() + theme_void() +
  geom_text(data = cnomes, aes(long, lat, label = nome.mun), size=2, color='black', 
            check_overlap = TRUE)

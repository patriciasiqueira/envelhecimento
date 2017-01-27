# Dissertação - Larissa ---------------------------------------------------
# Dados do Sul/Sudoeste de MG sobre envelhecimento
# AA usando variáveis originais e escores dos componentes principais

# Manipulação dos dados
# Criação do RDATA sul_dem

# Diretórios de trabalho --------------------------------------------------
setwd("/home/patricia/Dropbox/nupis/projetos/2015-larissa")

# Conjunto de dados .RData ------------------------------------------------

load("atlas_dem.rda")

# Dados Sul Sudoeste de Minas -------------------------------------------------------------------

library(dplyr)
glimpse(atlas_dem)  # para visualizar as colunas e os tipos de dados
names(atlas_dem)
sul_dem <- subset(atlas_dem,meso == "3110")
dim(sul_dem)

# salvar o data frame com 95 variáveis
sul_dem_completo <- as.data.frame(sul_dem)
row.names(sul_dem_completo) <- sul_dem$nome.mun
# salvar RData
save(sul_dem_completo, file = "sul_dem_completo.rda")

# variáveis --------------------------------------------------------------
# 1) Esperança de vida ao nascer
# 2) Taxa de fecundidade total
# 3) Mortalidade infantil
# 4) Mortalidade até 5 anos de idade
# 5) Razão de dependência
# 6) Probabilidade de sobrevivência até 40 anos
# 7) Probabilidade de sobrevivência até 60 anos
# 8) Taxa de envelhecimento
# 9) População total
names(sul_dem)
sul <- sul_dem[,c(1,2,5,7,10:17,95)]

# tirar coluna do nome do município
# sul <- sul[,-2]
# nomes dos municípios como rótulos
sul_dem <- as.data.frame(sul_dem)
sul <- as.data.frame(sul)
row.names(sul) <- sul_dem$nome.mun

# salvar RData
save(sul, file = "sul_dem_nomes_antigos.rda")

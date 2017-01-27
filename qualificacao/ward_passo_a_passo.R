# método de Ward - passo a passo
# renda e idade de seis indivíduos
renda <- c(9.6,8.4,2.4,18.2,3.9,6.4)
idade <- c(28,31,42,38,25,41)
dados <- data.frame(renda,idade)
Xb <- colMeans(dados)
S <- cov(dados)
sqrt(S)
X <- as.matrix(dados)


# Passo 1 -----------------------------------------------------------------
# obter as distâncias entre os grupos e escolher
# a fusão que minimize o valor da distância
G1 <- t(as.matrix(X[1,]))
G2 <- t(as.matrix(X[2,]))
G3 <- t(as.matrix(X[3,]))
G4 <- t(as.matrix(X[4,]))
G5 <- t(as.matrix(X[5,]))
G6 <- t(as.matrix(X[6,]))
# SQi
# G1
i = 1
Xb <- X[c(i),]
SQ <- t(X[i,]-Xb)%*%(X[i,]-Xb)
SQ
# idem para todos G2, ..., G6
# SQT = soma dos SQi
SQT <- 0
# distâncias
# G1G2
r <- 1
s <- 2
nr <- nrow(G1)
ns <- nrow(G2)
Xbr <- colMeans(G1)
Xbs <- colMeans(G2)
d12 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d12
# sem correção
d12e <- t(Xbr-Xbs)%*%(Xbr-Xbs)
# G1G3
r <- 1
s <- 3
nr <- nrow(G1)
ns <- nrow(G3)
Xbr <- colMeans(G1)
Xbs <- colMeans(G3)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G1G4
r <- 1
s <- 4
nr <- nrow(G1)
ns <- nrow(G4)
Xbr <- colMeans(G1)
Xbs <- colMeans(G4)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G1G5
r <- 1
s <- 5
nr <- nrow(G1)
ns <- nrow(G5)
Xbr <- colMeans(G1)
Xbs <- colMeans(G5)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G1G6
r <- 1
s <- 6
nr <- nrow(G1)
ns <- nrow(G6)
Xbr <- colMeans(G1)
Xbs <- colMeans(G6)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G2G3
r <- 2
s <- 3
nr <- nrow(G2)
ns <- nrow(G3)
Xbr <- colMeans(G2)
Xbs <- colMeans(G3)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G2G4
r <- 2
s <- 4
nr <- nrow(G2)
ns <- nrow(G4)
Xbr <- colMeans(G2)
Xbs <- colMeans(G4)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G2G5
r <- 2
s <- 5
nr <- nrow(G2)
ns <- nrow(G5)
Xbr <- colMeans(G2)
Xbs <- colMeans(G5)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G2G6
r <- 2
s <- 6
nr <- nrow(G2)
ns <- nrow(G6)
Xbr <- colMeans(G2)
Xbs <- colMeans(G6)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G3G4
r <- 3
s <- 4
nr <- nrow(G3)
ns <- nrow(G4)
Xbr <- colMeans(G3)
Xbs <- colMeans(G4)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G3G5
r <- 3
s <- 5
nr <- nrow(G3)
ns <- nrow(G5)
Xbr <- colMeans(G3)
Xbs <- colMeans(G5)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G3G6
r <- 3
s <- 6
nr <- nrow(G3)
ns <- nrow(G6)
Xbr <- colMeans(G3)
Xbs <- colMeans(G6)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G4G5
r <- 4
s <- 5
nr <- nrow(G4)
ns <- nrow(G5)
Xbr <- colMeans(G4)
Xbs <- colMeans(G5)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G4G6
r <- 4
s <- 6
nr <- nrow(G4)
ns <- nrow(G6)
Xbr <- colMeans(G4)
Xbs <- colMeans(G6)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# G5G6
r <- 5
s <- 6
nr <- nrow(G5)
ns <- nrow(G6)
Xbr <- colMeans(G5)
Xbs <- colMeans(G6)
d <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d
# menor distância foi G1G2, ou seja, unir 1 e 2
# sem a correção fica
d12
# para dar igual ao R
sqrt(d12)
# com correção
sqrt(d12e)

# Passo 2 -----------------------------------------------------------------
# obter as distâncias entre os grupos e escolher
# a fusão que minimize o valor da distância
G1 <- as.matrix(X[1:2,])
G2 <- t(as.matrix(X[3,]))
G3 <- t(as.matrix(X[4,]))
G4 <- t(as.matrix(X[5,]))
G5 <- t(as.matrix(X[6,]))
# unindo 12-3
r <- 1
s <- 2
nr <- nrow(G1)
ns <- nrow(G2)
Xbr <- colMeans(G1)
Xbs <- colMeans(G2)
d123 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d123
# unindo 12-4
r <- 1
s <- 3
nr <- nrow(G1)
ns <- nrow(G3)
Xbr <- colMeans(G1)
Xbs <- colMeans(G3)
d124 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d124
# unindo 12-5
r <- 1
s <- 4
nr <- nrow(G1)
ns <- nrow(G4)
Xbr <- colMeans(G1)
Xbs <- colMeans(G4)
d125 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d125
# unindo 12-6
r <- 1
s <- 5
nr <- nrow(G1)
ns <- nrow(G5)
Xbr <- colMeans(G1)
Xbs <- colMeans(G5)
d126 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d126
# unindo 3-4
r <- 2
s <- 3
nr <- nrow(G2)
ns <- nrow(G3)
Xbr <- colMeans(G2)
Xbs <- colMeans(G3)
d34 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d34
# unindo 3-5
r <- 2
s <- 4
nr <- nrow(G2)
ns <- nrow(G4)
Xbr <- colMeans(G2)
Xbs <- colMeans(G4)
d35 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d35
# unindo 3-6
r <- 2
s <- 5
nr <- nrow(G2)
ns <- nrow(G5)
Xbr <- colMeans(G2)
Xbs <- colMeans(G5)
d36 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d36
d36e <- t(Xbr-Xbs)%*%(Xbr-Xbs)
# unindo 4-5
r <- 3
s <- 4
nr <- nrow(G3)
ns <- nrow(G4)
Xbr <- colMeans(G3)
Xbs <- colMeans(G4)
d45 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d45
# unindo 4-6
r <- 3
s <- 5
nr <- nrow(G3)
ns <- nrow(G5)
Xbr <- colMeans(G3)
Xbs <- colMeans(G5)
d46 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d46
# unindo 5-6
r <- 4
s <- 5
nr <- nrow(G4)
ns <- nrow(G5)
Xbr <- colMeans(G4)
Xbs <- colMeans(G5)
d34 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d34
# 3-6 foram unidas, para dar igual ao R
sqrt(d36e)
# com correção
sqrt(d36)

# Passo 3 -----------------------------------------------------------------
G1 <- as.matrix(X[1:2,])
G2 <- as.matrix(X[c(3,6),])
G3 <- t(as.matrix(X[4,]))
G4 <- t(as.matrix(X[5,]))
# unindo 12-36
r <- 1
s <- 2
nr <- nrow(G1)
ns <- nrow(G2)
Xbr <- colMeans(G1)
Xbs <- colMeans(G2)
d1236 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d1236
# unindo 12-4
r <- 1
s <- 3
nr <- nrow(G1)
ns <- nrow(G3)
Xbr <- colMeans(G1)
Xbs <- colMeans(G3)
d124 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d124
# unindo 12-5
r <- 1
s <- 4
nr <- nrow(G1)
ns <- nrow(G4)
Xbr <- colMeans(G1)
Xbs <- colMeans(G4)
d125 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d125
d125e <- t(Xbr-Xbs)%*%(Xbr-Xbs)
# unindo 36-4
r <- 2
s <- 3
nr <- nrow(G2)
ns <- nrow(G3)
Xbr <- colMeans(G2)
Xbs <- colMeans(G3)
d364 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d364
# unindo 36-5
r <- 2
s <- 4
nr <- nrow(G2)
ns <- nrow(G4)
Xbr <- colMeans(G2)
Xbs <- colMeans(G4)
d365 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d365
# unindo 4-5
r <- 3
s <- 4
nr <- nrow(G3)
ns <- nrow(G4)
Xbr <- colMeans(G3)
Xbs <- colMeans(G4)
d45 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d45
# para dar igual ao R
sqrt(d125e)


# Passo 4 -----------------------------------------------------------------
G1 <- as.matrix(X[c(1:2,5),])
G2 <- as.matrix(X[c(3,6),])
G3 <- t(as.matrix(X[4,]))
# unindo 125-36
r <- 1
s <- 2
nr <- nrow(G1)
ns <- nrow(G2)
Xbr <- colMeans(G1)
Xbs <- colMeans(G2)
d12536 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d12536
# unindo 125-4
r <- 1
s <- 3
nr <- nrow(G1)
ns <- nrow(G3)
Xbr <- colMeans(G1)
Xbs <- colMeans(G3)
d1254 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d1254
# unindo 36-4
r <- 2
s <- 3
nr <- nrow(G2)
ns <- nrow(G3)
Xbr <- colMeans(G2)
Xbs <- colMeans(G3)
d364 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d364
d364e <- t(Xbr-Xbs)%*%(Xbr-Xbs)
# unindo 36-4
r <- 2
s <- 3
nr <- nrow(G2)
ns <- nrow(G3)
Xbr <- colMeans(G2)
Xbs <- colMeans(G3)
d364 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d364
# unindo 36-5
r <- 2
s <- 4
nr <- nrow(G2)
ns <- nrow(G4)
Xbr <- colMeans(G2)
Xbs <- colMeans(G4)
d365 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d365
# unindo 4-5
r <- 3
s <- 4
nr <- nrow(G3)
ns <- nrow(G4)
Xbr <- colMeans(G3)
Xbs <- colMeans(G4)
d45 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d45
# igual do R
sqrt(d364e)


# Passo 5 -----------------------------------------------------------------
G1 <- as.matrix(X[c(1:2,5),])
G2 <- as.matrix(X[c(3,6,4),])
# unindo 125-364
r <- 1
s <- 2
nr <- nrow(G1)
ns <- nrow(G2)
Xbr <- colMeans(G1)
Xbs <- colMeans(G2)
d125364 <- (nr*ns/(nr+ns))*t(Xbr-Xbs)%*%(Xbr-Xbs)
d125364
d125364e <- t(Xbr-Xbs)%*%(Xbr-Xbs)
sqrt(d125364e)


# resultado similar ao obtido pela Mingoti
dados <- data.frame(renda,idade)
dis <- dist(dados)^2
ward <- hclust(dis, method = "ward.D2")
plot(ward <- hclust(dis, method = "ward.D2"))
names(ward)
ward$merge
ward$height
ward$dist.method

# Ward.d, euclidiana ao quadrado
# igual da Mingoti
dados <- data.frame(renda,idade)
dis <- dist(dados)^2
ward <- hclust(dis, method = "ward.D")
plot(ward <- hclust(dis, method = "ward.D"))
names(ward)
ward$merge
ward$height
ward$dist.method

# clássico do R
# ward.D2, distância euclidiana
dados <- data.frame(renda,idade)
dis <- dist(dados)
ward <- hclust(dis, method = "ward.D2")
plot(ward <- hclust(dis, method = "ward.D2"))
names(ward)
ward$merge
ward$height
ward$height^2
ward$dist.method

#
clusters = cutree(hclust(dist(dados)^2), k=3)

# function to find medoid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}

sapply(unique(clusters), clust.centroid, dados, clusters)



# Carrega pacote para construção filogenética
library(ape)

# Carrega função para fazer gráfico
source("https://raw.githubusercontent.com/csdambros/R-functions/master/poncho.R")


## Simulando alguns dados com sinal filogenética (espécies próximas "gostam" do mesmo ambiente)

# 1. Cria uma filogenia aleatória (só para usar de exemplo)
phy<-rcoal(50)
plot(phy)

# 2.Cria um gradiente ambiental (só para exemplo)
env<-seq(10,40,length.out = 30)

# 3.Cria o ótimo ambiental de cada espécie ao longo do gradiente e coloca na ordem dos tips da phylo
envSp<-seq(1,30,length.out = 50)[order(phy$tip.label)]

# 4.Adiciona sinal filogenético na associação das espécies com o ambiente
phydist<-cophenetic(phy) # Calcula distância filog.
phydist.inv<-1-((phydist)/rowSums((phydist))) # Inverso da distância padronizada (quanto mais distante, menor o sinal)

phydist.col<-chol(phydist.inv)# Cholesky decomposition da matrix de distância inversa

# Recalcula a preferência ambiental de cada espécie com sinal filogenético
envSp2<-matrix(envSp,1,50)%*%phydist.col


## Relação da preferência ambiental e filogenia

# Sem sinal filogenético
plot(phy,tip.color = heat.colors(100)[cut(envSp,100)])

# Com forte sinal filogenético
plot(phy,tip.color = heat.colors(100)[cut(envSp2,100)])


# 5. Cria matrix locais x espécies 
envSp2mat<-t(matrix(envSp2,50,30))
envNorm<-dnorm(envSp2mat,env,3)

# Probabilidade de ocorrência de cada espécie ao longo do gradiente ambiental
matplot(env,envNorm,type="l")

# Probabilidade padronizada (entre zero e 1)
envNormRange<-envNorm/apply(envNorm,1,max)

# Cria matriz de presença e ausência
PA<-matrix(rbinom(envNorm,1,envNormRange),length(env))
colnames(PA)<-phy$tip.label


## Retomando:
# Existem 3 objetos:
# 1. Matriz de presença e ausência
# 2. Variável ambiental (gradiente)
# 3. Filogenia

# 6. Faz ordenação e plota a ocorrência das espécies ao longo do gradiente e da filogenia
poncho(PA,ylab.bottom = "PCA1") # Somente ordenação da ocorrência
poncho(PA,phy = phy)
poncho(PA,env = env,phy = phy)

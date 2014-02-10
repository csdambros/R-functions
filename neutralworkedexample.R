##############################################
#Worked examples of a metacommunity neutral model simulation using a network approach (semi analytical approach).
#See Economo and Keitt (2008)

#Developed by Dambros CS in 05 March 2013
#Last modified by Dambros CS in 06 March 2013

##############################################

#Install the anetwork function (first download 'anetwork.R' to your folder):

	source ('anetwork.R')

# Imagine three areas that are connected to each other trough migration
# What would be the species diversity in each area and the turnover of species among ares just under neutral processes? and incorporating environmental differences?
# These questions can be answered using the anetwork function!


#######################################################
# First example

#Creating a migration matrix where each row is an area receiving immigrants and columns are areas donating emigrants
#In the following example the matrix is simetrical (there is no source-sink dynamic). Lets start a model with arbitrary values just to exemplify

	M <- matrix(0.02,3,3) #set the migration rates among all ares to be 0.002
	diag(M) <- 1-rowSums(M)+diag(M) #set the non-migration probability (1 - migration probability)

	M #see how the migration matrix looks like

	drawnet(M)

#The probability of receiving immigrants + the probability of not receiving immigrants must be 1 (100%)
#Then the rowSums of this table MUST be all equals 1

	rowSums(M)

#######################################################

#set the number of individuals in each area to be the same
	
	N=c(100,100,100);N #Local community sizes
	N=100 # The same

#set the speciation rate

	v=0.001;v

######################################################

	model.1<-anetwork(N,M,v) # Run the model

# This model started with just one species in all areas. After some time, species were added by points speciation
# This means that in the beginning the diversity in each area is zero (probability of randomly picking 2 individuals of different species or Simpson's diversity index)

# The mean local diversity (alpha diversity) is shown in green and starts at zero
# Initially all areas have the same species, then the initial species turnover (beta diversity) is also zero

# The overall diversity (gamma diversity) can be understood as the sum of alpha and beta diversities (additive partitioning, see Jost 2006 and Tuomisto 2011)
# Gamma = alpha + beta

# Gamma (t0) = 0+0 = 0

# Gamma diversity is shown in red on the figure

#######################################################

# As species are added by speciation events, the local diversity (green) increases along generations (x axis)
# As the species that appear on a local community have a limited dispersion, this creates differences among local communities (beta diversity or species turnover)

# If species turnover is not zero, then the overall diversity (red) is always higher than local diversities (green), as shown on the figure

# After some time, the lost of species by drift (random process) and the gain by speciation reaches a steady-state.
# In other words, the local species diversity and turnover do not change anymore

# 1. The anetwork function provides the probability of randomly picking 2 individuals of different species from the same local community on each local community
# 2.The anetwork function provides the probability of randomly picking 2 individuals of different species from 2 different local communities for each pair of local communities

# With these information is possible to calculate many community parameters given migration, speciation, and the number of individuals in each local community, such as: 
# local species diversity (as simpson concentration index);
# species turnover (as morisita-horn overlap index);
# the expected number of species when all species have the same abundance (Hill numbers);
# which areas share more species
# which areas harbor the highest species diversity
# ...

# Visualize the expected local diversity in each local community:

	1-diag(model.1$finalF)

# As all communities have the same number of individuals in this initial model, receive the same amount of immigrants and have the same speciation rates
# The diversities are all the same for all areas. These values will be different if the areas have different immigration rates (M)


# Calculate the Morisita-Horn index of overlap among community pairs:

	extractMH(model.1$finalF)

# The overlap of a community with itself is 1 (100%) and is very high among communities in this model
# There is not 2 communities more similar to each other because all areas have exactly all the same parameters


###########################################################
# Second example: Some areas receive more immigrants

# Now let's imagine that we have areas that are more isolated than others, how can we include this in our model?
# We just need to change our migration matrix (M). For now let's keep the matrix symmetrical (non directional migration)
# We will simulate 3 areas again, but now one of those areas (area 3) is almost isolated from the others

	M <- matrix(0.02,10,10) #set the migration rates among all areas to be 0.002
	M [3,1:2] <-0.00001
	M [1:2,3] <-0.00001
	diag(M) <- 1-rowSums(M)+diag(M) #set the non-migration probability (1 - migration probability)

	M #see how the migration matrix looks like - Most of individuals in the area 3 are not immigrants
	drawnet(M)

	rowSums(M) # immigration rates + non-immigration rates = 1


######################################################

	model.2<-anetwork(N,M,v) # Run the model


# On the graph is possible to observe that the final overall diversity (red) is much higher than the local diversities (green)
# This indicates that now we have more species turnover (difference between red and green)

	1-diag(model.2$finalF)  # The local diversity is lower on the isolated area

	1-extractMH(model.2$finalF) # The areas 1 and 2 have much more in common to each other than with area 3


######################################################

#summarizing the species composition change using a PCA...

PCoA<-cmdscale(1-extractMH(model.2$finalF)) # subtract MH from 1 because a matrix of distance (dissimilarity) must be provided

plot(PCoA[,1],PCoA[,2],type='n',xlim=c(-0.4,0.7),ylim=c(-0.4,0.7))
text(PCoA[,1],PCoA[,2],1:3)

# As expected, areas 1 and 2 have a much more similar species composition

#######################################################
# Third example: Integrating simulation and reality

# Up to now, we saw how the model works when we provide the migration rates (M),the speciation rates (V) and the number of individuals (N)
# However, we usually don't know such parameters!!

# How can we then simulate our community and compare the results with empirical data?

# One way to do this is to test the best combination of the parameters M, V and N that generate our observed species diversity patterns

# For this example, we will use an optimization function that will look for those values for us

# We will start with the previous result from the anetwork function as if we had collected those values for local diversity and turnover
# We will then pretend that we don't know M, V and N and we will see if we can discover those values just looking at our diversities

	m=0.002
	v=0.001
	N=100

	M.structure<- matrix(c(0,1),3,3)

	M <- matrix(m,3,3)*M.structure #set the migration rates among all ares to be 0.002
	diag(M) <- 1-rowSums(M)+diag(M) #set the non-migration probability (1 - migration probability)

	M #see how the migration matrix looks like

#optimization process

	model.3<-anetwork(N,M,v)

	comm<-(1-model.3$finalF)

	diff.anetwork(comm,N,M,v,nruns=10,stepshow=1)

	ver<-optim(c(0.01,0.01),function(x){

			M <- matrix(x[1],3,3)*M.structure #set the migration rates among all ares to be 0.002
			diag(M) <- 1-rowSums(M)+diag(M) #set the non-migration probability (1 - migration probability)

			diff.anetwork(comm,100,M,x[2],nruns=1000,stepshow=1000)

			})
	
	estimated.m<-ver$par[1]
	estimated.v<-ver$par[2]

	M2<-matrix(estimated.m,nrow(comm),ncol(comm))*M.structure
	diag(M2)<- 1-rowSums(M2)+diag(M2)

	model.3.1<-anetwork(100,M2,estimated.v)

	PCoA3simu<-cmdscale(1-extractMH(model.3.1$finalF))
	PCoA3real<-cmdscale(1-extractMH(1-comm))

	plot(PCoA3real[,1],PCoA3real[,2])
	points(PCoA3simu[,1],PCoA3real[,2],pch=21,bg=1)

	plot(PCoA3real[,1],PCoA3simu[,1])
	text(PCoA3real[,1],PCoA3simu[,1],1:3)

	sum(sqrt((model.3.1$finalF-(1-comm))^2)) #total difference on observed and predicted diversity

#########################################################################

# Using real data and predicting species distribution

#coming soon


# Interesting literature:

# Chave et al. 2004
#










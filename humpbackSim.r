#Humpback whale population impact simulation
#Used in: Williams et al. "Gauging allowable harm limits to cumulative, sub-lethal 
# effects of human activities on wildlife populations"

#Author: Len Thomas

#All relevant driver code is at the bottom of the file, after the functions are defined.
#The function definitions are generic, but code has been modified to fit the humpback
# case study.

#Implementation note: makes extensive use of global variables to save passing 
# large structures into functions (since R does calling by value this would be
# very inefficient)

#----------------------------------------------------------------------------
# Function definitions
#----------------------------------------------------------------------------

do.sim<-function(){
#Purpose: runs the population dynamics simulation
# For each year, runs with and without prey reduction
#Inputs:
# Assumes global variables:
#  g.par -list of global parameters affecting the simulation
#  g.pop.M - array of numbers of males by age, year and simulation
#  g.pop.F - array of numbers of females by age, year and simultion
#  g.env - matrix of environment indices by year and simulation
#Outputs:
# Fills g.pop.M and g.pop.F with outcome of g.par$n.sims simulations

#Implementation note: In this simulation, the reference.year passed in to
# do.survival.and.ageing is always year 1.  Because of this, each "year"
# is effectively just a replicate -- you're always simulating forward from
# year 1 to year 2, n.years times
  
  for(sim in 1:g.par$n.sims){
    
    #simulate population sizes in year 1
    init.pop(sim,T)
    
    #simulate population sizes for the other years
    for(year in 2:g.par$n.years){
      #simulate forward from year 1 for no prey reduction pop
      do.survival.and.ageing(year,sim,reference.year=1)
      do.fecundity(year,sim)
      #simulate forward from year 1 for prey reduction pop
      do.survival.and.ageing(year,sim,reference.year=1,do.prey.red=T)
      do.fecundity(year,sim,do.prey.red=T)
    }
  }
  
}

#----------------------------------------------------------------------------

do.stableage.sim<-function(){
#Purpose: a version of the simulation to establish stable age structure
# for starting all the future simulations
#Inputs:
# Assumes global variables g.par, g.pop.M, g.pop.F and g.env
#Outputs:
# Fills g.pop.M and g.pop.F with outcome of g.par$n.sims simulations

#Implementation note: unlike the previous simulation, this time you simulate
# forward in time: year 1, year 2, year 3, etc -- this gives the population 
# time to reach its stable age structure, which is what you want to obtain
  
  for(sim in 1:g.par$n.sims){
    #simulate population sizes in year 1
    init.pop(sim,F)
    
    #simulate population sizes for the other years
    for(year in 2:g.par$n.years){
      #simulate forward one year
      do.survival.and.ageing(year,sim)
      do.fecundity(year,sim)
    }
  }
  
}

#----------------------------------------------------------------------------

init.pop<-function(sim,init.with.prey.red.pop){
#Purpose: initialize the population values in year 1 for simulation sim
# Requires g.par$stable.age.prop.M and .F
#Inputs: 
#  sim - simulation number
#  init.with.prey.red.pop - if T then sets g.pop.M.with.prey.red and g.pop.F.with.prey.red to be
#   the same as the no prey reduction population
#  Assumes same global vars as do.sim
#Outputs:
#  g.pop.M[,1,sim] and g.pop.F[,1,sim] have values in them
  
  #Initial number of females samples from a binomial
  n.females<-rbinom(1,g.par$init.pop,g.par$stable.age.proportion.females)
  n.males<-g.par$init.pop-n.females
  g.pop.M[1:g.par$n.ages,1,sim]<<-rmultinom(1,g.par$init.pop/2,g.par$stable.age.prop.M)
  g.pop.F[1:g.par$n.ages,1,sim]<<-rmultinom(1,g.par$init.pop/2,g.par$stable.age.prop.F)
  if(init.with.prey.red.pop){
    g.pop.M.with.prey.red[1:g.par$n.ages,1,sim]<<-g.pop.M[1:g.par$n.ages,1,sim]
    g.pop.F.with.prey.red[1:g.par$n.ages,1,sim]<<-g.pop.F[1:g.par$n.ages,1,sim]
  }
  
}

#----------------------------------------------------------------------------

do.survival.and.ageing<-function(year,sim,reference.year=year-1,do.prey.red=F){
#Purpose: simulates the process of getting one year older and surviving
# for all age classes of males and females in the specified year and simulation
# Note: prey.red can be between 0 and 1 -- the effect of this is to reduce
# the quality of the environment-- effectively
#  environment is env*(1-prey reduction)
#  Assumes same global vars as do.sim
#Inputs:
# year - year number (or replicate number if reference year is always 1)
# sim - sim number
# reference year - year to project forward from
# do.prey.red - if T then implement prey reduction
#Outputs:
# sets the global variables g.pop.M and g.pop.F for that year (or, if do.prey.red is T then
# sets g.pop.M.with.prey.red and g.pop.F.with.prey.red)

#Implementation note: prey reduction affects 1st year survival (i.e. calf to year 1)
# This info comes from Robbins (2007) PhD Thesis, Figure 4.6

#Implementation note: it's basically the same code twice, for males and
# females, which is not good programming practice -- done this way to
# increase code speed.  (Also could be subdivided into multiple functions...)
  
  #Work out the relevant survival vectors, given the environment, using linear
  # interpolation/extrapolation between the restricted and unrestricted survivals
  prey.red<-do.prey.red*g.par$prey.red
  surv.F<-g.par$F.surv.intercept+g.env[year,sim]*(1-prey.red)*g.par$F.surv.slope
  surv.M<-g.par$M.surv.intercept+g.env[year,sim]*(1-prey.red)*g.par$M.surv.slope
  #if you extrapolate, it's possible to get inadmissible values, so set to within range
  surv.F[surv.F>1]<-1;surv.F[surv.F<0]<-0
  surv.M[surv.M>1]<-1;surv.M[surv.M<0]<-0
  
  #Survival is a binomial process
  oldest.F<-g.pop.F[g.par$n.ages,reference.year,sim]
  oldest.M<-g.pop.M[g.par$n.ages,reference.year,sim]
  if(do.prey.red){
    #do those from ages 1-(n.ages-1)
    g.pop.F.with.prey.red[2:g.par$n.ages,year,sim]<<-rbinom(g.par$n.ages-1,g.pop.F[1:(g.par$n.ages-1),reference.year,sim],surv.F[-g.par$n.ages])
    g.pop.M.with.prey.red[2:g.par$n.ages,year,sim]<<-rbinom(g.par$n.ages-1,g.pop.M[1:(g.par$n.ages-1),reference.year,sim],surv.M[-g.par$n.ages])
    #do the last age group (those of n.ages)
    g.pop.F.with.prey.red[g.par$n.ages,year,sim]<<-g.pop.F.with.prey.red[g.par$n.ages,year,sim]+rbinom(1,oldest.F,surv.F[g.par$n.ages])
    g.pop.M.with.prey.red[g.par$n.ages,year,sim]<<-g.pop.M.with.prey.red[g.par$n.ages,year,sim]+rbinom(1,oldest.M,surv.F[g.par$n.ages])
  } else {
    #do those from ages 1-(n.ages-1)
    g.pop.F[2:g.par$n.ages,year,sim]<<-rbinom(g.par$n.ages-1,g.pop.F[1:(g.par$n.ages-1),reference.year,sim],surv.F[-g.par$n.ages])
    g.pop.M[2:g.par$n.ages,year,sim]<<-rbinom(g.par$n.ages-1,g.pop.M[1:(g.par$n.ages-1),reference.year,sim],surv.M[-g.par$n.ages])
    #do the last age group (those of n.ages)
    g.pop.F[g.par$n.ages,year,sim]<<-g.pop.F[g.par$n.ages,year,sim]+rbinom(1,oldest.F,surv.F[g.par$n.ages])
    g.pop.M[g.par$n.ages,year,sim]<<-g.pop.M[g.par$n.ages,year,sim]+rbinom(1,oldest.M,surv.F[g.par$n.ages])
  }
  
}
#----------------------------------------------------------------------------

do.fecundity<-function(year,sim,do.prey.red=F){
#Purpose: simulates the fecundity process
#  Assumes same global vars as do.sim
#Inputs: 
# year - year number (or replicate number if reference year is always 1)
# sim - sim number
# reference year - year to project forward from
# do.prey.red - if T then implement prey reduction
#Outputs:
# sets the global variables g.pop.M and g.pop.F for that year (or, if do.prey.red is T then
# sets g.pop.M.with.prey.red and g.pop.F.with.prey.red)

#Implementation note: no effect of prey reduction on fecundity in this case
# This info comes from Robbins (2007) PhD Thesis, Figure 4.6
  
  #Get environment-specific vector of calving probabilities
  prey.red<-do.prey.red*g.par$prey.red
  #Note - no effect of pred reduction on fecundity rate
  fecundity.rate<-p.calf(1:g.par$n.ages,g.env[year,sim])
  #Number of calves is a binomial random variable for each age class of females
  n.calves<-sum(rbinom(g.par$n.ages,g.pop.F[,year,sim],fecundity.rate))
  if(do.prey.red){
    #Assume expected sex ratio age 0.5 is 50:50, so number of females is Bernoulli(0.5)...
    g.pop.F.with.prey.red[1,year,sim]<<-rbinom(1,n.calves,0.5)
    #... and those not female are male
    g.pop.M.with.prey.red[1,year,sim]<<-n.calves-g.pop.F.with.prey.red[1,year,sim]
  } else {
    #Assume expected sex ratio age 0.5 is 50:50, so number of females is Bernoulli(0.5)
    g.pop.F[1,year,sim]<<-rbinom(1,n.calves,0.5)
    #... and those not female are male
    g.pop.M[1,year,sim]<<-n.calves-g.pop.F[1,year,sim]
  }  
  
}

#----------------------------------------------------------------------------

p.calf<-function(age,env){
#Purpose: Implements the fecundity schedule - in the case of humpbacks it is not
# dependent on the environment
# Assumes global var g.par with list elements $age.first.breed and $fecundity
#Inputs: age - age; env - environment value
#Outputs: probability of producing a calf
  
  #up to age 7 (i.e., 6 as first year is age 1), no breeding, thereafter a constant rate
  return(ifelse(age>=g.par$age.first.breed,g.par$fecundity,0))
}

#----------------------------------------------------------------------------

get.PBR<-function(N){
#Purpose: Returns Potential Biological Removal given N (can be a vector)
#Inputs:
# N population size (can be vector or matrix)
# Assumes g.par$PBR.C and $PBR.R.max have been defined
#outputs: 
# PBR value given inputs
  
  #Work out lognormal LCL given N CV and alpha
  LCL<-N/g.par$PBR.C
  
  return(LCL*0.5*g.par$PBR.R.max*g.par$PBR.f)
  
}

#----------------------------------------------------------------------------
# Driver code
#----------------------------------------------------------------------------

#Set simulation parameters
#-------------------------

g.par<-list(
  #number of age classes -- assume same for males and females, for simplicity
  n.ages=7,
  #number of time periods to simulate
  n.years=100,
  #number of simulations to do
  n.sims=50,
  #Good environment - value where unrestricted growth occurs. 
  #For humpbacks this is in terms of sand lance numbers in 1988
  env.good=3158,
  #Poor environment - value where restricted growth occurs
  #Sand lance numbers in 1990
  env.poor=19,
  #prey reduction metric - 0 is no prey reduction, 1 is total prey reduction
  #note that we are therefore extrapolating the environment
  prey.red=0.0,
  #Age at first breeding 
  age.first.breed=7, #Note this is 6 as first year (age 0) is denoted 1 here
  #fecundity
  fecundity=1/2.36,
  #Initial population size
  #must be an even number - half will be male and half female
  init.pop=1500,
  #Alpha level to use in calculating LCL for PBR calc
  PBR.alpha=0.4,
  #Rmax for use in calculting PBR
  PBR.R.max=0.04,
  #f (risk or recovery factor)
  PBR.f=0.5,
  #CV for use in calculating PBR
  PBR.CV.N=0.20,
  #Alpha level for prey reduction effect CI
  prey.red.alpha=0.05
)
#Put the PBR stuff together to calculate a C for use in working out the Lognormal
#LCL on N
PBR.var.ln.N<-log(1+g.par$PBR.CV.N^2)
g.par$PBR.C<-exp(qnorm(1-g.par$PBR.alpha/2) * sqrt(PBR.var.ln.N))

#Enter in the survival parameters
g.par$F.surv.restricted<-c(0.63,rep(0.96,6))
g.par$F.surv.unrestricted<-c(1.0,rep(0.96,6))
g.par$M.surv.restricted<-g.par$F.surv.restricted
g.par$M.surv.unrestricted<-g.par$F.surv.unrestricted

#Work out the intercept and slope for male survival and female survival as a function of environment, 
# assuming a linear relationship between survival and environment -- do
# it now to save working it out multiple times later
g.par$F.surv.slope<-(g.par$F.surv.unrestricted-g.par$F.surv.restricted)/(g.par$env.good-g.par$env.poor)
g.par$F.surv.intercept<-g.par$F.surv.restricted-g.par$F.surv.slope*g.par$env.poor
g.par$M.surv.slope<-(g.par$M.surv.unrestricted-g.par$M.surv.restricted)/(g.par$env.good-g.par$env.poor)
g.par$M.surv.intercept<-g.par$M.surv.restricted-g.par$M.surv.slope*g.par$env.poor


#Run initial simulations to get the stable age structure
#-------------------------------------------------------

#Temporarily store vars
n.years.bak<-g.par$n.years
n.sims.bak<-g.par$n.sims
g.par$n.years<-200
g.par$n.sims<-1

#Define a temporary stable age structure - same number in all categories
g.par$stable.age.prop.F<-rep(1/g.par$n.ages,g.par$n.ages)
g.par$stable.age.prop.M<-rep(1/g.par$n.ages,g.par$n.ages)
#Define a starting sex ratio - proportion of the population that is female
g.par$stable.age.proportion.females<-0.5
#Define structure to hold males and females for sim
g.pop.M<-g.pop.F<-array(0,c(g.par$n.ages,g.par$n.years,g.par$n.sims))
#Set environment to mean environment
g.env<-matrix(mean(c(g.par$env.good,g.par$env.poor)),g.par$n.years,g.par$n.sims)
#Run the simultion
do.stableage.sim()
#Store the stable age structure
g.par$stable.age.prop.F<-g.pop.F[,g.par$n.years,g.par$n.sims]/sum(g.pop.F[,g.par$n.years,g.par$n.sims])
g.par$stable.age.prop.M<-g.pop.M[,g.par$n.years,g.par$n.sims]/sum(g.pop.M[,g.par$n.years,g.par$n.sims])
#Restore settings
g.par$n.years<-n.years.bak; rm(n.years.bak)
g.par$n.sims<-n.sims.bak; rm(n.sims.bak)

#Run main simulation
#--------------------

#Define global structures that will be used through the simulation
#Arrays to hold males and females
#With and without prey reduction
g.pop.M<-g.pop.F<-array(0,c(g.par$n.ages,g.par$n.years,g.par$n.sims))
g.pop.M.with.prey.red<-g.pop.F.with.prey.red<-array(0,c(g.par$n.ages,g.par$n.years,g.par$n.sims))

#Environment - i.e., in this case Sand Lance
#sample values from from a uniform around the max and min
tmp<-g.par$n.years*g.par$n.sims
g.env<-matrix(runif(tmp,min=g.par$env.poor,max=g.par$env.good),g.par$n.years,g.par$n.sims)

prey.red.levels<-seq(0,1,by=0.1)
n.prey.red.levels<-length(prey.red.levels)
#store the prey reduction effect metric values and conf limits
prey.red.effect<-prey.red.effect.LCL<-prey.red.effect.UCL<-numeric(n.prey.red.levels)

for(i in 1:n.prey.red.levels){

  cat("Prey reduction",i,"\n")
  
  g.par$prey.red<-prey.red.levels[i]
  #Run the simulation
  do.sim()
  
  #Work out the number lost each year due to prey reduction
  #(discount the first year, as pops the same then)
  diff.pop.M<-g.pop.M[,-1,]-g.pop.M.with.prey.red[,-1,]
  diff.pop.F<-g.pop.F[,-1,]-g.pop.F.with.prey.red[,-1,]
  diff.pop<-diff.pop.M+diff.pop.F
  
  diff.pop.sumages<-apply(diff.pop,c(2,3),sum)
  #Here you have a matrix (years x sims)
  #Need to work out the PBR for each year (using g.pop.M)
  # and then can work out the prey reduction removal as a function of the PBR
  
  #work out total pop size, from year 2 onwards
  tot.pop<-apply(g.pop.M[,-1,],c(2,3),sum)+apply(g.pop.F[,-1,],c(2,3),sum)
  #... and so get the PBR
  PBR<-get.PBR(tot.pop)
  
  #Work out the difference in pop numbers due to prey reduction divided by the
  # potential biological removal
  diff.prop.pop<-diff.pop.sumages/PBR
  #work out mean
  mean.diff.prop.pop<-apply(diff.prop.pop,1,mean)
  
  #Plot, as a diagnostic tool - commented out in main simulation
  #plot(2:g.par$n.years,mean.diff.prop.pop,ylim=c(min(LCL.diff.prop.pop),max(UCL.diff.prop.pop)),type="l",main=paste("prey.red=",prey.red.levels[i],sep=""))
  #Can also work out CI and plot
  #LCL.diff.prop.pop<-apply(diff.prop.pop,1,quantile,0.025)
  #UCL.diff.prop.pop<-apply(diff.prop.pop,1,quantile,0.975)
  #lines(2:g.par$n.years,LCL.diff.prop.pop,lty=2)
  #lines(2:g.par$n.years,UCL.diff.prop.pop,lty=2)
  
  prey.red.effect[i]<-mean(mean.diff.prop.pop)
  prey.red.effect.LCL[i]<-quantile(mean.diff.prop.pop,g.par$prey.red.alpha/2)
  prey.red.effect.UCL[i]<-quantile(mean.diff.prop.pop,1-g.par$prey.red.alpha/2)

}

#Display results
#---------------

#Smooth the results to create a nicer plot, and enable estimation of intermediate points
library(mgcv)
prey.red.effect.smooth<-gam(prey.red.effect~s(prey.red.levels))
prey.red.effect.LCL.smooth<-gam(prey.red.effect.LCL~s(prey.red.levels))
prey.red.effect.UCL.smooth<-gam(prey.red.effect.UCL~s(prey.red.levels))

#plot
old.par<-par(no.readonly=T)
pdf(file="humpback.pdf",width=7,height=5)
plot.range<-range(c(predict(prey.red.effect.smooth),predict(prey.red.effect.LCL.smooth),
                    predict(prey.red.effect.UCL.smooth)))

plot(prey.red.levels,predict(prey.red.effect.smooth),type="n",ylim=plot.range,ylab="Population level effect",xlab="Prey reduction")
polygon(c(prey.red.levels,rev(prey.red.levels)),
        c(predict(prey.red.effect.LCL.smooth),rev(predict(prey.red.effect.UCL.smooth))),
        col=grey(0.8),border=grey(0.8))
lines(prey.red.levels,predict(prey.red.effect.smooth))
abline(h=1,lty=2)

#put graphics params back
par(old.par)
dev.off()


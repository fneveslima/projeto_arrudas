#### running simulation Gr4H 
 
setwd('C:/Data/fernando.lima/airGR/simulations_gr4h')
library("airGR")

### load of catchment data
inputs <- read.csv(file = "input_gr4h_simulationA.csv", header = TRUE, sep = ";")

# Initial Date 
dt1 <- as.POSIXct("2006/10/10", format = "%Y/%m/%d", tz = "UTC")
# Final Date 
n_days <- 3005 ### total of days of data INMET (oct 2006 - dec 2014)
n_secs <- n_days*24*3600
dt2 <- as.POSIXct(dt1) + n_secs

# Vector of date
dts <- seq(from = dt1, to = dt2 -1, by = "3600 secs")

### Input data
dts <- as.data.frame(dts)
dinput <- c(dts,inputs)
dinput <- as.data.frame(dinput)


## Fichier parameters
param <- read.table(file = "162parameters.txt", header = TRUE)
param <- as.matrix(param)

## run period selection
indrun_ini <-  as.POSIXct("11/10/2008 00:00", format = "%d/%m/%Y %H:%M", tz = "UTC")
indrun_final <- as.POSIXct("31/12/2014 23:00", format = "%d/%m/%Y %H:%M", tz = "UTC")
indrun <- seq(from = indrun_ini, to = indrun_final, by = "3600 secs")

### evporation data running
### difference between simulation time and running time
dfrt <- as.numeric(length(dts$dts) - length(indrun))
evp <- data.frame(dinput$ev..mm.h.)
evp2 <- data.frame(evp[-seq(1:dfrt),])


#### loading initial states 
load("iniprod162param.Rdata")
load("inirout162param.Rdata")

#### Simulating Events
### event date 	(day/month/year)
### A		 15/11/2012
### B		 03/12/2012	
### C		 12/12/2012	
### D		 07/01/2013	
### F    02/10/2013
### E2   09/04/2013

#### Line beggining of event
### search
lbegin <- which(as.character(indrun)=="2012-11-14 00:00:00")
lend <- which(as.character(indrun)=="2012-11-16 23:00:00")

## Get the level of reservoir production at start of the event
umid <- list()
for (k in 1:length(param[,1])){
  umid[[k]] <- iniprod[lbegin,k]
}


## Get the level of reservoir routing at start of the event
routt <- list()
for (k in 1:length(param[,1])){
  routt[[k]] <- inirout[lbegin,k]
}

### data input of event 
## Reading pluvio with thiessen for each subcatchment for each event
plu24 <- read.table("thiessen_plu_sub24_eventA.txt", header = TRUE)
## deleting the last line (not necessary)
plu24 <- plu24[-73,]

#### Extract evaporation from data input
evp_event <- evp2[lbegin:lend,]

# Initial Date
dtiA <- as.POSIXct("2012/11/14", format = "%Y/%m/%d", tz = "UTC")
# Final Date 
n_daysA <- 3 ### total of days of each event (3 days)
n_secsA <- n_daysA*24*3600
dtfA <- as.POSIXct(dtiA) + n_secsA

# Vector of date event 
dtsA <- seq(from = dtiA, to = dtfA -1, by = "3600 secs")


## Create a table with inputs of event
d24 <- data.frame(dtsA,plu24,evp_event)

### Period of running model - date event
Ind_Run_24 <- seq(which(format(d24$dtsA,format="%d/%m/%Y %H:%M")=="14/11/2012 00:00"),
                  which(format(d24$dtsA,format="%d/%m/%Y %H:%M")=="16/11/2012 23:00"))


## preparation of the InputsModel object
InputsModel24 <- CreateInputsModel(FUN_MOD = RunModel_GR4H,DatesR = d24$dtsA,
                                   Precip = d24$plu24,PotEvap = d24$evp_event)



### Table of contents
umide <- unlist(umid)
route <- unlist(routt)
cotr <- data.frame(param,umide,route)


## preparation of the RunOptions object and simulation
out_sim <- apply(cotr, 1, function(i){
  inrunoptions <- CreateRunOptions(FUN_MOD = RunModel_GR4H,
                                   InputsModel = InputsModel24,IndPeriod_Run = Ind_Run_24,
                                   IniResLevels = i[c("umide","route")])
  OutputsModel <- RunModel_GR4H(InputsModel = InputsModel24,RunOptions = inrunoptions, 
                                Param = c(i["X1"],i["X2"], i["X3"], i["X4"]))
  return(OutputsModel$Qsim)
})


# ### Rating Curve estimation using Manning Strickler equation
# #### Q = ks*(y^beta)*(I*0.5)
# ### Slope of cross section 24 is 0.6% or 0.006 m/m and beta in this case is equal 1.99 based in a HEC-RAS simulation
# slope24 <- 0.006
# ### ks = 1/manning manning varies 
# m24 <- seq(0.011,0.060,0.001)
# beta <- 2.02
# k24 <- 1/m24 


#### Rating Curve estimation using Manning Strickler equation
## reading parameters of station 24
p24 <- read.table(file = "param_station24.txt", header = TRUE)
m24 <- p24$manning
alpha <- p24$alpha
beta <- p24$beta


### Q = alpha*(y^beta)


## Reading water levels in event of subcatchment 24
w24 <- read.table(file = "depth_sub24hour_eventA.txt", header = TRUE)
wl24 <-  as.vector(w24$x)

# ### Creating a vector of flows in m3/s
# Q24 <- NULL
# 
# for (i in 1:length(k24)){
#   Q24[[i]] <- k24[i]*((wl24)^(beta))*(slope24^(0.5))
# }



#### Creating a vector of flows in m3/s

Q24 <- NULL

for (i in 1:length(m24)){
     Q24[[i]] <- alpha[i]*((wl24)^(beta[i]))
     }

#### Area of Subwatershed 24 in m²
area24 <- 48*(10^6)

### Converting flow in m3/s in mm/h dividing per area
facm <- 3600*1000/(area24) ### 3600 s in 1 hour and 1000 mm in 1 m

q24 <- NULL

for (k in 1:length(m24)){
  q24[[k]] <- Q24[[k]]*facm
  
}

### combination of possible paths
zz <- expand.grid(1:length(param[,1]),1:length(m24)) ### number of simulations x number rating curves

## number of simulations (number of rating curves x number os set parameters)
ns <- length(m24)*length(param[,1])

#### Calculating Nash criteria
### Function Nash
nse <- function(obs,sim){
  erz <- ((sim - obs)^2)
  nmra <- sum(erz)
  erobs <- (obs - (mean(obs)))^2
  dnmr <- sum(erobs)
  e1 <- 1 - (nmra/dnmr)
  return (e1)
}


nse24 <- NULL
for (i in 1:ns){
  nse24[i] <- nse(obs = q24[[zz[i,2]]], sim = out_sim[,zz[i,1]])
  
}

write.table (x = nse24, file = "nse_eventA_sub24_simulation162p.txt", quote = FALSE)


#### Calculating error
er24 <- NULL

for (i in 1:ns){
  er24[[i]] <- q24[[zz[i,2]]] - out_sim[,zz[i,1]]
}



##### Function that returns Root Mean Squared Error
rmse <- function(x)
{
  sqrt(mean((x^2),na.rm = TRUE))
}

rmse_list24 <- NULL

for (k in 1:ns){
  rmse_list24[[k]] <- rmse(er24[[k]])
}


write.table (x = rmse_list24, file = "rmse_eventA_sub24_simulation162p.txt", quote = FALSE)

### calculating rmse normalized
rmnor <- NULL

for (i in 1:ns){
  rmnor[i] <- rmse_list24[[i]]/ mean(q24[[zz[i,2]]])
}


write.table (x = rmnor, file = "rmnor_eventA_sub24_simulation162p.txt", quote = FALSE)



#### KGE package hydroGOF

library(hydroGOF)

kge24 <- NULL

for (i in 1:ns){
  kge24[i] <- KGE(obs = q24[[zz[i,2]]], sim = out_sim[,zz[i,1]])
}


write.table (x = kge24, file = "kge_eventA_sub24_simulation162p.txt", quote = FALSE)



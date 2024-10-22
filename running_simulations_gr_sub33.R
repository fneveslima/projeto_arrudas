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
param <- read.table(file = "27parameters.txt", header = TRUE)
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
load("iniprod27param.Rdata")
load("inirout27param.Rdata")

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
lbegin <- which(as.character(indrun)=="2013-01-06 00:00:00")
lend <- which(as.character(indrun)=="2013-01-08 23:00:00")

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
plu24 <- read.table("thiessen_plu_sub24_eventD.txt", header = TRUE)
## deleting the last line (not necessary)
plu24 <- plu24[-73,]

plu30 <- read.table("thiessen_plu_sub30_eventD.txt", header = TRUE)
## deleting the last line (not necessary)
plu30 <- plu30[-73,]

plu32 <- read.table("thiessen_plu_sub32_eventD.txt", header = TRUE)
## deleting the last line (not necessary)
plu32 <- plu32[-73,]


plu33 <- read.table("thiessen_plu_sub33_eventD.txt", header = TRUE)
## deleting the last line (not necessary)
plu33 <- plu33[-73,]



### Ponderation with areas 
#### Area of Subwatershed 24 in m²
area24 <- 48*(10^6)
#### Area of Subwatershed 30 in m²
area30 <- 54*(10^6)
#### Area of Subwatershed 32 in m²
area32 <- 23.5*(10^6)
#### Area of Subwatershed 33 in m²
area33 <- 38*(10^6)
### Area total
atotal <- area24 + area30 + area32 + area33

plup <- round(((plu24*area24)+(plu30*area30)+(plu32*area32)+(plu33*area33))/(atotal),2)

#### Extract evaporation from data input
evp_event <- evp2[lbegin:lend,]

# Initial Date
dtiA <- as.POSIXct("2013/01/06", format = "%Y/%m/%d", tz = "UTC")
# Final Date 
n_daysA <- 3 ### total of days of each event (3 days)
n_secsA <- n_daysA*24*3600
dtfA <- as.POSIXct(dtiA) + n_secsA

# Vector of date event 
dtsA <- seq(from = dtiA, to = dtfA -1, by = "3600 secs")


## Create a table with inputs of event
d33 <- data.frame(dtsA,plup,evp_event)

### Period of running model - date event
Ind_Run_33 <- seq(which(format(d33$dtsA,format="%d/%m/%Y %H:%M")=="06/01/2013 00:00"),
                  which(format(d33$dtsA,format="%d/%m/%Y %H:%M")=="08/01/2013 23:00"))


## preparation of the InputsModel object
InputsModel33 <- CreateInputsModel(FUN_MOD = RunModel_GR4H,DatesR = d33$dtsA,
                                   Precip = d33$plup,PotEvap = d33$evp_event)


### Table of contents
umide <- unlist(umid)
route <- unlist(routt)
cotr <- data.frame(param,umide,route)


## preparation of the RunOptions object and simulation
out_sim <- apply(cotr, 1, function(i){
  inrunoptions <- CreateRunOptions(FUN_MOD = RunModel_GR4H,
                                   InputsModel = InputsModel33,IndPeriod_Run = Ind_Run_33,
                                   IniResLevels = i[c("umide","route")])
  OutputsModel <- RunModel_GR4H(InputsModel = InputsModel33,RunOptions = inrunoptions, 
                                Param = c(i["X1"],i["X2"], i["X3"], i["X4"]))
  return(OutputsModel$Qsim)
})


### Rating Curve estimation using Manning Strickler equation
#### Q = alpha*(y^beta)
## reading parameters of station 33
p33 <- read.table(file = "param_station33.txt", header = TRUE)
alpha <- p33$alpha
beta <- p33$beta

## Reading water levels in event of subcatchment 33
w33 <- read.table(file = "depth_sub33hour_eventB.txt", header = TRUE)
wl33 <-  as.vector(w33$x)

### Creating a vector of flows in m3/s
Q33 <- NULL

for (i in 1:length(p33$manning)){
  Q33[[i]] <- alpha[i]*((wl33)^(beta[i]))
}



### Converting flow in m3/s in mm/h dividing per area
facm <- 3600*1000/(atotal) ### 3600 s in 1 hour and 1000 mm in 1 m

q33 <- NULL

for (k in 1:length(p33$manning)){
  q33[[k]] <- Q33[[k]]*facm
  
}

### combination of possible paths
zz <- expand.grid(1:length(param[,1]),1:length(p33$manning)) ### number of simulations x number rating curves

## number of simulations (number of rating curves x number os set parameters)
ns <- length(p33$manning)*length(param[,1])

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


nse33 <- NULL
for (i in 1:ns){
  nse33[i] <- nse(obs = q33[[zz[i,2]]], sim = out_sim[,zz[i,1]])
  
}

write.table (x = nse33, file = "nse_eventF_sub33_simulation162p.txt", quote = FALSE)


#### Calculating error
er33 <- NULL

for (i in 1:ns){
  er33[[i]] <- q33[[zz[i,2]]] - out_sim[,zz[i,1]]
}



##### Function that returns Root Mean Squared Error
rmse <- function(x)
{
  sqrt(mean((x^2),na.rm = TRUE))
}

rmse_list33 <- NULL

for (k in 1:ns){
  rmse_list33[[k]] <- rmse(er33[[k]])
}


write.table (x = rmse_list33, file = "rmse_eventF_sub33_simulation162p.txt", quote = FALSE)

### calculating rmse normalized
rmnor <- NULL

for (i in 1:ns){
  rmnor[i] <- rmse_list33[[i]]/ mean(q33[[zz[i,2]]])
}


write.table (x = rmnor, file = "rmnor_eventF_sub33_simulation162p.txt", quote = FALSE)


#### KGE package hydroGOF

library(hydroGOF)

kge33 <- NULL

for (i in 1:ns){
  kge33[i] <- KGE(obs = q33[[zz[i,2]]], sim = out_sim[,zz[i,1]] )
}


write.table (x = kge33, file = "kge_eventF_sub33_simulation162p.txt", quote = FALSE)


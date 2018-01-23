###########################################################################
# fp_vs_od.R                                                              #
# author: Pieter Coussement, Sofie Lodens                                 #
# license: GPLv3                                                          #
###########################################################################

(list = ls()) # clearing memory

# Setting Local Variables -------------------------------------------------
# Load custom functions
source("FunctionSet.R")
# Read in the data
data <- read.csv("~/20170421_gapdlenght+other_gain100/Aangepastedata.csv")
# read in the correlation between teacnOD and cuvetOD
CorLang <- read.table("~/YeastProject_SR_20150228_GrowthExperiment/FitLang.csv", quote="\"")
# set working directory
setwd('xxxx')

# Libraries ---------------------------------------------------------------
library(parallel)
library(iterators)
library(foreach)
library(doParallel)
library(ggplot2)
library(plyr)
library(reshape2)
library(grofit)
library(tidyverse)


# Preparing the data ------------------------------------------------------
# setting up for parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

data <- mutate(data, time_hour = time/3600) # time in second to hours 
data <- mutate(data, value = ifelse(signal == "ABS_600_600_NA", 
                                    CorLang[[1]] + CorLang[[2]] * value + CorLang[[3]] * value * value, 
                                    value
                                    ))

# recalculate reader OD to "cuvet OD"

# data$value[data$signal == "ABS_600_600_NA"] <-  (CorLang[[1]] + CorLang[[2]] * (data$value[data$signal == "ABS_600_600_NA"]) + CorLang[[3]] * ((data$value[data$signal == "ABS_600_600_NA"])^2))

# General OD graphs
GeneralPlotting(data, c("strain","medium"))

# Individual plotting and fitting of the correct data
data.splitWell = split(data, data$well)
temp_index = 0

for(df in 1:length(data.splitWell)){
  if ((data.splitWell[[df]])$signal[1] == "TEMP"){
    temp_index=as.numeric(df)
  }
}

temperature_df <- data.splitWell[temp_index]
data.splitWell[temp_index] <- NULL

# Start parallel processing to parse the data per well
rates <- foreach(counter=1:length(data.splitWell), .combine="rbind", .packages=c("grofit","ggplot2","doParallel")) %dopar%{
  # intialising the out dataframe
  out = data.frame(strain = character(),
                   signal = character(),
                   well = character(),
                   medium = character(),
                   parameter = character(),
                   value = numeric())
  # select the right data
  smallDataset = droplevels(data.splitWell[[counter]])
  # setting some variable for easy work afterwards
  strain = as.character(smallDataset$strain[1])
  well = as.character(smallDataset$well[1])
  medium = as.character(smallDataset$medium[1])
  # first determine the exponential phase using the OD600 values.
  growthExponentialTime = ExponentialGrowthFinalTime2(smallDataset, "ABS_600_600_NA","time_hour", 2)

  # working per signal, qx calculation
  smallDataset.split=split(smallDataset, smallDataset$signal)
  for (df in smallDataset.split){

    signal = (df$signal)[1]
    noDataSwitch = F

    if (growthExponentialTime != -1){
      #using only the data from the exponential phase
      df.exp = subset(df, time_hour < growthExponentialTime)

      #substracting the initial values
      df.exp$valueZero = df.exp$value - df.exp$value[2]
      df.exp$valueZero[df.exp$valueZero<0]=0

      med_setting = strsplit(as.character(df.exp$medium[1]), '_')[[1]][1]

      if (df.exp$medium == med_setting && df.exp$signal == "ABS_600_600_NA"){
        lowestValueIndex = which.min(df.exp$value)
        for (p in 1 : length(df.exp$value)){
          if (p < lowestValueIndex){
            df.exp$value[p] = df.exp$value[lowestValueIndex]
          }
        }
      }

      #do the fitting
      modelFit <- gcFitModel(df.exp$time_hour, df.exp$valueZero, gcID="undefinded", control=grofit.control(suppress.messages=TRUE,model.type = c("gompertz")))

      if(modelFit$fitFlag==TRUE){
        #intialize an dataOutput dataframe. Very short one, because initially only the fit parameters will be added.
        dataOutput=data.frame(parameter = character(),
                              value = numeric(),
                              stringsAsFactors = F)

        #extracting fit parameters together with their respective values.
        for(para in 1:length(modelFit$parameters)){
          for(i in 1:length(modelFit$parameters[[para]])){
            #extracting the name and the value. convert the names to something useful.
            paraName = paste(names(modelFit$parameters[para]), gsub('. ', '', names(modelFit$parameters[[para]][i])),sep="_")
            paraValue = modelFit$parameters[[para]][[i]]
            #write the info to a dataframe
            output = data.frame(parameter = as.character(paraName), value = as.numeric(paraValue))
            #add the created dataframe to an existing one.
            dataOutput = rbind(output, dataOutput)
          }
        }

        #adding additional information to the output like strain, medium, signal, well,...
        dataOutput$strain = strain
        dataOutput$signal = signal
        dataOutput$well = well
        dataOutput$medium = medium
      }else{
        #if no model could be fit
        noDataSwitch = T
      }
    }else{
      #if no exponential phase could detected
      noDataSwitch = T
    }
    if (noDataSwitch == T){
      dataOutput = data.frame(strain = strain,
                              signal = signal,
                              well = well,
                              medium = medium,
                              parameter = "mu_Estimate",
                              value = NA,
                              stringsAsFactors=FALSE
      )
    }

    #plot the data anyway.
    dir.create(file.path("./Figures/", strain), showWarnings = FALSE)

    wellSignalPlot <- ggplot(df) +
      geom_line(aes(x=time_hour, y= value)) +
      theme_bw() +
      ggtitle(paste(strain, well, medium, signal)) +
      xlab("Time [h]") +
      ylab(signal)
    ggsave(paste("./Figures/",strain,"/",well,"_",signal,".png", sep=""), width=15, height=10, units="cm")

    #data collection
    out = rbind( dataOutput, out)

    #end of per signal loop
  }


  #Calculate the fluorescence in function of the OD
  if(growthExponentialTime != -1){
    dfOdFluo = subset(smallDataset[, c("strain", "well", "time_hour","signal", "value")], time_hour < growthExponentialTime)
    # different fittings for all signals
    smallDataset$signal= factor(smallDataset$signal)
    signals= levels(smallDataset$signal)

    # extract ABS from the signals list

    for(sig in signals){
      # extracting the data
      Xval = dfOdFluo[dfOdFluo$signal == "ABS_600_600_NA",]
      Yval = dfOdFluo[dfOdFluo$signal == sig ,]

      # do the fitting
      lineair.fitting = lm(Yval$value ~ Xval$value)

      # collecting the data
      dataOutput = data.frame(strain = strain,
                              signal =  paste(as.character(sig), "vsOD", sep=""),
                              well = well,
                              medium = medium,
                              parameter = "FluoVsOD",
                              value = lineair.fitting$coefficients[[2]],
                              stringsAsFactors=FALSE
      )

      fluorescenceVsOD = ggplot() +
        geom_point(aes(x=Xval$value, y=Yval$value), color="red", size = 4) +
        geom_point(aes(x=smallDataset[smallDataset$signal == "ABS_600_600_NA","value"], y=smallDataset[smallDataset$signal == sig, "value"])) +
        ggtitle(paste(strain, medium, well, sep=" "))

      #adding the data
      out = rbind( dataOutput, out)

    }
  }

  return(out)
}
#stopping the cluster
stopCluster(cl)


#writing the raw data away
dir.create(file.path("./", "Results"), showWarnings = FALSE)
write.csv(rates, "./Results/FitResults.csv")


#making summary of the data in order to plot it.

growthRates =  subset(rates, parameter == "mu_Estimate" & signal == "ABS_600_600_NA" )
#length(growthRates$parameter)
growthRates.summ = ddply(growthRates, c("strain","medium","parameter"), summarize, N= length(value),
                         mean = mean(value, na.rm=T),
                         sd   = sd(value, na.rm=T),
                         se   = sd / sqrt(N) )

jpeg(filename= "./Figures/OverviewGrowth.jpg")
ggplot(growthRates.summ, aes(x=strain, y=mean, fill=medium)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ggtitle("Gemiddelde groeisnelheden voor Lang medium") +
  ylab("Gemiddelde groeisnelheid [h^-1]") +
  xlab("Geanalyseerde stammen") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


#Fluorescentie
fluorescentieRates =  subset(rates, parameter == "mu_Estimate" & signal == "FLUO_475_515_100" )

for (t in 1:length(fluorescentieRates$strain)){
  if (is.na(fluorescentieRates$value[t]) ){
    fluorescentieRates$value[t] = 0
  }
}

#length(growthRates$parameter)
fluorescentie.summ = ddply(fluorescentieRates, c("strain","medium","parameter"), summarize, N= length(value),
                           mean = mean(value, na.rm=T),
                           sd   = sd(value, na.rm=T),
                           se   = sd / sqrt(N) )

jpeg(filename= "./Figures/OverviewFluorescentie.jpg")
ggplot(fluorescentie.summ, aes(x=strain, y=mean, fill=medium)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ggtitle("Gemiddelde fluorescentieproductiesnelheden voor Lang medium") +
  ylab("Gemiddelde productiesnelheid [h^-1]") +
  xlab("Geanalyseerde stammen") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#Fluorescentie over OD door middel van de lineaireit tussen beide
fluorescentieVsODRates =  subset(rates, parameter == "FluoVsOD" & signal == "FLUO_475_515_100vsOD" )


#length(growthRates$parameter)
fluorescentieVsOD.summ = ddply(fluorescentieVsODRates, c("strain","medium","parameter"), summarize, N= length(value),
                               mean = mean(value, na.rm=T),
                               sd   = sd(value, na.rm=T),
                               se   = sd / sqrt(N) )

jpeg(filename= "./Figures/OverviewFluorescentieVsOD.jpg")
ggplot(fluorescentieVsOD.summ, aes(x=strain, y=mean, fill=medium)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ggtitle("Gemiddelde dGFP /dOD voor Lang medium") +
  ylab("Gemiddelde dGFP/dOD ") +
  xlab("Geanalyseerde stammen") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# Q calculation
Qdataset = dcast(rates, strain + well + medium ~ parameter + signal)[,c("strain","well","medium", "mu_Estimate_FLUO_475_515_100", "mu_Estimate_ABS_600_600_NA")]

Qdataset$Q = Qdataset$mu_Estimate_FLUO_45_515_100/Qdataset$mu_Estimate_ABS_600_600_NA

for (t in 1:length(Qdataset$strain)){
  if (is.na(Qdataset$Q[t]) & !is.na(Qdataset$mu_Estimate_ABS_600_600_NA[t]) && is.na(Qdataset$mu_Estimate_FLUO_475_515_100[t])){
    Qdataset$Q[t] = 0
  }
}

Qdataset.summ = ddply(Qdataset, c("strain","medium"), summarize, N= length(Q),
                      mean = mean(Q, na.rm=T),
                      sd   = sd(Q, na.rm=T),
                      se   = sd / sqrt(N) )

jpeg(filename= "./Figures/OverviewQ.jpg")
ggplot(Qdataset.summ, aes(x=strain, y=mean, fill=medium)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ggtitle("Gemiddelde specifiek productiesnelheden (q/µ) voor Lang medium") +
  ylab("Gemiddelde productiesnelheid q/µ") +
  xlab("Geanalyseerde stammen") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

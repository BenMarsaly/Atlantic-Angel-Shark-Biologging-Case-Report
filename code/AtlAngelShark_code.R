
setwd("~/Documents/Sciences/Publications/5_Angel shark ADL/Review 2/Github")

# import packages
library(plyr)           
library(tidyverse)
library(lubridate) 
library(circular)
library(ncdf4)
library(readr)
library(rstatix)
library(ggforce)
library(readr)

# clean environment 
rm(list = ls())



# 0 # PRELIMINARY DATA IMPORT AND FORMATING 
############################################################################

# 0.1 - Import datasets strored locally
#------------------------------------------------------------------


# import movement events summary table
Movement.Summary <- read_csv("data/Summary_Movement.csv")
# adjust time zone
Movement.Summary$Start_time <- force_tz(Movement.Summary$Start_time, "America/New_York")
Movement.Summary$End_time <- force_tz(Movement.Summary$End_time, "America/New_York")

# import stationary events summary table
Stationary.Summary <- read_csv("data/Summary_Stationary.csv")
Stationary.Summary$Start_time <- force_tz(Stationary.Summary$Start_time, "America/New_York")
Stationary.Summary$End_time <- force_tz(Stationary.Summary$End_time, "America/New_York")

# import spectrogram
Spectrogram <- read.csv2("data/Spectrogram.csv", sep = ",", dec = ".")

# import water height data from NOAA Lewes Station 8557380, include only a typical cycle for visual purpose only
Water.Height <- read_csv("data/Water.Height_NOAA8557380.csv")


# 0.2 - Import oceanographic data from Doppio ROMS using a thredds database (for section # 4 #)
#------------------------------------------------------------------


# Reference for Doppio: López AG, Wilkin JL, Levin JC (2020) Doppio - a ROMS (v3.6)-based circulation model
# for the Mid-Atlantic Bight and Gulf of Maine: configuration and comparison to integrated coastal 
# observing network observations. Geoscientific Model Development
# https://doi.org/10.5194/gmd-13-3709-2020


# 0.2.1 Extract current data from Doppio ROMS using a thredds database

library(ncdf4)

# open thredds database
f <- nc_open("https://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best")
f

# Convert the time in seconds and inform origin
tm <- ncvar_get(f, "time")*3600                         
tm <- as.POSIXct(tm, origin = "2017-11-01", tz = "UTC") 
range(tm)   # check the temporal range of the data available

# Inform the range of DateTime of the tag deployment
Start <- as.POSIXct("2023-08-11 18:00:00", tz = "UTC")  # first round hour before deployment start
End <- as.POSIXct("2023-08-15 11:00:00", tz = "UTC")    # first round hour after deployment start
tm_index <- which(tm >= Start & tm <= End)              # index of rows between start and end

# get latitude and longitude
lon <- ncvar_get(f, "lon_rho")  # extract longitude in rho coordinate system
lat <- ncvar_get(f, "lat_rho")  # extract latitude in rho coordinate system
z <- ncvar_get(f, "s_rho")      # extract vertical layers in rho coordinate system

# Find the good cell by finding the closest point to tagging location on grid
Tagging_coordinates <- c(38.857517, -75.080517)     # for tagging location
#Tagging_coordinates <- c(38.8268, -75.077)         # for tag release location (see Supplementary Materials)

# Identify the cell using lowest distance from Tagging_coordinates
dist <-sqrt( (lat - Tagging_coordinates[1])^2 + (lon - Tagging_coordinates[2])^2 ) # caclulate distances
min_index <- which(dist == min(dist), arr.ind = TRUE)     # index of the lowest distance
Lon_index <- min_index[1]                                 # index row of the cell
Lat_index <- min_index[2]                                 # index column of the cell
lon[Lon_index, Lat_index]                                 
lat[Lon_index, Lat_index]

# Extract u_eastward: dimensions are (x, y, z, t)
Current_X <- ncvar_get(f, "u_eastward",
                       start = c(Lon_index, Lat_index, 1, tm_index[1]),
                       count = c(1, 1, 1, length(tm_index)))

# Extract v_northward: dimensions are (x, y, z, t)
Current_Y <- ncvar_get(f, "v_northward", 
                       start = c(Lon_index, Lat_index, 1, tm_index[1]), 
                       count = c(1, 1, 1, length(tm_index)))

# calculate modeled current direction and magnitude
Bearing <- (atan2(Current_X, Current_Y) * 180 / pi) 

# Calculate current velocity as the norm of the vector 
Magnitude <- sqrt(Current_X^2 + Current_Y^2)

# Gather all the information into a dataframe
model.current <- data.frame(DateTime = tm[tm_index],
                            Bearing = Bearing,
                            Magnitude = Magnitude) 

# Convert to local time
model.current$DateTime <- with_tz(model.current$DateTime, tz = "America/New_York")

# Visual check bearing over time
ggplot(model.current, aes(x = DateTime, y = Bearing)) +
  geom_point() +
  geom_line() +
  theme_bw()

# Visual current velocity over time
ggplot(model.current, aes(x = DateTime, y = Magnitude)) +
  geom_point() +
  geom_line() +
  theme_bw()

# now add the current direction and velocity associated with each movement event in Movement.Summary
for (i in 1:nrow(Movement.Summary)) {
  index_current <- which((model.current$DateTime > (Movement.Summary$Start_time[i] - 30*60)) & (model.current$DateTime < (Movement.Summary$Start_time[i] + 30*60)))
  Movement.Summary$Current_direction[i] <- model.current[index_current, "Bearing"]
  Movement.Summary$Current_velocity[i] <- model.current[index_current, "Magnitude"]
}

# Same with the stationary events in Stationary.Summary
for (i in 1:nrow(Stationary.Summary)) {
  index_current <- which((model.current$DateTime > (Stationary.Summary$Start_time[i] - 30*60)) & (model.current$DateTime < (Stationary.Summary$Start_time[i] + 30*60)))
  Stationary.Summary$Current_direction[i] <- model.current[index_current, "Bearing"]
  Stationary.Summary$Current_velocity[i] <- model.current[index_current, "Magnitude"]
}



# 1 # RECOVERY PERIOD ESTIMATION
############################################################################


# 1.1 - Functions
#------------------------------------------------------------------


# 1.1.1 Splitting data into intervals and summarize


# Three functions that split data into time intervals of given durations and extract mean/count for three variables

# Function 1: Calculate the mean duration of the movements events 
#             included in each time interval (only if there is any movement)
Mean.Duration.Table <- function(Mov_Table, interval) {  # inform movement summary table and interval duration
  
  # assign a time bin to each movement event
  Time_class <- (Mov_Table$MinPostDeploy_Start %/% interval) * interval
  Table <- cbind.data.frame(Mov_Table, Time_class)                        
  
  # calculate mean duration for each time bin
  Mean.duration <- split(Table$Duration, as.factor(Table$Time_class)) %>%
    ldply(mean) %>%
    dplyr::rename(Time_class = .id, 
                  Mean_duration = V1)
  Mean.duration$Time_class <- as.numeric(Mean.duration$Time_class)    # convert to numeric
  return(Mean.duration)
}


# Function 2: Calculate the number of distinct movement events in each time interval
Mov.Freq.Table <- function(Mov_Table, interval) {       # inform movement summary table and interval duration
  
  # assign a time bin to each movement event
  Time_class <- (Mov_Table$MinPostDeploy_Start %/% interval)
  Table <- cbind.data.frame(Mov_Table, Time_class)
  
  # count number of movement event for each time bin
  Nb.mov <- table(Table$Time_class) %>%
    data.frame() %>%
    dplyr::rename(Time_class = Var1,
                  Nb_Movement = Freq)
  Nb.mov$Time_class <- as.numeric(levels(Nb.mov$Time_class))[Nb.mov$Time_class] # convert to numeric
  
  # add the time bins without any movement event, since they were not represented
  for (i in 0:max(Nb.mov$Time_class)) {
    if (!i %in% Nb.mov$Time_class) {
      Nb.mov <- add_row(Nb.mov, Time_class = i, Nb_Movement = 0, .after = i)
    }
  }
  Nb.mov$Time_class <- Nb.mov$Time_class*interval
  return(Nb.mov)
}


# Function 3: Calculate the mean tailbeat period (TBP) in each time interval
TBP.Table <- function(Mov_Table, Spectro, interval) {  # inform movement summary table, spectrogram from CWT, and interval duration
  
  # assign time bin to the spectrogram data
  Spectro$MinPostDeploy <- seq(1, nrow(Spectro))/60
  Spectro$Time_class <- Spectro$MinPostDeploy %/% interval
  
  # only keep the rows associated with movement state
  Spectro <- Spectro %>% filter(State == "Movement")
  
  Mean.TBP <- split(Spectro$DominantCycle, as.factor(Spectro$Time_class)) %>%   # split by bin
    ldply(mean) %>%                                                             # calculate mean for each bin
    dplyr::rename(Time_class = .id,                                             # rename columns
                  mean_TBP = V1)
  
  Mean.TBP$Time_class <- as.numeric(Mean.TBP$Time_class) * interval     # convert to numeric
  
  return(Mean.TBP)
}


# 1.1.2 Function for model fitting


# Three functions that  (1) fit an asymptotic regression
#                       (2) plot the data and the model
#                       (3) estimate T80% and R24h


# Function 4: Recovery model for the mean duration of movement events per interval
Mean.Duration.Recovery <- function(Table, threshold = 0.8, time = 24, interval, plot = TRUE) {
  # five arguments: - Table: data.frame with mean duration per interval (from the Mean.Duration.Table() function)
  #                 - threshold: % of change to consider the recovery. By defaut 0.8 for T80%
  #                 - time: time in hours after which we want to know the % change. By default 24 for R24h
  #                 - interval: duration of time intervals in minutes
  #                 - plot (optionnal): logical for ploting or not
  
  # non linear asymptotic regression
  Mod1 <- nls(formula = Mean_duration ~ a-(a-b)*exp(-c*Time_class),
              data = Table,
              start = list(a = 30, b = 100, c = 0.01),
              algorithm = "default", 
              trace = T,
              control = nls.control(minFactor = 0.001, maxiter= 100)) 
  
  # Calculate T80% (or TX% if different threshold))
  Recov.Threshold <-  round(-log(1-threshold)/coef(Mod1)[[3]], digit = 1)
  
  # Calculate R24h (or RXh if different time)
  Recov.Time  <- 100*round((coef(Mod1)[[2]] - (coef(Mod1)[[1]] - (coef(Mod1)[[1]]-coef(Mod1)[[2]])*exp(-coef(Mod1)[[3]]*time*60))) / (coef(Mod1)[[2]] - coef(Mod1)[[1]]), digit = 2 )
  
  if (plot == TRUE) {
    Plot1 <- ggplot(Table, aes(x = Time_class, y = Mean_duration)) +
      geom_rect(aes(xmin = 5*60, xmax = 16*60, ymin = -Inf, ymax = Inf), fill = "gray70", alpha = 0.02) +
      geom_rect(aes(xmin = 29*60, xmax = 40*60, ymin = -Inf, ymax = Inf), fill = "gray70", alpha = 0.02) +
      geom_rect(aes(xmin = 53*60, xmax = 64*60, ymin = -Inf, ymax = Inf), fill = "gray70", alpha = 0.02) +
      geom_rect(aes(xmin = 77*60, xmax = 88*60, ymin = -Inf, ymax = Inf), fill = "gray70",alpha = 0.02) +
      geom_point() +
      geom_smooth(method = "nls",
                  formula = y ~ a-(a-b)*exp(-c*x),
                  method.args = list(start = c(a = 30, b = 100, c = 0.01)),
                  se = FALSE,
                  col = "#006666") +
      geom_vline(xintercept = Recov.Threshold, col = "red", linetype = "dashed") + 
      scale_x_continuous(breaks = seq(0, 88*60, 12*60),
                         labels = seq(0, 88, 12)) +
      xlab("Hour post-release") +
      ylab("Mean duration (s)") +
      ggtitle(bquote(paste(.(interval), ' min window, T'['80%']*' = ',.(round(Recov.Threshold/60, digits = 1))," hrs, ", 'R'['24h']*' = ',.(Recov.Time),"%"))) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15))
    print(Plot1) 
  }
  return(c(Recov.Threshold, Recov.Time))
}


# Function 5: Recovery model for the number of movement events per interval (i.e. movement frequency)
Mov.Freq.Recovery <- function(Table, threshold, time, interval, plot = TRUE) {
  # five arguments: - Table: data.frame with number of movement events per interval (from the Mov.Freq.Table() function)
  #                 - threshold: % of change to consider the recovery. By defaut 0.8 for T80%
  #                 - time: time in hours after which we want to know the % change. By default 24 for R24h
  #                 - interval: duration of time intervals in minutes
  #                 - plot (optionnal): logical for ploting or not
  
  # non linear asymptotic regression
  Mod2 <- nls(formula = Nb_Movement ~ a-(a-b)*exp(-c*Time_class),
              data = Table,
              start = list(a = 5, b = 15, c = 0.001),
              algorithm = "default", trace = T) 
  
  # Calculate T80% (or TX% if different threshold))
  Recov.Threshold <-  round(-log(1-threshold)/coef(Mod2)[[3]], digit = 1)
  
  # Calculate R24h (or RXh if different time)
  Recov.Time  <- 100*round((coef(Mod2)[[2]] - (coef(Mod2)[[1]] - (coef(Mod2)[[1]]-coef(Mod2)[[2]])*exp(-coef(Mod2)[[3]]*time*60))) / (coef(Mod2)[[2]] - coef(Mod2)[[1]]), digit = 2 )
  
  if (plot == TRUE) {
    Plot2 <- ggplot(Table, aes(x = Time_class, y = Nb_Movement)) +
      geom_rect(aes(xmin = 5*60, xmax = 16*60, ymin = -Inf, ymax = Inf), fill = "gray70", alpha = 0.02) +
      geom_rect(aes(xmin = 29*60, xmax = 40*60, ymin = -Inf, ymax = Inf), fill = "gray70", alpha = 0.02) +
      geom_rect(aes(xmin = 53*60, xmax = 64*60, ymin = -Inf, ymax = Inf), fill = "gray70", alpha = 0.02) +
      geom_rect(aes(xmin = 77*60, xmax = 88*60, ymin = -Inf, ymax = Inf), fill = "gray70",alpha = 0.02) +
      geom_point() +
      geom_smooth(method = "nls",
                  formula = y ~ a-(a-b)*exp(-c*x),
                  method.args = list(start = c(a = 5, b = 15, c = 0.001)),
                  se = FALSE,
                  col = "#990066") +
      geom_vline(xintercept = Recov.Threshold, col = "red", linetype = "dashed") + 
      scale_x_continuous(breaks = seq(0, 88*60, 12*60),
                         labels = seq(0, 88, 12)) +
      xlab("Hour post-release") +
      ylab(paste("Nb of movements per", interval, "min")) +
      ggtitle(bquote(paste(.(interval), ' min window, T'['80%']*' = ',.(round(Recov.Threshold/60, digits = 1))," hrs, ", 'R'['24h']*' = ',.(Recov.Time),"%"))) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15)) 
    print(Plot2)
  }
  return(c(Recov.Threshold, Recov.Time))
}


# Function 6: Recovery model for the number of movement events per interval (i.e. movement frequency)
Mean.TBP.Recovery <- function(Table, threshold, time, interval, plot = TRUE) {
  # five arguments: - Table: data.frame with the mean tailbeat period per interval (from the TBP.Table() function)
  #                 - threshold: % of change to consider the recovery. By defaut 0.8 for T80%
  #                 - time: time in hours after which we want to know the % change. By default 24 for R24h
  #                 - interval: duration of time intervals in minutes
  #                 - plot (optionnal): logical for ploting or not
  
  # non linear asymptotic regression
  Mod3 <- nls(formula = mean_TBP ~ a-(a-b)*exp(-c*Time_class),
              data = Table,
              start = list(a = 2, b = 10, c = 0.001),
              algorithm = "default", 
              trace = T,
              control = nls.control(minFactor = 0.0001, maxiter= 100)) 
  
  # Calculate T80% (or TX% if different threshold))
  Recov.Threshold <-  round(-log(1-threshold)/coef(Mod3)[[3]], digit = 1)
  
  # Calculate R24h (or RXh if different time)
  Recov.Time  <- 100*round((coef(Mod3)[[2]] - (coef(Mod3)[[1]] - (coef(Mod3)[[1]]-coef(Mod3)[[2]])*exp(-coef(Mod3)[[3]]*time*60))) / (coef(Mod3)[[2]] - coef(Mod3)[[1]]), digit = 2 )
  
  if(plot == TRUE) {
    Plot3 <- ggplot(Table, aes(x = Time_class, y = mean_TBP)) +
      geom_rect(aes(xmin = 5*60, xmax = 16*60, ymin = -Inf, ymax = Inf), fill = "gray70", alpha = 0.02) +
      geom_rect(aes(xmin = 29*60, xmax = 40*60, ymin = -Inf, ymax = Inf), fill = "gray70", alpha = 0.02) +
      geom_rect(aes(xmin = 53*60, xmax = 64*60, ymin = -Inf, ymax = Inf), fill = "gray70", alpha = 0.02) +
      geom_rect(aes(xmin = 77*60, xmax = 88*60, ymin = -Inf, ymax = Inf), fill = "gray70",alpha = 0.02) +
      geom_point() +
      geom_smooth(method = "nls",
                  formula = y ~ a-(a-b)*exp(-c*x),
                  method.args = list(start = c(a = 2, b = 10, c = 0.001),
                                     control = nls.control(maxiter = 100)),
                  se = FALSE,
                  col = "#FF6600") +
      geom_vline(xintercept = Recov.Threshold, col = "red", linetype = "dashed") + 
      scale_x_continuous(breaks = seq(0, 88*60, 12*60),
                         labels = seq(0, 88, 12)) +
      #   ylim(c(1.2,2.0)) +
      xlab("Hour post-release") +
      ylab(paste("Mean tailbeat period (s)")) +
      ggtitle(bquote(paste(.(interval), ' min window, T'['80%']*' = ',.(round(Recov.Threshold/60, digits = 1))," hrs, ", 'R'['24h']*' = ',.(Recov.Time),"%"))) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15))
    
    print(Plot3)
  }
  return(c(Recov.Threshold, Recov.Time))
}



# 1.2 - Recovery period estimation with 60 minutes intervals
#------------------------------------------------------------------

Mean.Duration.MovEvent <- Mean.Duration.Table(Movement.Summary, interval = 60)
Mean.Duration.Recovery(Mean.Duration.MovEvent, threshold = 0.8, time = 24, interval = 60, plot = TRUE)

Mean.Freq.MovEvent <- Mov.Freq.Table(Movement.Summary, interval = 60)
Mov.Freq.Recovery(Mean.Freq.MovEvent, threshold = 0.8, time = 24, interval = 60)

Mean.TBP.MovEvent <- TBP.Table(Mov_Table = Movement.Summary,
                               Spectro = Spectrogram,
                               interval = 60)
Mean.TBP.Recovery(Mean.TBP.MovEvent, threshold = 0.8, time = 24, interval = 60)



# 1.3 Sensitivity analysis
#------------------------------------------------------------------


#Estimating recovery period using interval ranging from 15 to 180 min

# create a data.frame with one row per interval, that will incorporate T80% estimates for each variable
Sensitivity <- data.frame(Interval = seq(15,180), 
                          Duration =  integer(166), 
                          ODBA = integer(166), 
                          TBP = integer(166))

# loop to sequencially full the table with increasing interval duration
for (i in 15:180) {
  tryCatch( {             # prevents the loop to stop if the model does not converge
    
    Mean.Duration.MovEvent <- Mean.Duration.Table(Movement.Summary, interval = i)
    Sensitivity[which(Sensitivity$Interval == i), "Duration"] <- Mean.Duration.Recovery(Mean.Duration.MovEvent, 
                                                                                        threshold = 0.8, 
                                                                                        time = 24, 
                                                                                        interval = i, 
                                                                                        plot = FALSE)[1]
    
    Mean.Freq.MovEvent <- Mov.Freq.Table(Movement.Summary, interval = i)
    Sensitivity[which(Sensitivity$Interval == i), "ODBA"] <- Mov.Freq.Recovery(Mean.Freq.MovEvent, 
                                                                               threshold = 0.8, 
                                                                               time = 24, 
                                                                               interval = i, 
                                                                               plot = FALSE)[1]
    
    Mean.TBP.MovEvent <- TBP.Table(Mov_Table = , Spectro = Spectrogram, interval = i)
    Sensitivity[which(Sensitivity$Interval == i), "TBP"] <- Mean.TBP.Recovery(Mean.TBP.MovEvent, 
                                                                              threshold = 0.8, 
                                                                              time = 24, 
                                                                              interval = i, 
                                                                              plot = FALSE)[1]
  }, error = function(e) {})
}

# replace 0s by NAs (when the model did not converge)
Sensitivity[which(Sensitivity$Duration == 0), "Duration"] <- NA 
Sensitivity[which(Sensitivity$ODBA == 0), "ODBA"] <- NA 
Sensitivity[which(Sensitivity$TBP == 0), "TBC"] <- NA 

# convert the T80 estimates from minutes to hours
Sensitivity$Duration <- Sensitivity$Duration / 60
Sensitivity$ODBA <- Sensitivity$ODBA / 60
Sensitivity$TBP <- Sensitivity$TBP / 60

# plot how T80% estimates vary as a function of interval duration
ggplot(Sensitivity, aes(x = Interval, y = ODBA)) +
  geom_line(col = "#990066") +
  geom_line(aes(y = TBP), col = "#FF6600") +
  #  geom_line(aes(y = Duration), col = "#006666") + # the estimates for movement duration is not displayed because unrealistically low
  xlab("Time window (min)") +
  ylim(c(18,33))+
  ylab(bquote(paste('T'['80%'],' (hrs)'))) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 15)) 





# 2 # ANALYSIS OF THE MOVEMENT EVENTS 
############################################################################



# 2.1 - Remove the recovery period
#------------------------------------------------------------------

Recovery.Period <- 31     # inform recovery period in hours

Movement.Summary.PR <- Movement.Summary %>% filter(HourPostDeploy >= Recovery.Period)

# add an hour post recovery period
Movement.Summary.PR$HourPostRecovery <- Movement.Summary.PR$HourPostDeploy - Recovery.Period


# 2.2 - Overall movement metrics
#------------------------------------------------------------------

# How many movements?
nrow(Movement.Summary.PR)

#How long are these movements?
mean(Movement.Summary.PR$Duration) ; sd(Movement.Summary.PR$Duration)             # mean and sd
min(Movement.Summary.PR$Duration) ; max(Movement.Summary.PR$Duration)             # min and max
hist(Movement.Summary.PR$Duration,                                                # distribution
     breaks = 20,
     xlab = "Movement duration (s)")

# proportion of time spent moving
sum(Movement.Summary.PR$Duration) / (56.97639*60*60) 

# duration of the stationary events
stationary.duration <- Movement.Summary.PR$Start_time[2:nrow(Movement.Summary.PR)] - Movement.Summary.PR$End_time[1:(nrow(Movement.Summary.PR)-1)]
min(stationary.duration)  # shortest duration for a stationary event
max(stationary.duration)  # longest duration for a stationary event


# mean and sd of tailbeat period
mean(Movement.Summary.PR$Mean_TBP) # What is the avrage tailbeat period? 
sd(Movement.Summary.PR$Mean_TBP) # How does the average tailbeat period vary across movement events? 
mean(Movement.Summary.PR$SD_TBP) # What is the avrage variation of tailbeat period within movement events?

# plot of mean tailbeat period within each movement event over time
# Figure 2A:
ggplot(Movement.Summary.PR, aes(x = Start_time, Mean_TBP)) +
  geom_point() +
  ylab("Mean tailbeat period (s)") +
  xlab("") +
  geom_point(aes(x = Start_time[2], y = Mean_TBP[2]), color = "#1b9e77", size = 3) +
  geom_point(aes(x = Start_time[129], y = Mean_TBP[129]), color = "#d95f02", size = 3) +
  geom_point(aes(x = Start_time[108], y = Mean_TBP[108]), color = "#7570b3", size = 3) +
  theme_bw()


# difference in instantaneous ODBA across types of movements
Lim_TBP <- 1.1

# Long TBP
ind_long <- which(Movement.Summary.PR$Mean_TBP >= Lim_TBP & 
                    Movement.Summary.PR$Duration < 378)  # index of events with TBP>1.1
length(ind_long)                                         # how many movement events with TBP>1.1 (1 exception)
mean(Movement.Summary.PR$Duration[ind_long])             # mean duration of movement events with TBP>1.1
sd(Movement.Summary.PR$Duration[ind_long])               # sd duration of movement events with TBP>1.1
mean(Movement.Summary.PR$Average_ODBA[ind_long])         # mean ODBA of movement events with TBP>1.1
sd(Movement.Summary.PR$Average_ODBA[ind_long])           # sd ODBA of movement events with TBP>1.1
mean(Movement.Summary.PR$Max_absVV[ind_long])            # mean of max vertical velocity observed in movement events with TBP>1.1
sd(Movement.Summary.PR$Max_absVV[ind_long])              # sd of max vertical velocity observed in movement events with TBP>1.1

# short TBP
ind_short <- which(Movement.Summary.PR$Mean_TBP < Lim_TBP)  # index of events with TBP<1.1
length(ind_short)                                         # how many movement events with TBP<1.1
mean(Movement.Summary.PR$Duration[ind_short])             # mean duration of movement events with TBP<1.1
sd(Movement.Summary.PR$Duration[ind_short])               # sd duration of movement events with TBP<1.1
mean(Movement.Summary.PR$Average_ODBA[ind_short])         # mean ODBA of movement events with TBP<1.1
sd(Movement.Summary.PR$Average_ODBA[ind_short])           # sd ODBA of movement events with TBP<1.1
mean(Movement.Summary.PR$Max_absVV[ind_short])         # mean vertical velocity of movement events with TBP<1.1
sd(Movement.Summary.PR$Max_absVV[ind_short])           # sd vertical velocity of movement events with TBP<1.1



# 2.3 - Diel patterns of activity
#------------------------------------------------------------------

# 2.3.1 Diel patterns in movement duration


# table of the mean and sd duration for each hour post recovery
PR.Hourly.Mean.Duration <- split(Movement.Summary.PR$Duration, as.factor(Movement.Summary.PR$HourPostRecovery)) %>% # split the data by hour post recovery
  ldply(c(mean, sd)) %>%                                        # calculate mean for each bin
  dplyr::rename(HourPR = .id,                                   # rename columns
                Mean.Duration = V1,                             # rename columns
                SD.Duration = V2)                               # rename columns
PR.Hourly.Mean.Duration$HourPR <- as.numeric(PR.Hourly.Mean.Duration$HourPR) # convert hours to numeric
PR.Hourly.Mean.Duration$HourOfTheDay <- (as.numeric(PR.Hourly.Mean.Duration$HourPR) + Movement.Summary.PR$HourOfTheDay[1]) %% 24 # add column for hour of the day 

# table of the mean duration for each hour of the day
Hourly.Mean.Duration <- split(PR.Hourly.Mean.Duration$Mean.Duration, as.factor(PR.Hourly.Mean.Duration$HourOfTheDay)) %>% # split the data by hour of the day
  ldply(mean) %>%                                                 # calculate mean for each bin
  dplyr::rename(HourOfTheDay = .id,                               # rename column
                Mean.Duration = V1) %>%                           # rename column
  add_row(HourOfTheDay = "6", Mean.Duration = 0, .after = 6) %>%  # add 0s for hours of the day without any movements
  add_row(HourOfTheDay = "19", Mean.Duration = 0, .after = 19)
Hourly.Mean.Duration$HourOfTheDay <- as.numeric(Hourly.Mean.Duration$HourOfTheDay) # convert hours to numeric

# circular plot mean duration of movement vs hour of the day
Hourly.Mean.Duration$HourOfTheDay <- Hourly.Mean.Duration$HourOfTheDay + 0.5  # trick to align the bars correctly
# Figure 3A:
ggplot(Hourly.Mean.Duration, aes(x = HourOfTheDay, y = Mean.Duration)) +
  coord_polar() +
  geom_rect(aes(xmin = 20, xmax = 24, ymin = -Inf, ymax = Inf, alpha = 0.1), fill = "gray88", alpha = 0.02) +
  geom_rect(aes(xmin = 0, xmax = 6, ymin = -Inf, ymax = Inf), fill = "gray88", alpha = 0.02) +
  geom_col(aes(x = HourOfTheDay, y = Mean.Duration),
           col = "black", fill = "#006666",width = 1) +
  scale_x_continuous(breaks = seq(0, 23, 1),
                     labels = seq(0,23)) +
  scale_y_continuous(limits = c(0,120),
                     breaks = seq(0, 120, 30),
                     labels = seq(0, 120, 30)) +
  xlab("") +
  ylab("Mean movement duration (s)") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank())


# Test day vs night difference with Kruskall-Wallis test
kruskal.test(Duration ~ Diel, data = Movement.Summary.PR)     # test
table(Movement.Summary.PR$Diel)                               # sample size 
kruskal_effsize(Duration ~ Diel, data = Movement.Summary.PR)  # effect size


# 2.3.2 Diel patterns in movement frequency


# table number of movement per hour post recovery
Nb.mov <- table(Movement.Summary.PR$HourPostRecovery) %>%   # count occurence of hour post recovery in Movement.Summary.PR
  data.frame() %>%                        # convert to data.frame
  dplyr::rename(Time_class = Var1,        # rename column
                Nb_Movement = Freq)       # rename column
Nb.mov$Time_class <- as.numeric(levels(Nb.mov$Time_class))[Nb.mov$Time_class] # convert to numeric

# add a row for the hours during which no movement occured
for (i in 0:max(Nb.mov$Time_class)) {
  if (!i %in% Nb.mov$Time_class) {
    Nb.mov <- add_row(Nb.mov, Time_class = i, Nb_Movement = 0, .after = i)
  }
}

# add an Hour of the Day column
Nb.mov$HourOfTheDay <- (Nb.mov$Time_class + Movement.Summary.PR$HourOfTheDay[1]) %% 24

# table mean number of movement events per hour of the day
Hourly.Mean.NbMov <- split(Nb.mov$Nb_Movement, as.factor(Nb.mov$HourOfTheDay)) %>% # split data by hour of the day
  ldply(mean) %>%                         # calculate mean
  dplyr::rename(Hour = .id,               # rename column
                Mean.NumberofMov = V1)    # rename column
Hourly.Mean.NbMov$Hour <- as.numeric(Hourly.Mean.NbMov$Hour)  # convert to numeric

# circular plot mean number of movement events vs hour of the day
Hourly.Mean.NbMov$Hour <- as.numeric(Hourly.Mean.NbMov$Hour) + 0.5 # trick to align the bars correctly
# Figure 3B:
ggplot(Hourly.Mean.NbMov, aes(x = Hour, y = Mean.NumberofMov)) +
  coord_polar() +
  geom_rect(aes(xmin = 20, xmax = 24, ymin = -Inf, ymax = Inf, alpha = 0.1), fill = "gray88", alpha = 0.02) +
  geom_rect(aes(xmin = 0, xmax = 6, ymin = -Inf, ymax = Inf), fill = "gray88", alpha = 0.02) +
  geom_col(aes(x = Hour, y = Mean.NumberofMov),
           col = "black", fill = "#990066",width = 1) +
  scale_x_continuous(breaks = seq(0, 23, 1), 
                     labels = seq(0,23)) +
  ylim(c(0,12)) +
  xlab("") +
  ylab("Mean number of movement events") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank())


# Test day vs night difference with Kruskall-Wallis test
# add Diel categorical variable Day/Night to Nb.mov
Nb.mov <- Nb.mov %>% 
  mutate(Diel = case_when(HourOfTheDay >= 20 | HourOfTheDay <= 5 ~ "Night",
                          HourOfTheDay >= 6 & HourOfTheDay <= 19 ~ "Day")) 

# Kruskal-Wallis test
kruskal.test(Nb_Movement ~ Diel, data = Nb.mov)     # test
table(Nb.mov$Diel)                                  # sample size 
kruskal_effsize(Nb_Movement ~ Diel, data = Nb.mov)  # effect size


# 2.3.3 Diel patterns in ODBA


# table sum of ODBA for each hour post recovery
Hourly.Sum.ODBA <- split(Movement.Summary.PR$Sum_ODBA, as.factor(Movement.Summary.PR$HourPostRecovery)) %>% # spit the data by hour post recovery
  ldply(sum) %>%                              # calculate sum
  dplyr::rename(HourPR = .id,                 # rename column
                Sum.ODBA = V1)                # rename column
Hourly.Sum.ODBA$HourPR <- as.numeric(Hourly.Sum.ODBA$HourPR) # convert to numeric

# add 0s for the hours without any movement
for (i in 0:max(Hourly.Sum.ODBA$HourPR)) {      # for each hour post recovery
  if (!i %in% Hourly.Sum.ODBA$HourPR) {         # if i not represented
    Hourly.Sum.ODBA <- add_row(Hourly.Sum.ODBA, HourPR = i, Sum.ODBA = 0, .after = i) # add a row for hour i with a 0
  }
}

# add a hour of the day column
Hourly.Sum.ODBA$HourOfTheDay <- (Hourly.Sum.ODBA$HourPR + Movement.Summary.PR$HourOfTheDay[1]) %% 24

# table mean sum of ODBA per hour of the day
Hourly.Mean.Sum.ODBA <- split(Hourly.Sum.ODBA$Sum.ODBA, as.factor(Hourly.Sum.ODBA$HourOfTheDay)) %>% # split the data by hour of the day
  ldply(mean) %>%                             # calculate mean
  dplyr::rename(HourOfTheDay = .id,           # rename column
                Mean.Sum.ODBA = V1)           # rename column
Hourly.Mean.Sum.ODBA$HourOfTheDay <- as.numeric(Hourly.Mean.Sum.ODBA$HourOfTheDay) # convert to numeric

# circular plot mean sum ODBA vs hour of the day
Hourly.Mean.Sum.ODBA$HourOfTheDay <- Hourly.Mean.Sum.ODBA$HourOfTheDay + 0.5 # trick to align the bars correctly
# Figure 3C:
ggplot(Hourly.Mean.Sum.ODBA, aes(x = HourOfTheDay, y = Mean.Sum.ODBA)) +
  coord_polar() +
  geom_rect(aes(xmin = 20, xmax = 24, ymin = -Inf, ymax = Inf, alpha = 0.1), fill = "gray88", alpha = 0.02) +
  geom_rect(aes(xmin = 0, xmax = 6, ymin = -Inf, ymax = Inf), fill = "gray88", alpha = 0.02) +
  geom_col(aes(x = HourOfTheDay, y = Mean.Sum.ODBA),
           col = "black", fill = "#FF6600",width = 1) +
  scale_x_continuous(breaks = seq(0, 23, 1),
                     labels = seq(0,23)) +
  scale_y_continuous(limits = c(0, 141000),
                     breaks = seq(0, 140000, 35000),
                     labels = seq(0, 140000, 35000)) +
  xlab("") +
  ylab(expression('Mean sum ODBA (m.s'^-2*')')) +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank())


# Test day vs night difference with Kruskall-Wallis test
# add Diel categorical variable Day/Night to Hourly.Sum.ODBA
Hourly.Sum.ODBA <- Hourly.Sum.ODBA %>% 
  mutate(Diel = case_when(HourOfTheDay >= 20 | HourOfTheDay <= 5 ~ "Night",
                          HourOfTheDay >= 6 & HourOfTheDay <= 19 ~ "Day"))

# Kruskal-Wallis test
kruskal.test(Sum.ODBA ~ Diel, data = Hourly.Sum.ODBA)     # test
table(Hourly.Sum.ODBA$Diel)                               # sample size
kruskal_effsize(Sum.ODBA ~ Diel, data = Hourly.Sum.ODBA)  # effect size



# 2.4 - Tidal patterns of activity
#------------------------------------------------------------------


# 2.4.1 Tidal patterns in movement duration


# table of the mean duration for each hour post first higher high tide (ref = first of the deployement)
Mean.Duration.PostFirstHH <- split(Movement.Summary.PR$Duration, as.factor(Movement.Summary.PR$TimePostFirstHH)) %>% # split the data by hour post first higher high tide
  ldply(mean) %>%                                       # calculate mean for each bin
  dplyr::rename(Hour.Post.First.HH = .id,               # rename column
                mean_duration = V1)                     # rename column

# add column for hour post higher high tide 
Mean.Duration.PostFirstHH$Hour.Post.HH <- (as.numeric(Mean.Duration.PostFirstHH$Hour.Post.First.HH)) %% 25 

# table of the mean movement duration for each hour post higher high tide
Mean.Duration.Hour.PostHH <- split(Mean.Duration.PostFirstHH$mean_duration, as.factor(Mean.Duration.PostFirstHH$Hour.Post.HH)) %>% # split the data by hour post higher high tide
  ldply(mean) %>%                                       # calculte mean
  dplyr::rename(HourPostHH = .id,                       # rename column
                mean_duration = V1)                     # rename column
Mean.Duration.Hour.PostHH$HourPostHH <- as.numeric(Mean.Duration.Hour.PostHH$HourPostHH) # convert to numeric

# barplot of mean duration of movement vs hour post higher high tide
Mean.Duration.Hour.PostHH$HourPostHH <- Mean.Duration.Hour.PostHH$HourPostHH + 0.5 # trick to align the bars correctly
# Figure 3D
ggplot(Mean.Duration.Hour.PostHH, aes(x = HourPostHH, y = mean_duration)) +
  geom_col(fill = "#006666") +
  geom_line(data = Water.Height, aes(x = TimePostHH, y = (125/max(Verified))*Verified), 
            linetype = "solid", col = "black") +
  xlab("Hour post higher high tide") +
  ylab("Mean duration (s)") +
  scale_y_continuous(limits = c(0, 125),
                     breaks = seq(from = 0,to = 120, by = 40)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# 2.4.2 Tidal patterns in movement frequency


# table of the mean number of movement events for each hour post higher high tide (ref = first of the deployment)
HourPostFirstHH <- table(Movement.Summary.PR$TimePostFirstHH) %>%  #  count occurence of hour post first higher high tide in Movement.Summary.PR
  data.frame() %>%                            # convert to data.frame
  dplyr::rename(HourPostFirstHH = Var1,       # rename column
                Nb_Movement = Freq)           # rename column


# Correction for the last hour post higher high of the data that was represented by only 25 min 
# -> multiply the number of movements for this hour by 2.4
HourPostFirstHH[nrow(HourPostFirstHH), "Nb_Movement"] <- HourPostFirstHH[nrow(HourPostFirstHH), "Nb_Movement"] * 2.4

# Correction for the hour 24 that last less than 1 hour -> dive by 0.7 (or multiply by 1.4)
HourPostFirstHH[37, "Nb_Movement"] <- HourPostFirstHH[37, "Nb_Movement"]/0.7

# add column for hour post higher high tide
HourPostFirstHH$HourPostHH <- (HourPostFirstHH$HourPostFirstHH %>% 
                                 as.character() %>%
                                 as.numeric()) %% 25

# table mean number of movement events per hour post higher high tide
Mean.NbMov.Hour.PostHH <- split(HourPostFirstHH$Nb_Movement, as.factor(HourPostFirstHH$HourPostHH)) %>% # split the data by hour post higher high tide 
  ldply(mean) %>%                             # calculate mean
  dplyr::rename(HourPostHH = .id,             # rename column
                mean_number = V1)             # rename column
Mean.NbMov.Hour.PostHH$HourPostHH <- as.numeric(Mean.NbMov.Hour.PostHH$HourPostHH) # convert to numeric

# barplot of mean number of movement events vs hour post higher high tide
Mean.NbMov.Hour.PostHH$HourPostHH <- Mean.NbMov.Hour.PostHH$HourPostHH + 0.5 # trick to align the bars correctly
# Figure 3E
ggplot(Mean.NbMov.Hour.PostHH, aes(x = HourPostHH, y = mean_number)) +
  geom_col(fill = "#990066") +
  geom_line(data = Water.Height, aes(x = TimePostHH, y = (11/max(Verified))*Verified), 
            linetype = "solid", col = "black") +
  xlab("Hour from previous higher high water") +
  ylab("Mean number of movement events") +
  scale_y_continuous(limits = c(0, 11),
                     breaks = seq(0, 10, 2),
                     labels = seq(0, 10, 2)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# 2.4.3 Tidal patterns in ODBA

# table of sum ODBA for each hour post first higher high tide (ref = first of the deployment)
Sum.ODBA.PostFirstHH <- split(Movement.Summary.PR$Sum_ODBA, as.factor(floor(Movement.Summary.PR$TimePostFirstHH))) %>% # split the data by hour post first higher high tide
  ldply(sum) %>%                                # calculate sum
  dplyr::rename(Hour.Post.First.HH = .id,       # rename column
                Sum.ODBA = V1)                  # rename column

# Correction for the last hour post higher high of the data that was represented by only 25 min 
# -> multiply the number of movements for this hour by 2.4
Sum.ODBA.PostFirstHH[nrow(Sum.ODBA.PostFirstHH), "Sum.ODBA"] <- Sum.ODBA.PostFirstHH[nrow(Sum.ODBA.PostFirstHH), "Sum.ODBA"] * 2.4

# Correction for the hour 24 that last less than 1 hour -> dive by 0.7 (or multiply by 1.4)
Sum.ODBA.PostFirstHH[37, "Sum.ODBA"] <- Sum.ODBA.PostFirstHH[37, "Sum.ODBA"]/0.7

# add column for hour post higher high tide
Sum.ODBA.PostFirstHH$HourPostHH <- (as.numeric(Sum.ODBA.PostFirstHH$Hour.Post.First.HH)) %% 25

# table mean sum of ODBA for each hour post higher high tide
Mean.Sum.ODBA.PostHH <- split(Sum.ODBA.PostFirstHH$Sum.ODBA, as.factor(Sum.ODBA.PostFirstHH$HourPostHH)) %>% # split the data by hour post higher high tide 
  ldply(mean) %>%                               # calculate mean
  dplyr::rename(Hour.Post.HH = .id,             # rename column
                Mean.Sum.ODBA = V1)             # rename column
Mean.Sum.ODBA.PostHH$Hour.Post.HH <- as.numeric(Mean.Sum.ODBA.PostHH$Hour.Post.HH) # convert to numeric


# barplot of mean sum(ODBA) vs hour post higher high tide
Mean.Sum.ODBA.PostHH$Hour.Post.HH <- Mean.Sum.ODBA.PostHH$Hour.Post.HH + 0.5 # trick to align the bars correctly
# Figure 3F
ggplot(Mean.Sum.ODBA.PostHH, aes(x = Hour.Post.HH, y = Mean.Sum.ODBA)) +
  geom_col(fill = "#FF6600") +
  geom_line(data = Water.Height, aes(x = TimePostHH, y = (115000/max(Verified))*Verified), 
            linetype = "solid", col = "black") +
  xlab("Hour from previous higher high water") +
  ylab(expression('Mean sum ODBA (m.s'^-2*')')) +
  scale_y_continuous(limits = c(0, 120000),
                     breaks = seq(from = 0, to = 120000, by = 30000)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())





# 3 # BODY ORIENTATION 
############################################################################



# 3.1 Body orientation during movement events
#------------------------------------------------------------------


# 3.1.1 Movement events - Exploration


# visualize heading distribution
hist(Movement.Summary.PR$Mean_Heading, breaks = 50)

# assess variability in heading within movements
# proportion of movements that show a sd(heading)<39 
# (sd for which heading remain in a range of +-45° around the mean 75% of the time) 
length(which(Movement.Summary.PR$SD_Heading<39)) / nrow(Movement.Summary.PR)
hist(Movement.Summary.PR$SD_Heading, breaks = 50) ; abline(v = 30, col = "red", lwd = 2)

# sample size across tidal phases
table(Movement.Summary.PR$Tidal.Phase)

# Visualize difference in current direction across tidal phases
ggplot(Movement.Summary.PR, aes(x = as.factor(Tidal.Phase), y = Current_direction)) +
  geom_boxplot() +
  theme_bw()

# check how current velocity varies across tidal phases
ggplot(Movement.Summary.PR, aes(x = as.factor(Tidal.Phase), y = Current_velocity)) +
  geom_boxplot() +
  theme_bw()
# Slack tides show important variability in current direction, and generally weaker currents
# Analyses will only focus on the flood and ebb tide to avoid periods where current directions 
# and estimates from the model may not reliably describe what the animal experiences


# 3.1.2 Movement events - Orientation relative to tidal flow vs current velocity

# function to calculate the angular difference between two directions expressed in °
angle_diff <- function(a, b) {
  diff <- abs(a - b)                            # raw difference
  diff <- diff %% 360                           # wrap to 0–360
  diff <- ifelse(diff > 180, 360 - diff, diff)  # convert to 0–180
  return(diff) }

# Add a column for the angular difference between heading and current direction
Movement.Summary.PR$Angle_dif <- angle_diff(Movement.Summary.PR$Current_direction,
                                            Movement.Summary.PR$Mean_Heading)

# plot the difference between current direction and heading vs current velocity to see 
#if there is a velocity threshold after which the animal respond in terms of body orientation
# Figure 4A:
ggplot(Movement.Summary.PR, aes(x = Current_velocity, y = Angle_dif)) +
  geom_point(col = "#5ab4ac") +
  scale_y_continuous(limits = c(0, 180), 
                     breaks = seq(0, 180, 45)) +
  labs(x = "Current velocity (m/s)",
       y = "Angle Heading - Current direction (°)") +
  theme_bw(base_size = 14) 


# 3.1.3 Movement events - Circular plots with current direction and heading


# Make an flood and ebb datasets
Ebb.mov.PR <- Movement.Summary.PR %>% filter(Tidal.Phase == "Ebb")
Flood.mov.PR <- Movement.Summary.PR %>% filter(Tidal.Phase == "Flood")


# convert heading and current directions into horizontal and vertical components for circular ploting
Ebb.mov.PR <- Ebb.mov.PR %>%
  mutate(Headingx = sin(Mean_Heading * pi/180),
         Headingy = cos(Mean_Heading * pi/180),
         Currentx = sin(Current_direction * pi/180) * Current_velocity/max(Current_velocity),
         Currenty = cos(Current_direction * pi/180) * Current_velocity/max(Current_velocity))

Flood.mov.PR <- Flood.mov.PR %>%
  mutate(Headingx = sin(Mean_Heading * pi/180),
         Headingy = cos(Mean_Heading * pi/180),
         Currentx = sin(Current_direction * pi/180) * Current_velocity/max(Current_velocity),
         Currenty = cos(Current_direction * pi/180) * Current_velocity/max(Current_velocity))


# define unit circle
circle_data <- data.frame(x0 = 0, y0 = 0, r = 1)  # center at (0,0) with radius 1
r <- circle_data$r[1]

# First component of the movement plots with the arrows in the circle
# Ebb tide
p_arrows.ebb <- ggplot() +
  geom_circle(data = circle_data, aes(x0 = x0, y0 = y0, r = r), 
              color = "blue", fill = "white") +       # unit circle
  geom_segment(aes(x = -r, y = 0, xend = r, yend = 0), 
               color = "grey60", linewidth = 0.6) +   # horizontal line
  geom_segment(aes(x = 0, y = -r, xend = 0, yend = r), 
               color = "grey60", linewidth = 0.6) +   # vertical line
  geom_segment(data = Ebb.mov.PR, aes(x = 0, y = 0, xend = Headingx, yend = Headingy), 
               arrow = arrow(length = unit(0.5, "cm"))) + # plot heading 
  geom_segment(data = Ebb.mov.PR, aes(x = 0, y = 0, xend = Currentx, yend = Currenty), 
               arrow = arrow(length = unit(0.5, "cm")), col = "red")  + # plot current direction
  theme_minimal()
# visualize
p_arrows.ebb

# Flood tide
p_arrows.flood <- ggplot() +
  geom_circle(data = circle_data, aes(x0 = x0, y0 = y0, r = r),
              color = "blue", fill = "white") +       # unit circle
  geom_segment(aes(x = -r, y = 0, xend = r, yend = 0), 
               color = "grey60", linewidth = 0.6) +   # horizontal line
  geom_segment(aes(x = 0, y = -r, xend = 0, yend = r), 
               color = "grey60", linewidth = 0.6) +    # vertical line
  geom_segment(data = Flood.mov.PR, aes(x = 0, y = 0, xend = Headingx, yend = Headingy), 
               arrow = arrow(length = unit(0.5, "cm"))) + # heading arrows
  geom_segment(data = Flood.mov.PR, aes(x = 0, y = 0, xend = Currentx, yend = Currenty), 
               arrow = arrow(length = unit(0.5, "cm")), col = "red")  + # current direction arrows
  theme_minimal()
# visualize
p_arrows.flood

# Second component of the movement plots: circular barplot showing the heading circular distribution
bin.size <- 10                                            # bin size in degrees
bins <- seq(from = 0, to = 360 - bin.size, by = bin.size) # define bins
count.ebb <- rep(0, length(bins))                         # ebb empty count number of event per bin
count.flood <- rep(0, length(bins))                       # ebb empty count number of event per bin
hscale <- 0.04   # adjust to make bars taller/shorter

# Fill the count vectors
for(j in 1:length(bins)) {
  count.ebb[j] <- length(which(Ebb.mov.PR$Mean_Heading >= bins[j] & Ebb.mov.PR$Mean_Heading < bins[j] + bin.size))
  count.flood[j] <- length(which(Flood.mov.PR$Mean_Heading >= bins[j] & Flood.mov.PR$Mean_Heading < bins[j] + bin.size))
}

# Convert to table and include a unit circle buffer in the middle
Count.ebb <- data.frame(bins = bins, count.ebb) %>% 
  mutate(start = bins * pi/180,
         end   = (bins + bin.size) * pi/180,
         r0 = 1,                               # inner radius = unit circle
         r1 = 1 + count.ebb * hscale)        # outer radius

Count.flood <- data.frame(bins = bins , count.flood)  %>% 
  mutate(start = bins * pi/180,
         end   = (bins + bin.size) * pi/180,
         r0 = 1,                               # inner radius = unit circle
         r1 = 1 + count.flood * hscale)        # outer radius

# barplot with buffer for the ebb tide
ggplot() + 
  geom_arc_bar(data = Count.ebb, aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end),
               fill = "#5ab4ac",
               color = "black") +
  coord_fixed() +
  theme_void()

# barplot with buffer for the flood tide
ggplot() + 
  geom_arc_bar(data = Count.flood, aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end),
               fill = "#5ab4ac",
               color = "black") +
  coord_fixed() +
  theme_void()

# Combine arrow plot and circular barplot
# Figure 4B
p_final.ebb <- p_arrows.ebb + 
  geom_arc_bar(data = Count.ebb, aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end),
               fill = "#5ab4ac",
               color = "black") +
  coord_fixed() +
  theme_void()
p_final.ebb

# Figure 4C
p_final.flood <- p_arrows.flood + 
  geom_arc_bar(data = Count.flood, aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end),
               fill = "#5ab4ac",
               color = "black") +
  coord_fixed() +
  theme_void()
p_final.flood


# 3.1.4 Movement events - descriptive circular statistics


# ebb tide
Ebb.mov.PR$Angle_dif <- circular(Ebb.mov.PR$Angle_dif, units = "degrees")  # convert to circular
mean.circular(Ebb.mov.PR$Angle_dif)           # calculate circular mean
sd.circular(Ebb.mov.PR$Angle_dif) *180/pi      # calculate circular sd

# distribution of ebb heading
hist(Ebb.mov.PR$Mean_Heading, breaks = 50) ; abline(v = c(280,360), col = "red", lwd = 2)

# proportion of movement events with northwestward heading (280°-360°) against the current
length(which(Ebb.mov.PR$Mean_Heading > 280)) / nrow(Ebb.mov.PR)


# flood tide
Flood.mov.PR$Angle_dif <- circular(Flood.mov.PR$Angle_dif, units = "degrees")  # convert to circular
mean.circular(Flood.mov.PR$Angle_dif)           # calculate circular mean
sd.circular(Flood.mov.PR$Angle_dif) *180/pi      # calculate circular sd

# distribution of flood heading
hist(Flood.mov.PR$Mean_Heading, breaks = 50) ; abline(v = c(90,150), col = "red", lwd = 2)

# proportion of movement events with southeastward heading (90°-150°) against the current
length(which(Flood.mov.PR$Mean_Heading > 90 & Flood.mov.PR$Mean_Heading < 150 )) / nrow(Flood.mov.PR)



# 3.2 Body orientation during stationary events
#------------------------------------------------------------------


# 3.2.1 Stationary events - Remove the recovery period 


Start_PR <- as.POSIXct("2023-08-12 21:00:00", tz = "America/New_York")     # inform recovery period as DateTime

# This modifies the start_time of the stationary event overlapping the recovery time
if( any(Start_PR %within% interval(Stationary.Summary$Start_time, Stationary.Summary$End_time) == TRUE)) {
  Stationary.Summary[which(Start_PR %within% interval(Stationary.Summary$Start_time, Stationary.Summary$End_time)), "Start_time"] <- Start_PR
  warning(paste("The recovery time falls in the middle of movement #", which(Start_PR %within% interval(Stationary.Summary$Start_time, Stationary.Summary$End_time)), ", Start time has been modified to Recovery time"))
}

# Remove the recovery period
Stationary.Summary.PR <- Stationary.Summary[-which(Stationary.Summary$Start_time - Start_PR < 0),]

# add an hour post recovery period
Stationary.Summary.PR$HourPostRecovery <- difftime(Stationary.Summary.PR$Start_time, 
                                                   floor_date(Start_PR, unit = "hour"), 
                                                   units = "hours") %>% 
  as.numeric() %>% 
  floor()


# 3.2.2 Stationary events - Orientation relative to tidal flow vs current velocity


# function to calculate the angular difference between two directions expressed in ° (already defined)
angle_diff <- function(a, b) {
  diff <- abs(a - b)                            # raw difference
  diff <- diff %% 360                           # wrap to 0–360
  diff <- ifelse(diff > 180, 360 - diff, diff)  # convert to 0–180
  return(diff) }

# Add a column for the angular difference between heading and current direction
Stationary.Summary.PR$Angle_dif <- angle_diff(Stationary.Summary.PR$Current_direction,
                                              Stationary.Summary.PR$Heading)

# plot the difference between current direction and heading vs current velocity to see 
# if there is a velocity threshold after which the animal respond in terms of body orientation when stationary
# Figure 4D:
ggplot(Stationary.Summary.PR, aes(x = Current_velocity, y = Angle_dif)) +
  geom_point(col = "#d8b365") +
  scale_y_continuous(limits = c(0, 180), 
                     breaks = seq(0, 180, 45)) +
  labs(x = "Current velocity (m/s)",
       y = "Angle Heading - Current direction (°)") +
  theme_bw(base_size = 14) 


# 3.2.3 Stationary events - Circular plots with current direction and heading


# Make an flood and ebb datasets
Ebb.stat.PR <- Stationary.Summary.PR %>% filter(Tidal.Phase == "Ebb")
Flood.stat.PR <- Stationary.Summary.PR %>% filter(Tidal.Phase == "Flood")


# convert heading and current directions into horizontal and vertical components for circular ploting
Ebb.stat.PR <- Ebb.stat.PR %>%
  mutate(Headingx = sin(Heading * pi/180),
         Headingy = cos(Heading * pi/180),
         Currentx = sin(Current_direction * pi/180) * Current_velocity/max(Current_velocity),
         Currenty = cos(Current_direction * pi/180) * Current_velocity/max(Current_velocity))

Flood.stat.PR <- Flood.stat.PR %>%
  mutate(Headingx = sin(Heading * pi/180),
         Headingy = cos(Heading * pi/180),
         Currentx = sin(Current_direction * pi/180) * Current_velocity/max(Current_velocity),
         Currenty = cos(Current_direction * pi/180) * Current_velocity/max(Current_velocity))


# define unit circle (already did it for the movement events)
circle_data <- data.frame(x0 = 0, y0 = 0, r = 1)  # center at (0,0) with radius 1
r <- circle_data$r[1]

# First component of the stationary plots with the arrows in the circle
# Ebb
p_arrows.ebb <- ggplot() +
  geom_circle(data = circle_data, aes(x0 = x0, y0 = y0, r = r),
              color = "blue", fill = "white") +         # unit circle
  geom_segment(aes(x = -r, y = 0, xend = r, yend = 0), 
               color = "grey60", linewidth = 0.6) +     # horizontal line
  geom_segment(aes(x = 0, y = -r, xend = 0, yend = r), 
               color = "grey60", linewidth = 0.6) +     # vertical line
  geom_segment(data = Ebb.stat.PR, aes(x = 0, y = 0, xend = Headingx, yend = Headingy), 
               arrow = arrow(length = unit(0.5, "cm"))) +  # heading arrows
  geom_segment(data = Ebb.stat.PR, aes(x = 0, y = 0, xend = Currentx, yend = Currenty), 
               arrow = arrow(length = unit(0.5, "cm")), col = "red")  + # current direction arrows
  theme_minimal()
p_arrows.ebb


# Flood
p_arrows.flood <- ggplot() +
  geom_circle(data = circle_data, aes(x0 = x0, y0 = y0, r = r),
              color = "blue", fill = "white") +         # unit circle
  geom_segment(aes(x = -r, y = 0, xend = r, yend = 0), 
               color = "grey60", linewidth = 0.6) +     # horizontal line
  geom_segment(aes(x = 0, y = -r, xend = 0, yend = r), 
               color = "grey60", linewidth = 0.6) +     # vertical line
  geom_segment(data = Flood.stat.PR, aes(x = 0, y = 0, xend = Headingx, yend = Headingy), 
               arrow = arrow(length = unit(0.5, "cm"))) +     # heading arrows
  geom_segment(data = Flood.stat.PR, aes(x = 0, y = 0, xend = Currentx, yend = Currenty), 
               arrow = arrow(length = unit(0.5, "cm")), col = "red")  + # current direction arrows
  theme_minimal()
p_arrows.flood


# Second component of the movement plots: circular barplot showing the heading circular distribution
bin.size <- 10                                            # bin size in degrees
bins <- seq(from = 0, to = 360 - bin.size, by = bin.size) # define bins
count.ebb.2 <- rep(0, length(bins))                         # ebb empty count number of event per bin
count.flood.2 <- rep(0, length(bins))                       # ebb empty count number of event per bin
hscale <- 0.04   # adjust to make bars taller/shorter

# Fill the count vectors
for(j in 1:length(bins)) {
  count.ebb.2[j] <- length(which(Ebb.stat.PR$Heading >= bins[j] & Ebb.stat.PR$Heading < bins[j] + bin.size))
  count.flood.2[j] <- length(which(Flood.stat.PR$Heading >= bins[j] & Flood.stat.PR$Heading < bins[j] + bin.size))
}

# Convert to table and include a unit circle buffer in the middle
Count.ebb.2 <- data.frame(bins = bins, count.ebb.2) %>% 
  mutate(start = bins * pi/180,
         end   = (bins + bin.size) * pi/180,
         r0 = 1,                               # inner radius = unit circle
         r1 = 1 + count.ebb.2 * hscale)        # outer radius

Count.flood.2 <- data.frame(bins = bins , count.flood.2)  %>% 
  mutate(start = bins * pi/180,
         end   = (bins + bin.size) * pi/180,
         r0 = 1,                               # inner radius = unit circle
         r1 = 1 + count.flood.2 * hscale)        # outer radius

# barplot with buffer for the ebb tide
ggplot() + 
  geom_arc_bar(data = Count.ebb.2, aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end),
               fill = "#d8b365",
               color = "black") +
  coord_fixed() +
  theme_void()

# barplot with buffer for the flood tide
ggplot() + 
  geom_arc_bar(data = Count.flood.2, aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end),
               fill = "#d8b365",
               color = "black") +
  coord_fixed() +
  theme_void()

# Combine arrow plot and circular barplot
# Figure 4E
p_final.ebb <- p_arrows.ebb + 
  geom_arc_bar(data = Count.ebb.2, aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end),
               fill = "#d8b365",
               color = "black") +
  coord_fixed() +
  theme_void()
p_final.ebb

# Figure 4F
p_final.flood <- p_arrows.flood + 
  geom_arc_bar(data = Count.flood.2, aes(x0 = 0, y0 = 0, r0 = r0, r = r1, start = start, end = end),
               fill = "#d8b365",
               color = "black") +
  coord_fixed() +
  theme_void()
p_final.flood


# 3.2.4 Stationary events - descriptive circular statistics


# ebb tide
Ebb.stat.PR$Angle_dif <- circular(Ebb.stat.PR$Angle_dif, units = "degrees")  # convert to circular
mean.circular(Ebb.stat.PR$Angle_dif)           # calculate circular mean
sd.circular(Ebb.stat.PR$Angle_dif) *180/pi      # calculate circular sd

# distribution of ebb heading
hist(Ebb.stat.PR$Heading, breaks = 50) ; abline(v = c(290,30), col = "red", lwd = 2)

# proportion of movement events with north-northwestward heading (290°-30°) against the current
length(which(Ebb.stat.PR$Heading > 290 | Ebb.stat.PR$Heading < 30)) / nrow(Ebb.stat.PR)


# flood tide
Flood.stat.PR$Angle_dif <- circular(Flood.stat.PR$Angle_dif, units = "degrees")  # convert to circular
mean.circular(Flood.stat.PR$Angle_dif)           # calculate circular mean
sd.circular(Flood.stat.PR$Angle_dif) *180/pi      # calculate circular sd

# distribution of flood heading
hist(Flood.stat.PR$Heading, breaks = 50) ; abline(v = c(110,190), col = "red", lwd = 2)

# proportion of movement events with south-southeastward heading (100°-190°) against the current
length(which(Flood.stat.PR$Heading > 110 & Flood.stat.PR$Heading < 190 )) / nrow(Flood.stat.PR)

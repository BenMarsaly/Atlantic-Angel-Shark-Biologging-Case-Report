# Fine-scale behavior of a benthic predator: first insights into the ecology of an angel shark (Squatinidae) revealed by tri-axial biologging

This repositery contains summary data from the deployment of a biologger (accelerometer, gyroscope, magnetometer, temperature, and depth sensors) on an Atlantic angel shark (*Squatina dumeril*). This repository also includes the code for the analyses included in the publication:
> Marsaly B, Pinti J, Oliver M, Fox D, Hale E, Carlisle A (2026) Fine-scale behavior of a benthic predator: first insights into the ecology of an angel shark (Squatinidae) revealed by tri-axial biologging. *Animal Biotelemetry*.

## Data
Four data files are available in the `data/` folder: `Summary_Movement.csv`, `Summary_Stationary.csv`, `Spectrogram.csv`, and `Water.Height_NOAA8557380.csv`.

### Data summary of movement events

The file `Summary_Movement.csv` contains summary data for all the sequences classified as movement events for the whole deployment.

The table below details the variables included in the file:

| Variable | Description |
|----------|-------------|
| ID_event | Unique number ID assigned to each movement event |
| Nb_points | Number of data points included in the movement event (at 20 Hz, 20 points per seconds) |
| First_pt | Start index of event |
| Last_pt | End index of event |
| Start_time | Start Datetime of event (local time) |
| End_time | End Datetime of event (local time) |
| Duration | Event duration in seconds |
| Diel | Diel phase (day or night) based on sunrise and sunset and Start_time |
| Tidal.Phase | Tidal phase based on high and low tides from NOAA Station 8557380. High and Low tidal phases include one hour before and after the times of high and low tides |
| HourOfTheDay | Hour of the day |
| HourPostDeploy | Hour since start of the deployment |
| MinPostDeploy_Start | Event start time since start of the deployment in minutes |
| MinPostDeploy_End | Event end time since start of the deployment in minutes |
| TimePostFirstHH | Hour since a reference higher high tide (before start of the deployment) |
| TimePostHH | Hour since the last higher high tide |
| Average_ODBA | Mean instantaneous ODBA during the movement event |
| Sum_ODBA | Max instantaneous ODBA during the movement event |
| Max_absVV | Max absolute vertical velocity during the movement event |
| Mean_TBP | Mean tailbeat period during the movement event |
| SD_TBP | Standard deviation of tailbeat period during the movement event |
| Mean_Heading | Circular mean heading value during the movement event |
| SD_Heading | Circular standard deviation of heading during the movement event |


### Data summary of stationary events

The file `Summary_Stationary.csv` contains summary data for all the sequences classified as stationary events for the whole deployment.

The table below details the variables included in the file:

| Variable | Description |
|----------|-------------|
| ID_event | Unique number ID assigned to each stationary event |
| Nb_points | Number of data points included in the stationary event (at 20 Hz, 20 points per seconds) |
| First_pt | Start index of event |
| Last_pt | End index of event |
| Start_time | Start Datetime of event (local time) |
| End_time | End Datetime of event (local time) |
| Duration | Event duration in seconds |
| Diel | Diel phase (day or night) based on sunrise and sunset and Start_time |
| Tidal.Phase | Tidal phase based on high and low tides from NOAA Station 8557380. High and Low tidal phases include one hour before and after the times of high and low tides |
| HourOfTheDay | Hour of the day |
| HourPostDeploy | Hour since start of the deployment |
| Heading | Heading during the stationary event |


### Spectrogram of the tailbeat period

The file `Spectrogram.csv` is the time series resulting from the Continuous Wavelet Transform conducted on Angular Velocity in the dorsoventral axis.

The table below details the variables included in the file:

| Variable | Description |
|----------|-------------|
| DominantCycle | Dominant period in the spectrogram (highest signal strength) |
| DominantCycleAmplitude | Amplitude associated with the dominant period |
| DominantHz | Dominant frequency in the spectrogram (1/DominantCycle) |
| State | Behavioral state of the anima (Movement vs. Stationary) |


### Local water hight data 

The file `Water.Height_NOAA8557380.csv` contains a water height data from the NOAA station 8557380 located in Lewes, DE, USA. This table contains a representative water height time series over a tidal cycle on 2023-08-13, and is only use for visualization purposes.

The table below details the variables included in the file:

| Variable | Description |
|----------|-------------|
| Date | Date (UTC) |
| Time | Time (UTC) |
| DateTime | Datetime (EDT) |
| Predicted | Predicted water height |
| Predicted | Predicted water height |
| TimePostHH | Time since higher high tide |


## Code

The code used to analyze the biologging data is available in the `code/` folder. 

This script is divided in three main parts:

The first part is dedicated to the estimation of the recovery period based on the duration of movement, the frequency at which they occur, and Overall Dynamic Body Acceleration (ODBA). It also includes a sensitivity analysis regarding the temporal scale at which these three variables were integrated.

The second part focuses on the movement events post recovery, and provides overall movemeent metrics as well as investigates diel and tidal patterns in movement. 

The third part compares the animal's heading with current direction to investigate rheotaxis behavior during movement and stationary phases. Bottom current direction and velocities were estimated by a ROMS circulation model (Doppio, López et al. 2020) and accessed via a thredds database (https://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best).

This script also generates Figures 3 and 4 from the publication.



## Contact information
For any question regarding the data and code provided in this repositery, or for any additionnal data inquiries, please contact Benjamin Marsaly at marsaly.benjamin@gmail.com

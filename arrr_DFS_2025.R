#This is a series of functions for dendrochronological analysis
#By Stockton Maxwell, DFS Director
#Dendrochronology Field School
#Check out the Treeringist website for tutorial videos - https://sites.google.com/view/treeringist/dendro-help 
#For those new to R, check out Essential R - https://online.stat.psu.edu/stat484/lesson/essential-r-notes-learning-r
#For those learning dendro analysis in R, check out OpenDendro - https://opendendro.org/r/#introduction 

##############################################################################

#GETTING STARTED

#install libraries (run only once)
install.packages("dplR")
install.packages("treeclim")
install.packages("TRADER")
install.packages("DendroSync")
install.packages("dendroTools")
install.packages("burnr")
install.packages("dfoliatR")
install.packages("graphics")
install.packages("utils")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")
install.packages("tidyverse")
install.packages("utils")

#load libraries (run every R session)
library(dplR)
library(treeclim)
library(TRADER)
library(DendroSync)
library(dendroTools)
library(burnr)
library(dfoliatR)
library(graphics)
library(utils)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(utils)

#Set working directory or use Session menu
utils::choose.dir() #interactively select your directory or use the setwd() function or use the Session menu
#setwd("G:/My Drive/Research Projects/Code and Programs/rDFS")#make sure those slashes face the correct way

########################################################################################
#Dendrochronology Program Library in R 
#Cofecha type stuff in dplR - this section of the package helps the user crossdate tree ring series

#read in raw ring widths, change file name to run your own data
# Example: grow.rwl <- read.tucson(fname = "my_tree_data.rwl") #.rwl is a standard tree-ring width file format
#grow.rwl <- read.tucson(fname = "WAR_QUAL2017.rwl")
grow.rwl <- read.tucson(fname = file.choose())
rwl.stats(grow.rwl) #summary and stats of raw ring width file
seg.plot(grow.rwl) #plot of series time spans
spag.plot(grow.rwl)#spaghetti plot of raw ring widths
sens1(grow.rwl) #Calculate mean sensitivity
rwl.report(grow.rwl)

#crossdating- While dplR has crossdating functions, I highly recommend Andy Bunn's xdater app for interactive crossdating:
#https://andybunn.shinyapps.io/xDateR/
#The dplR functions below can be used for further analysis or specific tasks.
corr.rwl.seg(rwl = grow.rwl, seg.length = 50, bin.floor = 100, n = NULL, prewhiten = TRUE, pcrit = 0.05, 
             biweight = TRUE, method = c("spearman"), make.plot = TRUE, label.cex = 1, floor.plus1 = FALSE,
             master = NULL) #cofecha essentially
series.rwl.plot(grow.rwl, series = "WAR07B", series.yrs = as.numeric(names(series)), #look at an individual series
                seg.length = 50, bin.floor = 100, n = NULL,
                prewhiten = TRUE, biweight = TRUE, floor.plus1 = FALSE)
interseries.cor(grow.rwl, n = NULL, prewhiten = TRUE, biweight = TRUE, method = "spearman")#calculate interseries correlations for each series

markers <- pointer(grow.rwl)#Calculate marker rings from raw ring width series


###################################################################################
#ARSTAN stuff in dplR - this section of the package detrends or standardizes series into a site chronology
#I highly recommend Andy Bunn's iDetrend app to explore detrending:
#https://andybunn.shinyapps.io/iDetrend/
#The dplR detrending functions are shown below.

#interactive detrending - this allows you to explore curve fits for each tree ring series
grow.rwi.int <- i.detrend(rwl = grow.rwl, nyrs = NULL, f = 0.5,pos.slope = FALSE) #allows you to see a variety of fits
spag.plot(rwl = grow.rwi.int, zfac = 1, useRaster = FALSE, res = 300) #again but with the detrended series

#detrend all series at once - after you know which option is best for your data. Just adjust the method.
grow.rwi <- detrend(rwl = grow.rwl, method = c("spline"), nyrs = NULL, f = 0.5, pos.slope = FALSE) 
rwi.stats(grow.rwi) #stats for entire crn
#running stats can help you see when your common signal fades as sample size decreases
stat_out <- rwi.stats.running(grow.rwi, method = c("spearman"), prewhiten = FALSE, window.length = 50, window.overlap = 49) #running stats - time periods can be adjusted, see help

#building crn without AR model, this produces a standardized crn
grow.crn <- chron(x = grow.rwi, prefix = "XXX", biweight = TRUE, prewhiten = FALSE)
#plot crn
plot.crn(x = grow.crn, add.spline = TRUE, nyrs = 20)

#building crn with AR model, this produces a residual crn
grow.crn <- chron(x = grow.rwi, prefix = "XXX", biweight = TRUE, prewhiten = TRUE)
#plot crn
plot.crn(x = grow.crn[2:3], add.spline = TRUE, nyrs = 20)
#use sample depth cutoff if you fancy
grow.crn.trunc <- subset(grow.crn, samp.depth > 10)
plot(grow.crn.trunc[2:3],add.spline=T,nyrs=30)

#building crn using ARSTAN model that retains pooled autocorrelation
grow.crn.ars <- chron.ars(grow.rwi, biweight=TRUE, maxLag=10, firstAICmin=TRUE, 
          verbose=TRUE, prewhitenMethod=c("ar.yw","arima.CSS-ML"))
#plot crn
plot.crn(x = grow.crn.ars, add.spline = TRUE, nyrs = 20)

#building a variance stabilized crn
grow.crn.stab <- chron.stabilized(grow.rwi, winLength=101, biweight = TRUE, running.rbar = FALSE)
yrs <- time(grow.rwl)
plot.crn(x = grow.crn.stab, add.spline = TRUE, nyrs = 20)


# Wavelet transform - this allows you to look at frequencies or temporal patterns in your crn. It's good for paleoclimatology.
years <- time(grow.crn)
rings <- grow.crn[, 1]
tubular <- morlet(y1 = rings, x1 = years, p2 = NULL, dj = 0.25, siglvl = 0.95)
# --- Start of Plotting and Saving Block ---
# 1. Open the JPEG graphics device
# You can choose other formats like png(), pdf(), tiff(), etc.
# width and height are in pixels by default for jpeg/png
jpeg("waves.jpg", width = 800, height = 600, units = "px", res = 100) # res = dots per inch
# 2. Generate the plot
# The output of this will now go to "waves.jpg"
wavelet.plot(tubular, useRaster = NA)
# 3. Close the graphics device
dev.off()
# --- End of Plotting and Saving Block ---
message("Wavelet plot saved as 'waves.jpg' in your working directory.")


#####################################################################################
#Treeclim - this package allow for the assessment of growth-climate relationships

#bring in data frame from dplR, run a summary on it
summary(grow.crn)

#bring in PRISM climate data and format, change file name to run your own data
#https://prism.oregonstate.edu/explorer/; This site has monthly climate data for the USA
#clim <- read.table(file = "MTL_PRISM_ppt_tmean.csv", skip = 10, header = TRUE, sep = ",") #skips reading the header
clim <- read.table(file = file.choose(), skip = 10, header = TRUE, sep = ",") #skips reading the header
head(clim)

#---Transform PRISM data using custom function (separates Year and Month)---
PRISM.ym <- function(PRISM.data){
  PRISM.data <- PRISM.data %>% separate(Date, sep="-", into=c("Year", "Month"))
  PRISM.data$Year <- as.numeric(PRISM.data$Year)
  PRISM.data$Month <- as.numeric(PRISM.data$Month)
  return(PRISM.data)
}
#---Apply the function---
clim3 <- PRISM.ym(clim)
#Rename columns
colnames(clim3) <- (c("Year", "Month", "PPT", "TMEAN"))
head(clim3)


#Correlation & Response function analysis in treeclim - modeled after Dendroclim2002.
#Can take dynamic = "static", "moving", "evolving"
grow.crn.res <- grow.crn[,-1] #this will pull out the residual crn
grow.crn.std <- grow.crn[,-2] #this will pull out the standard crn
#this is the main function
resp <- dcc(chrono = grow.crn.std, climate = clim3, selection = -5:10,
            method = "correlation", dynamic = "static", win_size = 25, win_offset = 1, start_last = FALSE,
            timespan = NULL, var_names = NULL, ci = 0.05, boot = "stationary", sb = FALSE) #this is the main function in treeclim
resp
plot(resp) #plot the model coefficients
#traceplot(resp, variables = c("PPT.curr.jun", "PPT.curr.jul"), facet = FALSE) #shows correlations over time if moving or evolutionary selected

#More advanced climate-growth analysis
#evaluate recon skill with split calibration, requires 2 climate variables/months
#Example of using .sum() and .mean() in selection:
#selection = .sum(6:7, "PPT") + .mean(7:8, "TMEAN")  # Sum of PPT for months 6-7 and mean of TMEAN for months 7-8
calib <- dcc(chrono = grow.crn.res, climate = clim3, selection = .sum(7:8, "PPT") + .mean(7:8, "TMEAN"), #use a selection with recon variable of interest - modifiers like .mean or .sum can be used to average across months
             method = "response")
calib #show model results
plot(calib)
calib_coef <- coef(calib) #save coefficients
#this does the evaluation
skillz <- skills(object = calib, model = "ols", calibration = "50%", timespan = NULL)
plot(skillz)
skillz #show model results

#seasonal correlation - developed from Dave Meko originally in Matlab
seas <- seascorr(grow.crn.res, clim3, var_names = NULL, timespan = NULL, complete = 9,
         season_lengths = c(1, 3, 6), primary = "PPT", secondary = "TMEAN", ci = 0.05)#this is the main function
plot(seas)
seas


############################################################################################
#TRADER Code for Growth Release Detection

#run a quick summary to check your data
rwl.stats(grow.rwl) #summary and stats of raw ring width file

#Abrams and Nowacki technique - this is just one of many techniques
growthAveragingALL(grow.rwl, releases = NULL, m1 = 15, m2 = 15,buffer = 10, drawing = TRUE, criteria = 0.25, criteria2 = 0.50,prefix = "ga", gfun = mean, length = 5, storedev = jpeg)
#see output in your working directory

#plot raw data
spag.plot(grow.rwl)
raw.crn <- chron(grow.rwl, prefix = "XXX", prewhiten=FALSE, biweight = FALSE)
plot(raw.crn)

#This function calculates the synchronous growth changes (sgc), semi synchronous growth changes (ssgc) and the length of the compared overlap for a given set of tree-ring records.
changes <- sgc(grow.rwl,overlap = 50, prob = TRUE)
mean(changes$sgc_mat, na.rm = TRUE)
mean(changes$ssgc_mat, na.rm = TRUE)


##############################################################################
#Basal Area Increment calculation in dplR

#run a quick summary to check your data
rwl.stats(grow.rwl) #summary and stats of raw ring width file
basal <- bai.out(grow.rwl, diam = NULL)
spag.plot(basal[5], zfac = 1, useRaster = FALSE, res = 300)
write.csv(basal, "bai.csv")


###############################################################################

#BURNR
#How about a little fire history graphing - there's some superposed epoch analysis in there too if you fancy
#From Chris Gentry

library(burnr)
#Zion <- read_fhx('path/to/Zion.fhx') 
Zion <- read_fhx(file.choose()) #Or use file.choose()
#Sites <- read.csv('ZionSiteIDs.csv') 
Sites <- read.csv(file.choose())
facetplot <- plot_demograph(Zion, facet_group = Sites$SiteID, facet_id = Sites$series, plot_legend = TRUE)
print(facetplot)
rugplot <- plot_demograph(Zion, composite_rug = TRUE, plot_legend = TRUE)
compositerug <- rugplot + annotate('rect', xmin = 1721, xmax = 1723, ymin = 0, ymax = 21, alpha = 0.4) + annotate('rect', xmin = 1734, xmax = 1736, ymin = 0, ymax = 21, alpha = 0.4) + annotate('rect', xmin = 1748, xmax = 1750, ymin = 0, ymax = 21, alpha = 0.4) + annotate('rect', xmin = 1777, xmax = 1779, ymin = 0, ymax = 21, alpha = 0.4) + annotate('rect', xmin = 1793, xmax = 1795, ymin = 0, ymax = 21, alpha = 0.4) + scale_x_continuous(limits=c(1450, 2005), breaks = seq(1450,2005,25)) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(compositerug)


#############################################################################
#dfoliatR
#this package will allow you to detect insect outbreak events using host and non-host tree-ring chronologies
#this is an example from the Intro to doliatR vignette

#load libraries
#library(devtools) #devtools allow you to download packages not on CRAN
#install_github("chguiterman/dfoliatR") #dfoliatR is now on CRAN but you can get it from github as well
library(dfoliatR)

#load data
data(dmj_h) #host data is Douglas fir with western spruce budworm damage.
data(dmj_nh) #non-host data is ponderosa pine from a nearby control site

#analysis - uses the OUTBREAK program with defaults
dmj_defol <- defoliate_trees(host_tree = dmj_h, nonhost_chron = dmj_nh, 
                             duration_years = 8, max_reduction = -1.28, list_output = FALSE)
plot(dmj_defol)
defol_stats(dmj_defol)
dmj_obr <- outbreak(dmj_defol, filter_perc = 25, filter_min_series = 3)
plot_outbreak(dmj_obr)
dmj_obr_stats <- outbreak_stats(dmj_obr)
head(dmj_obr_stats)
dmj_interv <- diff(dmj_obr_stats$start)

# All return intervals
dmj_interv

# Mean interval
mean(dmj_interv)

# Median interval
median(dmj_interv)


####That's the end of our regularly scheduled programming.

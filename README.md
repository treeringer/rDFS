# rDFS - Dendrochonology Field School Tutorial
This is a tutorial demonstrating many of the packages created for dendrochronological analysis in R. It is meant to facilitate learning at the Dendrochronology Field School. I will update as needed.

To get started, you can download a zipped file with all of the files to run in Rstudio on your own machine. Then, open the rmarkdown document "DFS_2024.rmd". This will walk you through how to load data and run packages. Alternatively, you could run the R script and load in your own data if you are comfortable using R.

Packages include:
+**burnr**<br></br>
Steven Malevich (2019). burnr: Fire-History Analysis in R. R package version 0.3.1. https://CRAN.R-project.org/package=burnr
+**dendroTools**<br></br>
Jevsenak J. and Levanic T., 2018. dendroTools: R package for studying linear and nonlinear responses between tree-rings and daily environmental data. _Dendrochronologia_, *48*, 32-39. doi 10.1016/j.dendro.2018.01.005
Jevsenak J., Levanic T. and Dzeroski S., 2018. Comparison of an optimal regression method for climate reconstruction with the compare_methods() function from the dendroTools R package. _Dendrochronologia_, *52*, 96-104.
doi 10.1016/j.dendro.2018.10.001

+**DendroSync**<br></br>
Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios and Jordi Voltas (2019). DendroSync: A Set of Tools for Calculating Spatial Synchrony Between Tree-Ring Chronologies. R package version 0.1.3.
https://CRAN.R-project.org/package=DendroSync

+**dplR**<br></br>
Bunn AG (2008). “A dendrochronology program library in R (dplR).” _Dendrochronologia_, *26*(2), 115-124. ISSN 1125-7865, doi:10.1016/j.dendro.2008.01.002 (URL: http://doi.org/10.1016/j.dendro.2008.01.002).
Bunn AG (2010). “Statistical and visual crossdating in R using the dplR library.” _Dendrochronologia_, *28*(4), 251-258. ISSN 1125-7865, doi: 10.1016/j.dendro.2009.12.001 (URL: http://doi.org/10.1016/j.dendro.2009.12.001).
Andy Bunn, Mikko Korpela, Franco Biondi, Filipe Campelo, Pierre Mérian, Fares Qeadan, Christian Zang, Darwin, Pucha-Cofrep and Jakob Wernicke (2018). dplR: Dendrochronology Program Library in R. R package version 1.6.9.
https://CRAN.R-project.org/package=dplR
  
+**TRADER**<br></br>
Altman J, Fibich P, Dolezal J & Aakala T (2014). TRADER: a package for Tree Ring Analysis of Disturbance Events in R. _Dendrochonologia_ *32*: 107-112, URL http://www.sciencedirect.com/science/article/pii/S1125786514000058 .

**treeclim**<br></br>
Zang C, Biondi F (2015). “treeclim: an R package for the numerical calibration of proxy-climate relationships.” _Ecography_, *38*(4), 431-436. ISSN 1600-0587, doi: 10.1111/ecog.01335 (URL: http://doi.org/10.1111/ecog.01335).

Test datasets:
+ va013.rwl - Picea rubens chronology collected by Ed Cook in the 1980s near Mountain Lake, Virginia
+ Zion.fhx - Fire history collection from Chris Gentry and Peter Brown.
+ MTL_PRISM_ppt_tmean.csv - monthly precipitation and mean temperature data for Mountain Lake, Virginia. Downloaded from the PRISM Climate Group.
+ MTL_daily_PRISM_ppt_tmean.csv - daily precipitation and mean temperature data for Mountain Lake, Virginia. Dowloaded from the PRISM Climate Group.

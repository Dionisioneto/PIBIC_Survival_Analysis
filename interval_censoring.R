## ------
## Interval-censored data in R
## ------

## ---
## Libraries
## ---

if(!require(pacman)) install.packages('pacman'); library(pacman)

p_load(R2OpenBUGS, R2WinBUGS)

## Suppose that the data consist of four observations 
## (exact event time, right-censored, left-censored, and interval-censored)

## NA: censored

list(logt = c(1.0986, NA, NA, NA),
     logt.low = c(NA, 1.3863, NA, 0.6931),
     logt.upp = c(NA, NA, 0.6931, 1.0986))

















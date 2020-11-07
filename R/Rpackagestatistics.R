################Tests implemented as functions######################

#Packages
library(rio)
library(tidyverse)
library(exuber)
library(dynlm)


########################Bhargava Statistic#############################



bhargava_statistic <- function(x, tau_null = 0.2){


  #Variables for the loop
  T <- length(x)
  starting <- floor(tau_null*T)
  taus <- seq(starting, (T-starting), by = 1)
  statistics <- rep(0, length(taus))

  #For-loop

  for (i in 1:length(taus)){

    #Break Point
    break_point <- taus[i]

    #Take sample
    run_sample <- x[break_point:T]

    #Sample variance C
    sample_variance <- 1/(T-break_point) * sum(diff(run_sample)^2)

    #sample statistic
    statistic <- 1/(sample_variance * (T-break_point)^2)  * (sum((run_sample - run_sample[1])^2))
    statistics[i] <- statistic
  }

  maximum_statistic <- max(statistics)
  #names(maximum_statistic) <- "Bhargava-Statistic"
  print(maximum_statistic)
  #plot(statistics, type = "l")

}




#########################Busetti-Taylor Statistic#################################

busettitaylor_statistic <- function(x, tau_null = 0.2){

  #Initialize values
  T <- length(x)
  starting <- round(tau_null*T)
  taus <- seq(starting, (T - starting), by = 1)
  statistics <- rep(0, length(taus))


  #SChätzung der Varianz des Fehlerterms
  sample_variance <- 1/(T-1) * sum(diff(x)^2, na.rm = TRUE)

  #Berechnung der jeweiligen Statistiken
  for (i in 1:length(taus)){
    break_point <- taus[i]
    run_sample <- x[break_point:T]
    endpoint <- length(run_sample)

    #sample statistic
    statistic <- (sample_variance * (T-break_point)^2)^(-1)  * (sum((run_sample[endpoint] - run_sample)^2))
    statistics[i] <- statistic
  }


  maximum_statistic <- max(statistics)
  names(maximum_statistic) <- "Busetti-Taylor Statistic"
  print(maximum_statistic)

  data <- data.frame(taus, statistics)
  ggplot2::ggplot(data, ggplot2::aes(x = taus, y = statistics)) + ggplot2::geom_line()


}



################Kim Statistic#######################

kim_statistic <- function(x, tau_null = 0.2){

  #Initialize values
  T <- length(x) # L?nge der Zeitreihe
  starting <- round(tau_null*T) #(bei welchem starting point soll angefangen werden), beschreibt den Teil der sample, der auf explosiveness getestet wird
  taus <- seq(starting, (T-starting), 1)
  statistics <- rep(0, length(taus)) #Nullvektor, in den die Statistiken eingesetzt werden
  nominator <- rep(0, length(taus))
  denominator <- rep(0, length(taus))

  for (i in 1:length(taus)){
    break_point <- taus[i]

    #Berechnung des Zählers:
    values_nominator <- x[(break_point+1):T]
    nominator[i] <- 1/((T-break_point)^(2)) * sum((values_nominator - x[break_point])^2)

    #Berechnung des Nenners:

    values_denominator <- x[1:break_point]
    denominator[i] <- ((break_point)^(-2)) * sum((values_denominator - x[1])^2)

    #Berechnung der Teststatistik
    kim <- nominator[i]/denominator[i]

    statistics[i] <- kim
  }

  maximum_statistic <- max(statistics)
  names(maximum_statistic) <- "Kim Statistic"
  print(maximum_statistic)
 plot(statistics, type = "l")


}


########################################
#### SUPADF STATISTIC #################



supadf_statistic <- function(x, tau_null = 0.1){

  #Initiate values
  T <- length(x) # L?nge der Zeitreihe
  taus <- seq(round(tau_null * T), (T-tau_null), by = 1)
  statistics <- rep(0, length(taus)) #Nullvektor, in den die Statistiken eingesetzt werden


  #For loop
  for (i in 1:length(taus)){

    data_run <- x[1:taus[i]]

    #make a time series object
    data_run <- ts(data_run)

    #Run the regression
    DFregression <- dynlm::dynlm(diff(data_run) ~ L(data_run, 1) - 1)

    statistic <- summary(DFregression)$coefficients[1,3]
    statistics[i] <- statistic
  }

  maximum_statistic <- max(statistics)
  names(maximum_statistic) <- "Supremum ADF Statistic"
  print(maximum_statistic)
  plot(statistics, type = "l")


}



##################DFC Statistic############################

dfc_statistic <- function(x, tau_null = 0.1){



  #Initiate the for-loop
  T <- length(x)
  starting <- floor(tau_null * T)
  taus <- seq(from = starting , to = (T-starting), by = 1)
  statistics <- rep(0, length(taus))

  #Make a zero vector
  indicator_variable <- rep(0, times = T)

  #Make x a time series object for dynlm to work
  values_vector <- ts(x)


  for ( i in (1:length(taus))){
    break_point_index <- taus[i]

    #Make a zero vector for the indicator function
    indicator_variable <- rep(0, times = T)
    indicator_variable[break_point_index:T] <- 1
    indicator_variable <- ts(indicator_variable)

    #Regression
    lm <- dynlm(diff(values_vector) ~ L(values_vector,1):indicator_variable - 1)

    statistic <- summary(lm)$coefficients[1, 3]

    #Result
    statistics[i] <- statistic

  }


  maximum_statistic <- max(statistics)
  names(maximum_statistic) <- "Chow Type Dickey-Fuller Statistic"
  print(maximum_statistic)
  plot(statistics, type = "l")

}

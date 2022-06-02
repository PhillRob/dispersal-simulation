my_simulation <- function(parameters) {
  # landscape
  landscapes = ImportedLandscape(LandscapeFile = "climate_suitability.asc", 
                                 Resolution = 1, 
                                 HabPercent = TRUE, 
                                 K_or_DensDep = parameters$carrying_capacity, 
                                 SpDistFile = "news_location.asc", 
                                 SpDistResolution = 1)
  # Demography
  demos <- Demography(Rmax = parameters$max_offspring, 
                      bc = parameters$competition_coefficient, ReproductionType = 0)
  
  # Dispersal
  dist = matrix(c(parameters$short_range_distance, parameters$long_range_distance, parameters$short_range_probability), ncol = 3, byrow = T)
  disps <-  Dispersal(Emigration = Emigration(EmigProb = parameters$emigration_probability),
                      Transfer =  DispersalKernel(DoubleKernel = T, Distances = dist)) + 
    # settlement
    Settlement(MaxSteps = 0, Settle = 2) 
  
  # Initialization
  inits <- Initialise(InitType = 1, #  initialisation from a loaded species distribution map
                      SpType = 0,# all suitable cells within all distribution presence cells
                      InitDens = 2,
                      IndsHaCell = parameters$no_introduced_per_cell) 
  
  # Simulation parameters
  sims <- Simulation(Simulation = 1,
                     Years = parameters$Year,
                     Replicates = parameters$Species,
                     OutIntRange = 1,
                     OutIntPop = 1,
                     OutIntOcc = 1)
  
  # simulate
  simulate <- RSsim(seed = 13976) + landscapes + demos + disps + sims + inits
  validateRSparams(simulate) ## check parameter validity
  
  set.seed(13976)
  RunRS(RSparams = simulate, dirpath = paste0(getwd(), "/"))
  
  # summary of results
  dirpath = paste0(getwd(), "/")
  range_df = readRange(simulate, dirpath)
  
  range_1000_years = table(range_df$Rep) |> data.frame() |> rename(Replicate = Var1, Max_Years = Freq) |> mutate(Max_Years = Max_Years -1)
  range_1000_years2 = range_1000_years[order(range_1000_years$Max_Years, decreasing = T), ]
  write_delim(range_1000_years2, paste0("SummaryResults-Y", parameters$Year, "-S", parameters$Species, ".csv"), delim = ",")
  
  # lag detection
  # devtools::install_github("PhillRob/PopulationGrowthR", force = F)
  library(PopulationGrowthR)
  jpeg(paste0("Abundance-Y", parameters$Year, "-S", parameters$Species, ".png"))
  Abundance = plotAbundance(range_df, sd = T, replicates = F, main = "Abundance")
  dev.off()
  
  # data prep
  range_df <- range_df[, c(1, 2, 4)]
  colnames(range_df)<-c("Species","Year","Frequency")
  freYear <-
    aggregate(Frequency ~ Year, range_df, function(x)
      cumsum(x))
  
  fdata1 <- cbind(range_df, as.vector(unlist(freYear[2])))
  colnames(fdata1) <- c("Species", "Year", "Frequency", "Specimens")
  yeardata1 <- aggregate(Frequency ~ Year, fdata1, function(x)
    sum(x))
  
  colnames(yeardata1)<-  c("Year", "Specimens")
  
  # detection
  fitall = lagfit(data = fdata1, yeardata = yeardata1)
  lagresults<-fitall[["fitdata"]]
  
  `Mean Lag End` = max(0, mean(lagresults$Laglength, na.rm=T))
  lagspecies = table(lagresults$Lag) |> as.data.frame()
  lag_count = sum(lagspecies$Freq) - lagspecies[lagspecies$Var1 == FALSE, "Freq"]
  
  # output table for lag assessment
  parameters_table = as.tibble(parameters)
  lag_table = mutate(parameters_table, LagSpecies = lag_count, `Mean Lag End` = `Mean Lag End`, 
                     `% Lag Species` = 100* lag_count/sum(lagspecies$Freq), 
                     `% Lag Length of Years Assessed` = 100 * `Mean Lag End`/Year)
  
  if(file.exists("lag_table.csv") == TRUE) {
    lag_table_0 = read_delim("lag_table.csv", delim = ",", show_col_types = FALSE)
    lag_table = rbind(lag_table_0, lag_table)
    write_delim(lag_table, "lag_table.csv", delim = ",")
  } 
  
  write_delim(lag_table, "lag_table.csv", delim = ",")
  
  # aggregate lag by emigration probability
  aggregate(`% Lag Species` ~ emigration_probability, lag_table, FUN = mean)
}

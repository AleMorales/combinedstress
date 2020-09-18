
# Load libraries and settings ---------------------------------------------
library(SeedSorter)
library(readxl)
library(purrr)
library(furrr)
library(dplyr)
library(caret)

# Retrieve all the files and their IDs ------------------------------------


filenames = list.files(path = "Data/SeedSorting", pattern = "extra_ch0_prf\\.txt$", 
                       recursive = TRUE, full.names = FALSE, include.dirs = TRUE)

extract_name = function(x) {
  pieces = strsplit(x = x, split = "/",fixed = TRUE)[[1]]
  name = pieces[length(pieces)]
  name = strsplit(name, split = ".", fixed = TRUE)[[1]][1]
  name_pieces = strsplit(name, "[_-]",)[[1]]
  if(length(name_pieces) != 7) warning(paste("The file", x, "does not have a regular name! length =", length(name_pieces)))
  name_pieces = as.numeric(name_pieces[1:4])
  names(name_pieces) = c("ID", "Block", "Sample", "Rep")
  name_pieces
}

sample_information = plyr::ldply(as.list(filenames), extract_name) %>% mutate(filename = filenames)

# Compare to fitness file
fitness = read_excel("Data/Fitness.xlsx")
fitness = mutate(fitness, Sample = ID)

sample_information = left_join(sample_information, select(fitness, Block, Sample, Genotype, Treatment), 
                               by = c("Block", "Sample"))

sample_information = filter(sample_information, !is.na(Genotype))
N = nrow(sample_information)

data =  map(1:nrow(sample_information), 
     function(x) {
       filename = sample_information[x, "filename"]
       profile_file = gsub(".txt", ".fst",sample_information[x, "filename"])
       main_file = gsub("_ch0_prf", "",sample_information[x, "filename"])
       getPredictionData(main_file = main_file, profile_file = profile_file, datadir = paste0("Data/SeedSorting"))
     })


# For each filename, classify the class of particle -----------------------

classifiers = readRDS("Data/modelsCoarse.rds")
names(classifiers) = c("An-1", "Bay-0", "Bur-0", "Col-0", "Lp2-6")


# lda ---------------------------------------------------------------------
lda_data = tibble(nSeeds = rep(NA, N),
                  nNonSeeds = rep(NA, N),
                  medSize = rep(NA, N),
                  madSize = rep(NA, N))
for(i in 1:nrow(sample_information)) {
  Genotype = sample_information[i,"Genotype"]
  temp = classifySeeds(classifiers[[Genotype]][["lda"]], as.data.frame(data[[i]][,1:7]))$data$response
  data[[i]] = mutate(data[[i]], ldaClass = temp)
  temp = filter(data[[i]], ldaClass == "S") %>% 
    summarise(nSeeds = n(), medSize = median(Size,na.rm = TRUE), madSize = mad(Size,na.rm = TRUE))
  lda_data[i,-2] = select(temp, nSeeds, medSize, madSize)
  lda_data[i,2] = nrow(data[[i]]) - temp$nSeeds[[1]]
}
sample_information = mutate(sample_information, 
                            nSeeds_lda = lda_data$nSeeds, 
                            nNonSeeds_lda = lda_data$nNonSeeds,
                            medSize_lda = lda_data$medSize,
                            madSize_lda = lda_data$madSize)

# qda ---------------------------------------------------------------------
qda_data = tibble(nSeeds = rep(NA, N),
                  nNonSeeds = rep(NA, N),
                  medSize = rep(NA, N),
                  madSize = rep(NA, N))
for(i in 1:nrow(sample_information)) {
  Genotype = sample_information[i,"Genotype"]
  temp = classifySeeds(classifiers[[Genotype]][["qda"]], as.data.frame(data[[i]][,1:7]))$data$response
  data[[i]] = mutate(data[[i]], qdaClass = temp)
  temp = dplyr::filter(data[[i]], qdaClass == "S") %>% 
    summarise(nSeeds = n(), medSize = median(Size), madSize = mad(Size))
  qda_data[i,-2] = select(temp, nSeeds, medSize, madSize)
  qda_data[i,2] = nrow(data[[i]]) - temp$nSeeds[[1]]
}
sample_information = mutate(sample_information, 
                            nSeeds_qda = qda_data$nSeeds, 
                            nNonSeeds_qda = qda_data$nNonSeeds,
                            medSize_qda = qda_data$medSize,
                            madSize_qda = qda_data$madSize)

# extinction ---------------------------------------------------------------------
extinction_data = tibble(nSeeds = rep(NA, N),
                         nNonSeeds = rep(NA, N),
                         medSize = rep(NA, N),
                         madSize = rep(NA, N))
for(i in 1:nrow(sample_information)) {
  Genotype = sample_information[i,"Genotype"]
  temp = classifySeeds(classifiers[[Genotype]][["extinction"]], as.data.frame(data[[i]][,1:7]))$data$response
  data[[i]] = mutate(data[[i]], extinctionClass = temp)
  temp = dplyr::filter(data[[i]], extinctionClass == "S") %>% 
    summarise(nSeeds = n(), medSize = median(Size), madSize = mad(Size))
  extinction_data[i,-2] = select(temp, nSeeds, medSize, madSize)
  extinction_data[i,2] = nrow(data[[i]]) - temp$nSeeds[[1]]
}
sample_information = mutate(sample_information, 
                            nSeeds_extinction = extinction_data$nSeeds, 
                            nNonSeeds_extinction = extinction_data$nNonSeeds,
                            medSize_extinction = extinction_data$medSize,
                            madSize_extinction = extinction_data$madSize)


# knn ---------------------------------------------------------------------
knn_data = tibble(nSeeds = rep(NA, N),
                  nNonSeeds = rep(NA, N),
                  medSize = rep(NA, N),
                  madSize = rep(NA, N))
for(i in 1:nrow(sample_information)) {
  Genotype = sample_information[i,"Genotype"]
  temp = classifySeeds(classifiers[[Genotype]][["knn"]], as.data.frame(data[[i]][,1:7]))$data$response
  data[[i]] = mutate(data[[i]], knnClass = temp)
  temp = dplyr::filter(data[[i]], knnClass == "S") %>% 
    summarise(nSeeds = n(), medSize = median(Size), madSize = mad(Size))
  knn_data[i,-2] = select(temp, nSeeds, medSize, madSize)
  knn_data[i,2] = nrow(data[[i]]) - temp$nSeeds[[1]]
}
sample_information = mutate(sample_information, 
                            nSeeds_knn = knn_data$nSeeds, 
                            nNonSeeds_knn = knn_data$nNonSeeds,
                            medSize_knn = knn_data$medSize,
                            madSize_knn = knn_data$madSize)

# naiveBayes ---------------------------------------------------------------------
naiveBayes_data = tibble(nSeeds = rep(NA, N),
                         nNonSeeds = rep(NA, N),
                         medSize = rep(NA, N),
                         madSize = rep(NA, N))
for(i in 1:nrow(sample_information)) {
  Genotype = sample_information[i,"Genotype"]
  temp = classifySeeds(classifiers[[Genotype]][["naiveBayes"]], as.data.frame(data[[i]][,1:7]))$data$response
  data[[i]] = mutate(data[[i]], naiveBayesClass = temp)
  temp = dplyr::filter(data[[i]], naiveBayesClass == "S") %>% 
    summarise(nSeeds = n(), medSize = median(Size), madSize = mad(Size))
  naiveBayes_data[i,-2] = select(temp, nSeeds, medSize, madSize)
  naiveBayes_data[i,2] = nrow(data[[i]]) - temp$nSeeds[[1]]
}
sample_information = mutate(sample_information, 
                            nSeeds_naiveBayes = naiveBayes_data$nSeeds,
                            nNonSeeds_naiveBayes = naiveBayes_data$nNonSeeds,
                            medSize_naiveBayes = naiveBayes_data$medSize,
                            madSize_naiveBayes = naiveBayes_data$madSize)


# logistic ---------------------------------------------------------------------
logistic_data = tibble(nSeeds = rep(NA, N),
                       nNonSeeds = rep(NA, N),
                       medSize = rep(NA, N),
                       madSize = rep(NA, N))
for(i in 1:nrow(sample_information)) {
  Genotype = sample_information[i,"Genotype"]
  temp = classifySeeds(classifiers[[Genotype]][["logistic"]], as.data.frame(data[[i]][,1:7]))$data$response
  data[[i]] = mutate(data[[i]], logisticClass = temp)
  temp = dplyr::filter(data[[i]], logisticClass == "S") %>% 
    summarise(nSeeds = n(), medSize = median(Size), madSize = mad(Size))
  logistic_data[i,-2] = select(temp, nSeeds, medSize, madSize)
  logistic_data[i,2] = nrow(data[[i]]) - temp$nSeeds[[1]]
}
sample_information = mutate(sample_information, 
                            nSeeds_logistic = logistic_data$nSeeds, 
                            nNonSeeds_logistic = logistic_data$nNonSeeds,
                            medSize_logistic = logistic_data$medSize,
                            madSize_logistic = logistic_data$madSize)


# xgboost ---------------------------------------------------------------------
xgboost_data = tibble(nSeeds = rep(NA, N),
                      nNonSeeds = rep(NA, N),
                      medSize = rep(NA, N),
                      madSize = rep(NA, N))
for(i in 1:nrow(sample_information)) {
  Genotype = sample_information[i,"Genotype"]
  temp = classifySeeds(classifiers[[Genotype]][["xgboost"]], as.data.frame(data[[i]][,1:7]))$data$response
  data[[i]] = mutate(data[[i]], xgboostClass = temp)
  temp = dplyr::filter(data[[i]], xgboostClass == "S") %>% 
    summarise(nSeeds = n(), medSize = median(Size), madSize = mad(Size))
  xgboost_data[i,-2] = select(temp, nSeeds, medSize, madSize)
  xgboost_data[i,2] = nrow(data[[i]]) - temp$nSeeds[[1]]
}
sample_information = mutate(sample_information, 
                            nSeeds_xgboost = xgboost_data$nSeeds, 
                            nNonSeeds_xgboost = xgboost_data$nNonSeeds,
                            medSize_xgboost = xgboost_data$medSize,
                            madSize_xgboost = xgboost_data$madSize)


# Calculate median predictions for each type of variable ------------------

sample_information = mutate(sample_information, 
                            nSeeds = plyr::aaply(cbind(nSeeds_lda, nSeeds_qda, nSeeds_knn, nSeeds_extinction, 
                                                       nSeeds_naiveBayes, nSeeds_logistic, nSeeds_xgboost), 1, median),
                            
                            nNonSeeds = plyr::aaply(cbind(nNonSeeds_lda, nNonSeeds_qda, nNonSeeds_knn, nNonSeeds_extinction, 
                                                       nNonSeeds_naiveBayes, nNonSeeds_logistic, nNonSeeds_xgboost), 1, median),
                            
                            medSize = plyr::aaply(cbind(medSize_lda, medSize_qda, medSize_knn, medSize_extinction, 
                                                        medSize_naiveBayes, medSize_logistic, medSize_xgboost), 1, median),
                            
                            madSize = plyr::aaply(cbind(madSize_lda, madSize_qda, madSize_knn, madSize_extinction, 
                                                        madSize_naiveBayes, madSize_logistic, madSize_xgboost), 1, median))


# Rename the table to match other trait information

treatment = c(C = "C", Ht = "HT", D1 = "D1", D2 = "D2", HtD = "HTD", Flww = "S", FlD = "PSD")
  
sample_information = mutate(sample_information, 
                            TreatmentOld = Treatment, 
                            Treatment = treatment[TreatmentOld],
                            XP = "Fitness")

# Save the data for later
saveRDS(file = "Intermediate/SeedSortingAll.rds", object = sample_information)


# Modify the data in preparation for trait analysis
seed_sorting = select(sample_information, XP, Block, Rep, Treatment, Genotype, nSeeds, medSize, madSize)
seed_sorting = group_by(sample_information, XP, Block, Rep, Treatment, Genotype) %>% 
  summarise(medSize = mean(medSize, na.rm = TRUE),
            madSize = mean(madSize, na.rm = TRUE), 
            nSeeds = mean(nSeeds, na.rm = TRUE),
            nNonSeeds = mean(nNonSeeds, na.rm = TRUE)) %>% ungroup()


saveRDS(file = "Intermediate/SeedSorting.rds", object = seed_sorting)

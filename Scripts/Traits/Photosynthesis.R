# 01 Load packages and data -----------------------------------------------
library(tidyverse)
library(magrittr)


# Process Photosynthesis Data ---------------------------------------------

# Column names and folders
folder = "Data/LiCOR/"
colnames = as.character(read_delim(paste0(folder,"/f2/2018-04-16-0953_F2-21-C2-27"),
                                   skip = 55, n_max = 1, delim = "\t", col_names = FALSE,
                                   col_types = cols(.default = col_character(), X82 = col_logical())))


# Function to load the data
readLi6800 = function(file) {
  temp = suppressMessages(suppressWarnings(read_delim(
    paste0(folder,file), skip = 57,col_names = colnames, delim = "\t")))
  temp$date = as.POSIXlt(temp$date,format="%Y%m%d %H:%M:%S")
  temp$TOD = sapply(temp$date, function(x) x$hour + x$min/60 + x$sec/3600)
  list(temp)
}


# 9 files from F2
Li6800 =     tibble(Batch = "F2", GT = 21, Treatment = "C", MT = 27, Replicate = 1, SWC = 0.8, SLA = 598, N = 6.23, Chl = 95.05, rawdata = readLi6800("/f2/2018-04-16-0953_F2-21-C2-27"))
Li6800 %<>% add_row(Batch = "F2", GT = 21, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 598, N = 6.23, Chl = 95.05, rawdata = readLi6800("/f2/2018-04-16-1345_f2-21-c2-21"))
Li6800 %<>% add_row(Batch = "F2", GT = 21, Treatment = "Su7", MT = 21, Replicate = 1, SWC = 0.8, SLA = 545, N = 7.07, Chl = 114.17, rawdata = readLi6800("/f2/2018-04-17-0914_f2-21-Su"))
Li6800 %<>% add_row(Batch = "F2", GT = 21, Treatment = "PSD7", MT = 21, Replicate = 1, SWC = 0.31, SLA = 505, N = 3.58, Chl = 166.68, rawdata = readLi6800("/f2/2018-04-17-1139_f2-21-psd2"))
Li6800 %<>% add_row(Batch = "F2", GT = 21, Treatment = "C", MT = 21, Replicate = 2, SWC = 0.8, SLA = 580, N = 5.58, Chl = 164.76, rawdata = readLi6800("/f2/2018-04-17-1402_f2-21-c3"))
Li6800 %<>% add_row(Batch = "F2", GT = 21, Treatment = "D7", MT = 21, Replicate = 1, SWC = 0.09, SLA = 414, N = 5.34, Chl = 220.65, rawdata = readLi6800("/f2/2018-04-17-1607_f2-21-d2"))
Li6800 %<>% add_row(Batch = "F2", GT = 27, Treatment = "C", MT = 27, Replicate = 1, SWC = 0.8, SLA = 636, N = 5.44, Chl = 78.96, rawdata = readLi6800("/f2/2018-04-18-0911_f2-27-ht-27"))
Li6800 %<>% add_row(Batch = "F2", GT = 27, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 636, N = 5.44, Chl = 78.96, rawdata = readLi6800("/f2/2018-04-18-1357_f2-27-ht-21"))
Li6800 %<>% add_row(Batch = "F2", GT = 27, Treatment = "D7", MT = 27, Replicate = 1, SWC = 0.07, SLA = 557, N = 6.26, Chl = 143.26, rawdata = readLi6800("/f2/2018-04-18-1136_f2-27-htd"))

# 9 files from F3
Li6800 %<>% add_row(Batch = "F3", GT = 27, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 591, N = 5.70, Chl = 162.00, rawdata = readLi6800("/f3/2018-04-25-1506_f3-27-ht-21"))
# The following files removed because matching went wrong
#Li6800 %<>% add_row(Batch = "F3", GT = 27, Treatment = "D7", MT = 27, Replicate = 1, rawdata = readLi6800("/f3/2018-04-25-1217_f3-27-htd"))
#Li6800 %<>% add_row(Batch = "F3", GT = 27, Treatment = "C", MT = 27, Replicate = 1, rawdata = readLi6800("/f3/2018-04-25-0918_f3-27-ht-27"))
Li6800 %<>% add_row(Batch = "F3", GT = 21, Treatment = "C", MT = 27, Replicate = 1, SWC = 0.8, SLA = 591, N = 5.94, Chl = 182.22, rawdata = readLi6800("/f3/2018-04-24-1343_f3-21-c3-27"))
Li6800 %<>% add_row(Batch = "F3", GT = 21, Treatment = "D7", MT = 21, Replicate = 1, SWC = 0.18, SLA = 491, N = 6.01, Chl = 227.83, rawdata = readLi6800("/f3/2018-04-24-1133_f3-21-d2"))
Li6800 %<>% add_row(Batch = "F3", GT = 21, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 502, N = 5.94, Chl = 182.22, rawdata = readLi6800("/f3/2018-04-24-0949_f3-21-c3-21"))
Li6800 %<>% add_row(Batch = "F3", GT = 21, Treatment = "Su7", MT = 21, Replicate = 1, SWC = 0.8, SLA = 586, N = 4.32, Chl = 133.74, rawdata = readLi6800("/f3/2018-04-23-1224_f3-21-su"))
Li6800 %<>% add_row(Batch = "F3", GT = 21, Treatment = "PSD7", MT = 21, Replicate = 1, SWC = 0.8, SLA = 663, N = 4.55, Chl = 129.10, rawdata = readLi6800("/f3/2018-04-23-0926_f3-21-psd2"))
Li6800 %<>% add_row(Batch = "F3", GT = 21, Treatment = "C", MT = 27, Replicate = 2, SWC = 0.8, SLA = 573, N = 6.48, Chl = 173.38, rawdata = readLi6800("/f3/2018-04-21-1352_f3-21-c2-27"))
Li6800 %<>% add_row(Batch = "F3", GT = 21, Treatment = "C", MT = 21, Replicate = 2, SWC = 0.8, SLA = 573, N = 6.48, Chl = 173.38, rawdata = readLi6800("/f3/2018-04-21-0927_f3-21-c2-21"))
Li6800 %<>% add_row(Batch = "F3", GT = 21, Treatment = "Su1", MT = 21, Replicate = 1, SWC = 0.36, SLA = 633, N = 6.10, Chl = 118.22, rawdata = readLi6800("/f3/2018-04-21-1158_f3-21-psd1"))

# 15 files from F4
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "C", MT = 27, Replicate = 1, SWC = 0.8, SLA = 560, N = 5.58, Chl = 182.97, rawdata = readLi6800("/f4/2018-04-26-0919_f4-21-c2-27"))
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "C", MT = 27, Replicate = 2, SWC = 0.8, SLA = 555, N = 5.57, Chl = 197.94, rawdata = readLi6800("/f4/2018-04-30-1422_f4-21-c3-27"))
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "Su1", MT = 21, Replicate = 1, SWC = 0.8, SLA = 676, N = 5.95, Chl = 123.98, rawdata = readLi6800("/f4/2018-04-26-1151_f4-21-psd1"))
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 560, N = 5.58, Chl = 182.97, rawdata = readLi6800("/f4/2018-04-26-1402_f4-21-c2-21"))
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "C", MT = 21, Replicate = 2, SWC = 0.8, SLA = 555, N = 5.57, Chl = 197.94, rawdata = readLi6800("/f4/2018-04-30-0959_f4-21-c3-21"))
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "D4", MT = 21, Replicate = 1, SWC = 0.36, SLA = 581, N = 6.51, Chl = 201.20, rawdata = readLi6800("/f4/2018-04-27-0918_f4-21-d4r1"))
Li6800 %<>% add_row(Batch = "F4", GT = 27, Treatment = "D4", MT = 27, Replicate = 1, SWC = 0.36, SLA = 594, N = 5.37, Chl = 82.81, rawdata = readLi6800("/f4/2018-04-27-1245_f4-27-htd4-r1"))
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "D4", MT = 21, Replicate = 1, SWC = 0.36, SLA = 489, N = 5.41, Chl = 214.19, rawdata = readLi6800("/f4/2018-04-29-1322_f4-21-d4-r2"))
Li6800 %<>% add_row(Batch = "F4", GT = 27, Treatment = "D4", MT = 27, Replicate = 1, SWC = 0.29, SLA = NA, N = 4.93, Chl = 174.12, rawdata = readLi6800("/f4/2018-04-29-0951-f4-27-htd4-r2"))
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "D7", MT = 21, Replicate = 1, SWC = 0.22, SLA = 488, N = 4.78, Chl = 221.71, rawdata = readLi6800("/f4/2018-04-30-1213_f4-21-d2"))
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "PSD7", MT = 21, Replicate = 1, SWC = 0.28, SLA = 516, N = 3.51, Chl = 135.90, rawdata = readLi6800("/f4/2018-05-01-0924_f4-21-psd2"))
Li6800 %<>% add_row(Batch = "F4", GT = 21, Treatment = "Su7", MT = 21, Replicate = 1, SWC = 0.8, SLA = 613, N = 3.62, Chl = 166.64, rawdata = readLi6800("/f4/2018-05-01-1149_f4-21-su"))
Li6800 %<>% add_row(Batch = "F4", GT = 27, Treatment = "C", MT = 27, Replicate = 1, SWC = 0.8, SLA = 645, N = 5.25, Chl = 155.00, rawdata = readLi6800("/f4/2018-05-02-0921_f4-27-ht-27"))
Li6800 %<>% add_row(Batch = "F4", GT = 27, Treatment = "D7", MT = 27, Replicate = 1, SWC = 0.12, SLA = 522, N = 5.64, Chl = 116.07, rawdata = readLi6800("/f4/2018-05-02-1215_f4-27-htd"))
Li6800 %<>% add_row(Batch = "F4", GT = 27, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 645, N = 5.25, Chl = 155.00, rawdata = readLi6800("/f4/2018-05-02-1450_f4-27-ht-21"))

# 19 files from F5
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 586, N = 6.05, Chl = 195.80, rawdata = readLi6800("/f5/2018-05-03-0922_f5-21-c2-21"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "Su1", MT = 21, Replicate = 1, SWC = 0.8, SLA = 561, N = 5.58, Chl = 140.74, rawdata = readLi6800("/f5/2018-05-03-1233_f5-psd1"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "C", MT = 27, Replicate = 1, SWC = 0.8, SLA = 586, N = 6.05, Chl = NA, rawdata = readLi6800("/f5/2018-05-03-1445_f5-21-c2-27"))
Li6800 %<>% add_row(Batch = "F5", GT = 27, Treatment = "D4", MT = 27, Replicate = 1, SWC = 0.36, SLA = 568, N = 5.95, Chl = 73.30, rawdata = readLi6800("/f5/2018-05-04-0920_f5-27-htd4"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "D4", MT = 21, Replicate = 1, SWC = 0.37, SLA = 591, N = 6.32, Chl = 182.14, rawdata = readLi6800("/f5/2018-05-04-1208_f5-21-d4"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "D4", MT = 21, Replicate = 2, SWC = 0.38, SLA = 578, N = 5.71, Chl = 172.30, rawdata = readLi6800("/f5/2018-05-05-1001_f5-21-d4r2"))
Li6800 %<>% add_row(Batch = "F5", GT = 27, Treatment = "D4", MT = 27, Replicate = 2, SWC = 0.29, SLA = 705, N = 5.28, Chl = 153.02, rawdata = readLi6800("/f5/2018-05-05-1310_F5-htd4-r2"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "PSD4", MT = 21, Replicate = 1, SWC = 0.43, SLA = 588, N = 4.64, Chl = 116.02, rawdata = readLi6800("/f5/2018-05-06-1002_f5-21-psd4-r1"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "PSD4", MT = 21, Replicate = 2, SWC = 0.42, SLA = 658, N = 4.42, Chl = NA, rawdata = readLi6800("/f5/2018-05-06-1228_f5-21-psd4-r2"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "PSD4", MT = 21, Replicate = 3, SWC = 0.41, SLA = 583, N = 3.78, Chl = NA, rawdata = readLi6800("/f5/2018-05-06-1429_f5-21-psd4-r3"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 552, N = 6.00, Chl = 205.36, rawdata = readLi6800("/f5/2018-05-07-0925_f5-21-c3-21"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "C", MT = 21, Replicate = 2, SWC = 0.8, SLA = 549, N = 5.81, Chl = 205.28, rawdata = readLi6800("/f5/2018-05-07-1222_f5-21-c3E-21"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "C", MT = 27, Replicate = 2, SWC = 0.8, SLA = 549, N = 6.00, Chl = NA, rawdata = readLi6800("/f5/2018-05-07-1436_F5-21-c3-27"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "PSD7", MT = 21, Replicate = 1, SWC = 0.25, SLA = 603, N = 3.76, Chl = 138.73, rawdata = readLi6800("/f5/2018-05-08-0921_f5-21-psd2"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "PSD7", MT = 21, Replicate = 2, SWC = 0.19, SLA = 480, N = 3.45, Chl = 89.48, rawdata = readLi6800("/f5/2018-05-08-1140_f5-21-psd2e"))
Li6800 %<>% add_row(Batch = "F5", GT = 21, Treatment = "Su7", MT = 21, Replicate = 1, SWC = 0.8, SLA = 551, N = 3.36, Chl = 139.73, rawdata = readLi6800("/f5/2018-05-08-1327_f5-21-Su"))
Li6800 %<>% add_row(Batch = "F5", GT = 27, Treatment = "C", MT = 27, Replicate = 1, SWC = 0.8, SLA = 656, N = 5.95, Chl = 151.03, rawdata = readLi6800("/f5/2018-05-09-0929_f5-27-ht-27"))
Li6800 %<>% add_row(Batch = "F5", GT = 27, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 656, N = 5.95, Chl = 151.03, rawdata = readLi6800("/f5/2018-05-09-1421_f5-27-ht-21"))
Li6800 %<>% add_row(Batch = "F5", GT = 27, Treatment = "D7", MT = 27, Replicate = 1, SWC = 0.11, SLA = 589, N = 5.46, Chl = 177.22, rawdata = readLi6800("/f5/2018-05-09-1242_f5-27-htd"))

# 19 files from F6
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "C", MT = 27, Replicate = 1, SWC = 0.8, SLA = 554, N = 6.10, Chl = NA, rawdata = readLi6800("/f6/2018-05-10-0919_f6-21-c2-27"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "C", MT = 21, Replicate = 1, SWC = 0.8, SLA = 554, N = 6.10, Chl = NA, rawdata = readLi6800("/f6/2018-05-10-1408_f6-21-c2-21"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "Su1", MT = 21, Replicate = 1, SWC = 0.8, SLA = 649, N = 5.15, Chl = 151.63, rawdata = readLi6800("/f6/2018-05-10-1226_f6-21-psd1"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "D4", MT = 21, Replicate = 1, SWC = 0.43, SLA = 587, N = 6.07, Chl = 185.25, rawdata = readLi6800("/f6/2018-05-11-0919_f6-21-d4"))
Li6800 %<>% add_row(Batch = "F6", GT = 27, Treatment = "D4", MT = 27, Replicate = 1, SWC = 0.41, SLA = 598, N = 5.05, Chl = 161.49, rawdata = readLi6800("/f6/2018-05-11-1216_f6-27-htd4"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "Su1", MT = 21, Replicate = 2, SWC = 0.8, SLA = 636, N = 6.04, Chl = 144.11, rawdata = readLi6800("/f6/2018-05-11-1449_f6-21-psd1"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "PSD4", MT = 21, Replicate = 1, SWC = 0.42, SLA = 623, N = 4.37, Chl = 141.41, rawdata = readLi6800("/f6/2018-05-13-0940_f6-21-psd4r1"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "PSD4", MT = 21, Replicate = 2, SWC = 0.44, SLA = 573, N = 4.29, Chl = 148.98, rawdata = readLi6800("/f6/2018-05-13-1258_f6-21-psd4r2"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "D7", MT = 21, Replicate = 1, SWC = 0.19, SLA = 528, N = 5.60, Chl = 223.60, rawdata = readLi6800("/f6/2018-05-14-0922_f6-21-d2"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "D7", MT = 21, Replicate = 2, SWC = 0.17, SLA = 486, N = 5.47, Chl = 208.51, rawdata = readLi6800("/f6/2018-05-14-1207_f6-21-d2e"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "Su7", MT = 21, Replicate = 1, SWC = 0.8, SLA = 499, N = 3.83, Chl = 141.18, rawdata = readLi6800("/f6/2018-05-15-0850_f6-21-su"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "PSD7", MT = 21, Replicate = 1, SWC = 0.22, SLA = 567, N = 3.87, Chl = 165.24, rawdata = readLi6800("/f6/2018-05-15-1059_f6-21-psd2"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "PSD7", MT = 21, Replicate = 2, SWC = 0.29, SLA = 484, N = 3.16, Chl = 183.69, rawdata = readLi6800("/f6/2018-05-15-1230_f6-21-psd2e"))
Li6800 %<>% add_row(Batch = "F6", GT = 27, Treatment = "D7", MT = 27, Replicate = 1, SWC = 0.16, SLA = 592, N = 4.78, Chl = 214.55, rawdata = readLi6800("/f6/2018-05-16-1148_f6-27-htd"))
Li6800 %<>% add_row(Batch = "F6", GT = 27, Treatment = "C", MT = 27, Replicate = 1, SWC = 0.8, SLA = 625, N = 5.62, Chl = 152.20, rawdata = readLi6800("/f6/2018-05-16-1332_f6-27-ht-27"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "C", MT = 21, Replicate = 2, SWC = 0.8, SLA = 539, N = 5.98, Chl = 104.47, rawdata = readLi6800("/f6/2018-05-17-0853_f6-21-c3e1-21"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "C", MT = 27, Replicate = 2, SWC = 0.8, SLA = 539, N = 5.98, Chl = 104.47, rawdata = readLi6800("/f6/2018-05-17-1247-f6-21-c3e1-27"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "C", MT = 21, Replicate = 3, SWC = 0.8, SLA = 532, N = 5.92, Chl = 212.27, rawdata = readLi6800("/f6/2018-05-18-0858_f6-21-c3e2-21"))
Li6800 %<>% add_row(Batch = "F6", GT = 21, Treatment = "C", MT = 27, Replicate = 3, SWC = 0.8, SLA = 532, N = 5.92, Chl = 212.27, rawdata = readLi6800("/f6/2018-05-18-1301_f6-21-c3e2-27"))

Li6800 = mutate(Li6800, SLA = SLA, # cm2/g
                Nmol = N/100/(SLA/1e4)/14) # mol N/m2



# Extract gas exchange traits for each plant -------------------------------------------

# # Fit FvCB model to low light and to part of the ACi curve (avoid TPU region)
# # Use Kmc and Kmo for each temperature Arabidopsis as calculated by Morales et al. (2018)
# # Parameters are calculate at their corresponding temperatures rather than normalized
# Hkmo = 2.91e4
# Hkmc = 4.94e4
# Kmo27 = 200*exp((27 + 273.15 - 298.15)*Hkmo/(298.15*8.31*(27 + 273.15)))
# Kmo21 = 200*exp((21 + 273.15 - 298.15)*Hkmo/(298.15*8.31*(21 + 273.15)))
# Kmc27 = 261*exp((27 + 273.15 - 298.15)*Hkmc/(298.15*8.31*(27 + 273.15)))
# Kmc21 = 261*exp((21 + 273.15 - 298.15)*Hkmc/(298.15*8.31*(21 + 273.15)))
# curves = vector("list", nrow(Li6800))
# for(i in 1:length(curves)) {
#   data = Li6800$rawdata[[i]]
#   HT = data$Tleaf[1] > 25
#   Kmc = ifelse(HT, Kmc27, Kmc21)
#   Kmo = ifelse(HT, Kmo27, Kmo21)
#   data = rename(data, PAR = Qin) %>% dplyr::select(PAR, Ci, A)
#   data = round(data, 2)
#   data = data[c(1:4,11:21),]
#   curvettypes = c(rep("PAR", 4), rep("Ci", nrow(data) - 7))
#   curves[[i]] = fitACiLRC(data = data, curves = curvettypes,
#                           pars = initial_ACiLRC(TPU = 100, KmC = Kmc, KmO = Kmo),
#                           priors = priors_ACiLRC(Vcmax = prior(45,20), 
#                                                  Jmax = prior(90,30)),
#                           fixed = c("KmC", "KmO","TPU","theta"))
# }

library(photofit)

fit_photosynthesis = function(data) {
  data = data[1:9,c("Q","A")]
  names(data) = c("PAR","A")
  fit = fitNRH(data, priors = priors_NRH(alpha = prior(0.05,0.02),
                                         theta = prior(0.7, 0.3)))
  coef = apply(fit$posterior, 2, median)
}

# Point estimates of photosynthesis parameters
Li6800 = mutate(Li6800,
                photopars = map(Li6800$rawdata, fit_photosynthesis),
                LUE = map_dbl(photopars, ~.x[["alpha"]]),
                Rd  = map_dbl(photopars, ~.x[["Rd"]]))


# Values calculated directly from the data
Li6800 = mutate(Li6800,
                A   = LUE*120 - Rd,
                gsw = map_dbl(rawdata, ~mean(.x$gsw[3:4])),
                WUE = map_dbl(rawdata, ~mean(.x$A[3:4]/.x$E[3:4])/1000)) # mmol/mol

TreatmentCode = c(C_21 = "C", Su7_21 = "Su7", PSD7_21 = "PSD7",
                  D7_21 = "D7", C_27 = "HT", D7_27 = "HTD27",
                  Su1_21 = "Su1", D4_21 = "D4",
                  D4_27 = "HTD4", PSD4_21 = "PSD4")

Li6800 = mutate(Li6800,
                FullTreatment = TreatmentCode[paste(Treatment, GT, sep = "_")])

GasExchange = select(Li6800, Batch, GT, Replicate, MT, FullTreatment, N, Nmol, 
                     Chl, SWC, SLA, gsw, WUE, Rd, LUE, A)

# Save processed dataset --------------------------------------------------

saveRDS(GasExchange, file = "Intermediate/Photosynthesis/AllTraits.rds")


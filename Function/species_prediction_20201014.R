#Environment setup
library(ggplot2)

calculate_HL_one_sample_one_direction <- function(abund, sub_path, direction){
  #Find direction data
  sub_path <- paste(sub_path, 'DD_S-map', direction, sep = '/')
  analyzed_states <- dir(sub_path)
  analyzed_states <- analyzed_states[-1*(length(analyzed_states)-1):length(analyzed_states)]
  
  #Collect coefs
  coefs_loc <- paste(sub_path, 'coefs', sep = '/')
  coef_files <- dir(coefs_loc)
  
  #Calculate harmony level
  HL <- array(0, dim = c(nrow(abund), ncol(abund)))
  for (i in 1:length(coef_files)) {
    #Read analyzed state file
    analyzed_states_loc <- read.csv(paste(sub_path, analyzed_states[i], sep = '/'), 
                                    row.names = 1)
    
    #Process coefs files
    coefs <- read.csv(paste(coefs_loc, coef_files[i], sep = '/'), row.names = 1)
    interested_otu <- as.numeric(strsplit(coef_files[i], split='_')[[1]][1])
    target_otu <- as.numeric(strsplit(coef_files[i], split='_')[[1]][2])
    for (j in 1:(nrow(coefs)-1)) {
      if (abund[j, target_otu+1] > abund_threshold) {
        HL[analyzed_states_loc[j,1]+1, interested_otu+1] <- HL[analyzed_states_loc[j,1]+1, 
                                                               interested_otu+1] + sum(abund[analyzed_states_loc[j,1]+1, target_otu+1]*coefs[j,])
      }
    }
  }
  
  return(HL)
}

calculate_VCR_one_sample_one_direction <- function(abund, sub_path, direction){
  #Find direction data
  sub_path <- paste(sub_path, 'DD_S-map', direction, sep = '/')
  analyzed_states <- dir(sub_path)
  analyzed_states <- analyzed_states[-1*(length(analyzed_states)-1):length(analyzed_states)]
  
  #Collect coefs
  coefs_loc <- paste(sub_path, 'coefs', sep = '/')
  coef_files <- dir(coefs_loc)
  
  #Calculate harmony level
  VCR <- array(0, dim = c(nrow(abund), ncol(abund)))
  for (i in 1:length(coef_files)) {
    #Read analyzed state file
    analyzed_states_loc <- read.csv(paste(sub_path, analyzed_states[i], sep = '/'), 
                                    row.names = 1)
    
    #Process coefs files
    coefs <- read.csv(paste(coefs_loc, coef_files[i], sep = '/'), row.names = 1)
    interested_otu <- as.numeric(strsplit(coef_files[i], split='_')[[1]][1])
    target_otu <- as.numeric(strsplit(coef_files[i], split='_')[[1]][2])
    for (j in 1:nrow(coefs)) {
      if (abund[j, target_otu+1] > abund_threshold) {
        VCR[analyzed_states_loc[j,1]+1, interested_otu+1] <- VCR[analyzed_states_loc[j,1]+1, interested_otu+1] + coefs[j,target_otu+1]
      }
    }
  }
  
  return(VCR)
}

#Preparation
info <- array(0, dim = c(1,10))
abund_threshold = 1
change_threshold = .5
change_fold = 1

results <- array(NA, dim = c(19, 7))
colnames(results) <- c('Sample','Sampling_interval', 'HL', 'VCR', 'SI', 'AI', 'HR')
results[,1] <- c('55.39953_10.41479', '55.64929_12.05857', '56.21314_10.24247',
                 '56.56513_9.04216', '57.04516_10.04576', '57.04951_9.86479', '57.42126_9.97541',
                 'HK_ST', 'HK_SWHA', 'JX', 'yearXH_sample',
                 'HK_Daily', 'HK_M1', 'HK_M2', 'HK_M3', 'HK_M4', 'HK_M5', 'HK_M6', 'month_sample'
)
results[,2] <- c('Quarterly', 'Quarterly', 'Quarterly', 
                 'Quarterly', 'Quarterly', 'Quarterly', 'Quarterly', 
                 'Biweekly', 'Biweekly', 'Biweekly', 'Biweekly', 
                 'Daily', 'Daily', 'Daily', 'Daily', 'Daily', 'Daily', 'Daily', 'Daily')




#Find data folder
path <- 'E:/学习/研究生/研究生课题/大规模污泥群落数据统计/数据/outputs/taxa/5'
samples <- dir(path)

#Collect info for each sample
for (base in c('HL', 'VCR')) { #, 'SI', 'AI', 'HR'
  print(base)
  
  for (sample in samples) {
    print(sample)
    #Input abund file
    sub_path <- paste(path, sample, sep = '/')
    abund <- read.csv(paste(sub_path, '/abund.csv', sep = ''), row.names = 1)
    
    #Real change
    real_change <- array(0, dim = c(nrow(abund)-1, ncol(abund)))
    for (i in 1:nrow(real_change)) {
      for (j in 1:ncol(real_change)) {
        ###Show best prediction when abund_threshold and change_threshold = 2
        # if (abund[i+1,j] > abund[i,j]*change_threshold) {
        #   real_change[i,j] <- 1
        # }
        # if (abund[i+1,j] < abund[i,j]/change_threshold) {
        #   real_change[i,j] <- -1
        # }
        
        if (abund[i+1,j] > abund[i,j]+change_threshold) {
          real_change[i,j] <- 1
        }
        if (abund[i+1,j] < abund[i,j]-change_threshold) {
          real_change[i,j] <- -1
        }
      }
    }
    
    #Predict change direction
    predicted_change <- array(1000, dim = c(nrow(abund)-1, ncol(abund)))
    
    ##According to absolute HL difference
    if (base == 'HL') {
      #Calculate the harmony level at increase and decrease direction
      HL_i <- calculate_HL_one_sample_one_direction(abund, sub_path, 'increase')
      HL_d <- calculate_HL_one_sample_one_direction(abund, sub_path, 'decrease')
      
      for (i in 1:nrow(real_change)) {
        for (j in 1:ncol(real_change)) {
          if (abs(HL_i[i,j]) > abs(HL_d[i,j])*change_fold) {
            predicted_change[i,j] <- 1
          }
          if (abs(HL_d[i,j]) > abs(HL_i[i,j])*change_fold) {
            predicted_change[i,j] <- -1
          }
          if (abs(HL_d[i,j]) == abs(HL_i[i,j])*change_fold) {
            predicted_change[i,j] <- 0
          }
        }
      }
    }
    
    ##According to absolute VCR difference
    if (base == 'VCR') {
      #Calculate the VCR at increase and decrease direction
      VCR_i <- calculate_VCR_one_sample_one_direction(abund, sub_path, 'increase')
      VCR_d <- calculate_VCR_one_sample_one_direction(abund, sub_path, 'decrease')
      
      for (i in 1:nrow(real_change)) {
        for (j in 1:ncol(real_change)) {
          if (abs(VCR_i[i,j]) > abs(VCR_d[i,j])*change_fold) {
            predicted_change[i,j] <- 1
          }
          if (abs(VCR_d[i,j]) > abs(VCR_i[i,j])*change_fold) {
            predicted_change[i,j] <- -1
          }
          if (abs(VCR_d[i,j]) == abs(VCR_i[i,j])*change_fold) {
            predicted_change[i,j] <- 0
          }
        }
      }
    }
    
    
    ##Accuracy summarization
    accuracy_matrix <- array(0, dim = c(nrow(abund)-1, ncol(abund)))
    for (i in 1:nrow(accuracy_matrix)) {
      for (j in 1:ncol(accuracy_matrix)) {
        if (predicted_change[i,j] != real_change[i,j]) {
          accuracy_matrix[i,j] = 1
        }
      }
    }
    ##More accurate with higher change threshold and change fold
    accuracy <- length(which(accuracy_matrix == 1))/length(accuracy_matrix)*100
    print(accuracy)
    
    results[which(results[,1] == sample), which(colnames(results) == base)] <- accuracy
  }
  
}

#Plot
results <- as.data.frame(results)
for (i in 3:7) {
  results[,i] <- as.numeric(as.character(results[,i]))
}

for (i in 1:2) {
  p <- ggplot(results, aes(x=Sampling_interval, y=results[,i+2]))+
    theme(panel.grid.major =element_line(size = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.line = element_line(linetype = "solid", size = 0.5),
          panel.background = element_rect(color = "black", size = 0.5),
          panel.grid.minor = element_blank(),
          legend.position = 'none'
    )+
    geom_boxplot()+
    ylim(40, 100)+
    geom_hline(yintercept = 66, size=1.5, color='red', linetype='dashed')+
    scale_x_discrete(limits=c('Daily', 'Biweekly', 'Quarterly'))
  
  print(p)
}


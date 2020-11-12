#Objective
##Try to find the abundance-dependent pattern of HL and HR

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

#Parameter setting
abund_threshold = 1
change_threshold = .5
change_fold = 1

#Data location input
path <- 'E:/学习/研究生/研究生课题/大规模污泥群落数据统计/数据/outputs/manifold'

##Input abundance file
abund_file <- paste(path, 'abund.csv', sep = '/')
abund <- read.csv(abund_file, row.names = 1)

#Calculate all kinds of parameter
##Harmony level
HL_i <- calculate_HL_one_sample_one_direction(abund, path, 'increase')
HL_d <- calculate_HL_one_sample_one_direction(abund, path, 'decrease')

##Structural stability
VCR_i <- calculate_VCR_one_sample_one_direction(abund, path, 'increase')
VCR_d <- calculate_VCR_one_sample_one_direction(abund, path, 'decrease')

#Make input
for (j in 1:11) {
  input <- cbind(abund[,j], 0)
  # for (i in 1:nrow(abund)) {
  #   if (abs(HL_i[i,j]) > abs(HL_d[i,j])*change_fold) {
  #     input[i,2] <- 1
  #   }
  #   if (abs(HL_d[i,j]) > abs(HL_i[i,j])*change_fold) {
  #     input[i,2] <- -1
  #   }
  # }
  for (i in 1:nrow(abund)) {
    if (abs(VCR_i[i,j]) > abs(VCR_d[i,j])*change_fold) {
      input[i,2] <- 1
    }
    if (abs(VCR_d[i,j]) > abs(VCR_i[i,j])*change_fold) {
      input[i,2] <- -1
    }
  }
  
  colnames(input) <- c('Abundance', 'Change_Direction')
  input <- as.data.frame(input)
  input[,1] <- as.numeric(as.character(input[,1]))
  input[,2] <- as.character(input[,2])
  taxonomy <- strsplit(colnames(abund)[j], split = '__')[[1]]
  title = paste(sub('.D_4', '', taxonomy[5]), sub('.D_5', '', taxonomy[6]), sep = ' ')
  
  #Plot
  p <- 
  ggplot(input,aes(x=Abundance, fill=Change_Direction, color=Change_Direction))+
    theme(panel.grid.major =element_line(size = 1),
          axis.title = element_text(size = 16, colour = "black"),
          axis.text = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.line = element_line(linetype = "solid", size = 0.5),
          panel.background = element_rect(color = "black", size = 0.5),
          panel.grid.minor = element_blank(),
          legend.position = 'none'
    )+
    geom_density(alpha=.3, size=1,
                 position="stack")+
    # scale_fill_manual(values=c('#C52427','#87AED6'))+
    # scale_color_manual(values = c('#C52427','#87AED6'))+
    xlab('Abundances (%)')+
    ylab('')
  
  print(p)
  
  p <- 
  ggplot(input, aes(x=Abundance))+
    theme(panel.grid.major =element_line(size = 1),
          axis.title = element_text(size = 16, colour = "black"),
          axis.text = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.line = element_line(linetype = "solid", size = 0.5),
          panel.background = element_rect(color = "black", size = 0.5),
          panel.grid.minor = element_blank()
    )+
    geom_density(fill='lightblue', color='darkblue', size=1)+
    labs(title = title)+
    xlab('Abundances (%)')+
    ylab('')
  
  print(p)
}


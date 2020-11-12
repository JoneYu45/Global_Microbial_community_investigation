#Environment setup
library(ggplot2)

#Preparations
thts <- c(0.1, 1, 10, 2)
info <- array(0, dim = c(1,8))

#Data input
path <- 'D:/PythonProject/Density_dependency_S_map/DD_S-map_2/outputs/tht_wo_T1'
tht_samples <- dir(path)

for (j in 1:length(tht_samples)) {
  #Find WWTPs
  tht_sample <- tht_samples[j]
  sub_path <- paste(path, tht_sample, sep = '/')
  WWTPs <- dir(sub_path)
  
  #Collect data from each WWTP
  for (WWTP in WWTPs) {
    #Find fit results from both directions
    for (direction in c('increase','decrease')) {
      sub_path_2 <- paste(sub_path, WWTP, 'DD_S-map', direction, 'fit_result', sep = '/')
      fit_result_files <- dir(sub_path_2)
      
      #Collect info in sequence
      for (i in 1:length(fit_result_files)) {
        file_info <- unlist(strsplit(fit_result_files[i], split = '_'))
        workbook <- paste(sub_path_2, fit_result_files[i], sep = '/')
        fit_result <- read.csv(workbook)
        info <- rbind(info, cbind(file_info[1], thts[j], fit_result$X, 
                                  fit_result$RMSE_o, fit_result$RMSE_o/fit_result$Std_o, 
                                  fit_result$Test.set.score, fit_result$ymax_test, fit_result$ymin_test))
      }
    }
  }
}

##Make input
info <- info[-1,]
colnames(info) <- c('target_otu', 'Theta', 'State', 'RMSE', 
                    'RMSE/STD', 'Test_score', 'max', 'min')
info <- as.data.frame(info)
for (i in 1:ncol(info)) {
  info[,i] <- as.numeric(as.character(info[,i]))
}
info$Theta <- as.character(info$Theta)

##Plot
ggplot(info, aes(x=Theta, y=`RMSE/STD`))+
  theme(panel.grid.major =element_line(size = 1),
        axis.title = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.line = element_line(linetype = "solid", size = 0.5),
        panel.background = element_rect(color = "black", size = 0.5),
        panel.grid.minor = element_blank()
  )+
  geom_boxplot(fill=c('#C52427', '#87AED6', '#7BAC53')[1])+
  scale_x_discrete(limits=c('0.1', '1', '2', '10'))+
  ylim(0,2)+
  geom_hline(yintercept=1, linetype='dashed', color='red', size=1.5)
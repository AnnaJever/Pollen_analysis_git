### Analysis of manuscript with title:

# Palaeoecological data reveal long-term trajectories of non-native vegetation on islands globally
# Walentowitz et al.
# last edited on 16. November 2022



### packages
library(readxl)
library(tidyverse)
library(zoo) # moving average
library(lme4) # lmms
library(car) # anova
library(report)
library(redres) # residuals
library(segmented) # piecewise regression



### functions
star <- function(x){a<-""
  if(x>0.05){a<-"n.s."}
  if(x<=0.05){a<-"*"}	
  if(x<=0.01){a<-"**"}
  if(x<=0.001){a<-"***"}	
  a}



### get data
load("./Data/IslandPollenExtended.RData") # pollen data



### prepare table for results of wilcoxon test
box_table_lower <- as.data.frame(matrix(nrow = 29, ncol = 6))
colnames(box_table_lower) <- c("island", "pvalue", "mean_b", "mean_a", "human.arrival", "end_value_ali")
rownames(box_table_lower) <- env$island

box_table_upper <- as.data.frame(matrix(nrow = 29, ncol = 6))
colnames(box_table_upper) <- c("island", "pvalue", "mean_b", "mean_a", "human.arrival", "end_value_ali")
rownames(box_table_upper) <- env$island



### data preparation and single island analsis

for(i in 1:nrow(env)){
  if(!is.na(env$pollen_assign[i]) &&
     # these six islands were excluded due to data resolution etc. 
     !(env$island[i] %in% c("Robinson Crusoe",
                            "Nightingale Island",
                            "Iceland",
                            "New Caledonia",
                            "La Gomera",
                            "Hispaniola"))){
 
    # get data
    pollen <- as.data.frame(specs[[paste0(env$ID[i])]])
 
    # categorization
    nat           <- pollenstatus$pollen_taxa[which(pollenstatus$assigned_dataset == env$pollen_assign[i]
                                        & pollenstatus$status == "native")]

    ali           <- pollenstatus$pollen_taxa[which(pollenstatus$assigned_dataset == env$pollen_assign[i]
                                        & pollenstatus$status == "alien")]

    extended      <- pollenstatus$pollen_taxa[which(pollenstatus$assigned_dataset == env$pollen_assign[i]
                                        & pollenstatus$status %in% c("native", "mixed"))]

    crop          <- pollenstatus$pollen_taxa[which(pollenstatus$assigned_dataset == env$pollen_assign[i]
                                                    & pollenstatus$crop == "yes")]
    
    endemic       <- pollenstatus$pollen_taxa[which(pollenstatus$assigned_dataset == env$pollen_assign[i]
                                                    & pollenstatus$endemic == "yes")]
    
    
    # all taxa
    tax <- unique(colnames(pollen))
    tax_filtered <- unique(colnames(pollen[which(colnames(pollen) %in% c(extended, ali))]))
    
    # Re-scaling
    # final data ranges between 0 and 100
    
    # dataset_lower
    pollen_lower <- pollen[which(colnames(pollen) %in%
                             c(extended, ali))]
    pollen_lower <- data.matrix(sapply(pollen_lower, as.numeric))
    pollen_lower <- as.data.frame(pollen_lower)

    for(j in 1:nrow(pollen_lower)){
      pollen_lower[j,] <- (pollen_lower[j,]/rowSums(pollen_lower)[j])*100
    }
    
    # dataset_upper
    pollen_upper <- pollen[which(colnames(pollen) %in%
                              c(nat, ali))]
    pollen_upper <- data.matrix(sapply(pollen_upper, as.numeric))
    pollen_upper <- as.data.frame(pollen_upper)
    
    for(j in 1:nrow(pollen_upper)){
      pollen_upper[j,]<- (pollen_upper[j,]/rowSums(pollen_upper)[j])*100
    }
    
    # dataset raw
    pollen_raw <- pollen
    pollen_raw <- data.matrix(sapply(pollen_raw, as.numeric))
    pollen_raw <- as.data.frame(pollen_raw)
    
    for(j in 1:nrow(pollen_raw)){
      pollen_raw[j,]<- pollen_raw[j,]/rowSums(pollen_raw)[j]*100
    }
    
    # Pre-processing
    # build total sums of native and non-native pollen
    pollen_df <- as.data.frame(as.numeric(ages[[paste0(env$ID[i])]]))
    colnames(pollen_df) <- c("time")
    
    pollen_df$sum.inv.lower <- round(rowSums(pollen_lower[which(colnames(pollen_lower) %in% ali)]))
    pollen_df$sum.inv.upper <- round(rowSums(pollen_upper[which(colnames(pollen_upper) %in% ali)]))
    pollen_df$sum.endem     <- round(rowSums(pollen_raw[which(colnames(pollen_raw) %in% endemic)]))
    pollen_df$sum.crops     <- round(rowSums(pollen_raw[which(colnames(pollen_raw) %in% crop)]))
    

    # number of pollen taxa
    for(k in 1:nrow(pollen_df)){
 
      temp.ali <- pollen_lower[k, which(colnames(pollen_lower) %in% ali)]
      pollen_df$sp.ali.lower[k] <- length(temp.ali[which(temp.ali > 0)])
      
      temp.ext <- pollen_upper[k, which(colnames(pollen_upper) %in% ali)]
      pollen_df$sp.ali.upper[k] <- length(temp.ext[which(temp.ext > 0)])

    }
    
    # final dataset cut off data at 5000 years before now
    pollen_df <- pollen_df[order(pollen_df$time, decreasing = T),]
    pollen_df[,c(2:5)] <- round(pollen_df[,c(2:5)])
    pollen_df <- pollen_df[complete.cases(pollen_df),]
    pollen_df <- pollen_df[pollen_df$time <= 5000,]
    
    
  
    ### Single island figure
    # Figures with moving averages
    # Take the widow size so that ~ 25 averages are being made.
    # This means that each average covers about 200 years.
    # Additional condition: at least 2 data points must be averaged
    # k is thus individual for each island
    
    if(nrow(pollen_df) >= 25){
    k <- round(nrow(pollen_df)/25)
      if(k == 1){
      k <- 2
      }
    }else{
    k <- 2
    }
    
    pollen_df <- pollen_df %>%
      mutate(amb.av_mw.lower = rollmean(sum.inv.lower, k = k, fill = NA))

    pollen_df <- pollen_df %>%
      mutate(amb.av_mw.upper = rollmean(sum.inv.upper, k = k, fill = NA))
    
    # build categories (before and after human settlement)
    pollen_df$time <- round(pollen_df$time)
    pollen_df$category <- NA
    
    pollen_df$category[
      which(pollen_df$time
            %in% c(5000:env$human.arrival[i]))] <- "b" # before human arrival
    
    pollen_df$category[
      which(pollen_df$time <= env$human.arrival[i])] <- "a" # after human arrival
    
    pollen_df$category <- as.factor(pollen_df$category)
    pollen_df$category <- factor(pollen_df$category,
                                 levels = c("b", "a"))
    
    pollen_df_b <- pollen_df %>% 
      filter(category == "b")
    
    pollen_df_a <- pollen_df %>% 
      filter(category == "a")
    
    # Additional line for moving average
    # vertial black line indicating human arrival
    plot_mw_lower <- ggplot(pollen_df, aes(time, sum.inv.lower)) +
      geom_point(aes(time,sum.inv.lower), col = "white")+
      geom_point(data = pollen_df_a, aes(time,sum.inv.lower), col = "salmon") + 
      geom_point(data = pollen_df_b, aes(time,sum.inv.lower), col = "darkgrey")+
      geom_line(aes(time, amb.av_mw.lower),color="coral4") + 
      geom_line(aes(time, amb.av_mw.lower),color="gray30") +
      theme(panel.background = element_rect(fill = "NA", colour = "black")) +
      xlim(max(pollen_df$time),min(pollen_df$time))+
      geom_vline(xintercept=as.numeric(env$human.arrival[
        which(env$island == env$island[i])])) +
      ggtitle(env$island[i])+
      theme(axis.text = element_text(size = 15))+
      theme(plot.title = element_text(size = 15))+
      theme(plot.title = element_text(face = "bold"))+
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank())

    svg(file=paste0("Figures/trend/lower_trend_", env$island[i], ".svg"),
        height = 2, width = 1.7)
    print(plot_mw_lower)
    dev.off()
    
    
    
    plot_mw_upper <- ggplot(pollen_df, aes(time, sum.inv.upper)) +
      geom_point(aes(time,sum.inv.upper), col = "white")+
      geom_point(data = pollen_df_a, aes(time,sum.inv.upper), col = "salmon") + 
      geom_point(data = pollen_df_b, aes(time,sum.inv.upper), col = "darkgrey")+
      geom_line(aes(time, amb.av_mw.upper),color="coral4") + 
      geom_line(aes(time, amb.av_mw.upper),color="gray30") + 
      theme(panel.background = element_rect(fill = "NA", colour = "black")) +
      xlim(max(pollen_df$time),min(pollen_df$time))+
      geom_vline(xintercept=as.numeric(env$human.arrival[
        which(env$island == env$island[i])])) +
      ggtitle(env$island[i])+
      theme(axis.text = element_text(size = 15))+
      theme(plot.title = element_text(size = 15))+
      theme(plot.title = element_text(face = "bold"))+
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank())
    
    svg(file=paste0("Figures/trend/upper_trend_", env$island[i], ".svg"),
        height = 2, width = 1.7)
    print(plot_mw_upper)
    dev.off() 


    
    # boxplots: sum.inv before and after human arrival (ha)
    # Vanuatu is excluded here as ha happened before pollen series starts
    
    if(env$island[i] != "Vanuatu"){
      
      # two group test
      wil_lower <- wilcox.test(pollen_df$sum.inv.lower ~ pollen_df$category)
      wil_upper <- wilcox.test(pollen_df$sum.inv.upper ~ pollen_df$category)
      
    if(!is.na(wil_lower$p.value)){
        
    boxplot_lower <- ggplot(aes(x=category, y=sum.inv.lower, fill = category), data = pollen_df) +
      geom_boxplot() +
      scale_fill_manual(values=c("darkgrey", "salmon"))+
      geom_jitter(color="black", size=0.4, alpha=0.9) +
      theme(panel.background = element_rect(fill = "NA", colour = "NA")) +
      theme(legend.position="none")+
      annotate("text", x = 2, y = max(pollen_df$sum.inv.lower),
               label= paste0(star(wil_lower$p.value)))+
      ggtitle(env$island[i])+
      theme(axis.text = element_text(size = 15, colour = "white"))+
      theme(plot.title = element_text(size = 15, colour = "white"))+
      theme(plot.title = element_text(face = "bold"))+
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank())
    
      svg(file=paste0("Figures/trend/lower_boxp_", env$island[i], ".svg"),
          height = 2, width = 0.7)
      print(boxplot_lower)
      dev.off()
      

      box_table_lower$island[i] <- env$island[i]
      box_table_lower$pvalue[i] <- wil_lower$p.value
      box_table_lower$mean_b[i] <- mean(pollen_df_b$sum.inv.lower)
      box_table_lower$mean_a[i] <- mean(pollen_df_a$sum.inv.lower)
      box_table_lower$human.arrival[i]   <- env$human.arrival[i]
      box_table_lower$end_value_ali[i] <- pollen_df$sum.inv.lower[which(pollen_df$time == min(pollen_df$time))]
    }
      
      
    if(!is.na(wil_upper$p.value)){
      
      boxplot_upper <- ggplot(aes(x=category, y=sum.inv.upper, fill = category), data = pollen_df) +
        geom_boxplot() +
        scale_fill_manual(values=c("darkgrey", "salmon"))+
        geom_jitter(color="black", size=0.4, alpha=0.9) +
        theme(panel.background = element_rect(fill = "NA", colour = "NA")) +
        theme(legend.position="none")+
        annotate("text", x = 2, y = max(pollen_df$sum.inv.upper),
                 label= paste0(star(wil_upper$p.value)))+
        ggtitle(env$island[i])+
        theme(axis.text = element_text(size = 15, colour = "white"))+
        theme(plot.title = element_text(size = 15, colour = "white"))+
        theme(plot.title = element_text(face = "bold"))+
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank())
      
      svg(file=paste0("Figures/trend/upper_boxp_", env$island[i], ".svg"),
          height = 2, width = 0.7)
      print(boxplot_upper)
      dev.off()    
      
      
    box_table_upper$island[i] <- env$island[i]
    box_table_upper$pvalue[i] <- wil_upper$p.value
    box_table_upper$mean_b[i] <- mean(pollen_df_b$sum.inv.upper)
    box_table_upper$mean_a[i] <- mean(pollen_df_a$sum.inv.upper)
    box_table_upper$human.arrival[i]   <- env$human.arrival[i]
    box_table_upper$end_value_ali[i] <- pollen_df$sum.inv.upper[which(pollen_df$time == min(pollen_df$time))]
  }
    }
    
    
    
    #### prepare multi-island data
    

    if(i == 1){
      pollen_df$island <- env$island[i]
      pollen_all <- pollen_df
      pollen_tax <- as.data.frame(tax)
      pollen_tax$island <- env$ID[i]
      pollen_tax_filtered <- as.data.frame(tax_filtered)
      pollen_tax_filtered$island <- env$ID[i]
      pollen_all_raw <- pollen_raw 
    }else{
      pollen_df$island <- env$island[i]
      pollen_all <- rbind(pollen_all, pollen_df)
      temp_tax <- as.data.frame(tax)
      temp_tax$island  <- env$ID[i]
      pollen_tax <- rbind(pollen_tax, temp_tax)
      temp_tax_filtered <- as.data.frame(tax_filtered)
      temp_tax_filtered$island  <- env$ID[i]
      pollen_tax_filtered <- rbind(pollen_tax_filtered, temp_tax_filtered)
    }
    
      print(paste(env$island[i], "is done."))
    
  }
}

# cleaning
rm(temp.ali,temp.all,i, k,j,temp.exe, wil_lower, wil_upper, tax_filtered)


# descriptive statistics single island analysis
# lower
box_table_lower <- box_table_lower[which(complete.cases(box_table_lower)),]
box_table_lower$change <- round(box_table_lower$mean_a - box_table_lower$mean_b, digits = 2)
mean(box_table_lower$mean_a) #5%
mean(box_table_lower$mean_b) #1%
mean(box_table_lower$end_value_ali) #8%

# upper
box_table_upper <- box_table_upper[which(complete.cases(box_table_upper)),]
box_table_upper$change <- round(box_table_upper$mean_a - box_table_upper$mean_b, digits = 2)
mean(box_table_upper$mean_a) #16%
mean(box_table_upper$mean_b) #6%
mean(box_table_upper$end_value_ali) #25%



# multi-island analysis

# LMM
# build LMM to explain sum of non-native pollen with time and with islands as random effect

# remove Foa and Ha'afeva from multi-island analysis. The pollen data had a very low temporal distribution

pollen_all <- pollen_all[which(pollen_all$island != c("Ha'afeva", "Foa")),]

# data preparation for LMM
hist(pollen_all$sum.inv.lower, freq = FALSE)
hist(pollen_all$sum.inv.upper, freq = FALSE)
pollen_all$island <- as.factor(pollen_all$island)

# model non-natives
all_pollen_glmm1a <- lmer(log1p(sum.inv.lower) ~ time + (1|island), data = pollen_all)
all_pollen_glmm1b <- lmer(log1p(sum.inv.upper) ~ time + (1|island), data = pollen_all)

  print(all_pollen_glmm1a, corr=F)
  summary(all_pollen_glmm1a)
  Anova(all_pollen_glmm1a)
  ranef(all_pollen_glmm1a)
  rc_resids1a <- compute_redres(all_pollen_glmm1a)
  hist(rc_resids1a)
  table.nonnata <- as.data.frame(report_table(all_pollen_glmm1a))
  write.csv2(table.nonnata, "table.nonnat_lower.csv")

  print(all_pollen_glmm1b, corr=F)
  summary(all_pollen_glmm1b)
  Anova(all_pollen_glmm1b)
  ranef(all_pollen_glmm1b)
  rc_resids1b <- compute_redres(all_pollen_glmm1b)
  hist(rc_resids1b)
  table.nonnatb <- as.data.frame(report_table(all_pollen_glmm1b))
  write.csv2(table.nonnatb, "table.nonnat_upper.csv")



### piecewise regression
# reduced to last 2000 years
  
pollen_2000 <- pollen_all[which(pollen_all$time <= 2000),]

#fit simple linear regression model

fit1 <- lm(sum.inv.lower ~ time, data=pollen_2000)
fit2 <- lm(sum.inv.upper ~ time, data=pollen_2000)

# fit piecewise regression model to original model, estimating a breakpoint
segmented.fit1 <- segmented(fit1, seg.Z = ~time, control = seg.control(display = FALSE))
dat1 = data.frame(x = pollen_2000$time, y = broken.line(segmented.fit1)$fit)

segmented.fit2 <- segmented(fit2, seg.Z = ~time, control = seg.control(display = FALSE))
dat2 = data.frame(x = pollen_2000$time, y = broken.line(segmented.fit2)$fit)

#view summary of segmented model
summary(segmented.fit1) # 575 
plot(segmented.fit1)

summary(segmented.fit2) # 280 
plot(segmented.fit2)

table.nonnat_piecwise1 <- as.data.frame(report_table(segmented.fit1))
write.csv2(table.nonnat_piecwise1, "table.nonnat_lower.csv")

table.nonnat_piecwise2 <- as.data.frame(report_table(segmented.fit2))
write.csv2(table.nonnat_piecwise2, "table.nonnat_upper.csv")

mean(pollen_all$sum.inv.lower[which(pollen_all$time < 575)])# 5%
mean(pollen_all$sum.inv.lower[which(pollen_all$time > 575)])# 2%

mean(pollen_all$sum.inv.upper[which(pollen_all$time < 280)])# 20%
mean(pollen_all$sum.inv.upper[which(pollen_all$time > 280)])# 7%



# correlation pollen percentages and pollen taxa
cor(pollen_all$sum.inv.lower, pollen_all$sp.ali.lower) # 0.5
cor.test(pollen_all$sum.inv.lower, pollen_all$sp.ali.lower)

cor(pollen_all$sum.inv.upper, pollen_all$sp.ali.upper) # 0.4
cor.test(pollen_all$sum.inv.upper, pollen_all$sp.ali.upper)



# lmm endemics
all_pollen_glmm2 <- lmer(log1p(sum.endem) ~ time + (1|island), data = pollen_all)

print(all_pollen_glmm2, corr=F)
summary(all_pollen_glmm2)
Anova(all_pollen_glmm2)

rc_resids2 <- compute_redres(all_pollen_glmm2)
hist(rc_resids2)

table.endemics <- as.data.frame(report_table(all_pollen_glmm2)) 
write.csv2(table.endemics, "table.endemics.csv")

# model cultivars
all_pollen_glmm3 <- lmer(log1p(sum.crops) ~ time + (1|island), data = pollen_all)

print(all_pollen_glmm3, corr=F)
summary(all_pollen_glmm3)
Anova(all_pollen_glmm3)

rc_resids3 <- compute_redres(all_pollen_glmm3)
hist(rc_resids3)

table.cultivars <- as.data.frame(report_table(all_pollen_glmm3)) 
write.csv2(table.cultivars, "table.cultivars.csv")



### multi-island plots
# plot pollen abundance 5000 years (affinity 300*300)
plot_abu <- ggplot(pollen_all) +
  geom_point(aes(time,sum.inv.lower), col = "indianred4", alpha = 0.3, size = 0.5) +
  geom_point(aes(time,sum.inv.upper), col = "salmon", alpha = 0.3, size = 0.5) +
  geom_smooth(aes(time,sum.inv.upper), span = 1,
              se = T, colour = "salmon")+
  geom_smooth(aes(time,sum.inv.lower), span = 1,
              se = T, colour = "indianred4")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  labs(x = "Time [Cal. years BP]", y = "Pollen [%]")+
  ylim(0,50)+
  xlim(5000,-70)

plot_abu


# breakpoint (affinity 300*300)
plot_bp <- ggplot(pollen_all) +
  geom_point(aes(time,sum.inv.lower), col = "indianred4", alpha = 0.3, size = 0.5) +
  geom_point(aes(time,sum.inv.upper), col = "salmon", alpha = 0.3, size = 0.5) +
  geom_line(aes(x = x, y = y), data = dat2, color = 'salmon', size = 2)+
  geom_line(aes(x = x, y = y), data = dat1, color = 'indianred4', size = 2)+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Time [Cal. years BP]", y = "Non-native pollen [%]")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  xlim(2000,-70)+
  ylim(0,50)

plot_bp



# corr species and abundance (affinity 300*300)
plot_cor.lower <- ggplot(pollen_all) +
  geom_point(aes(sum.inv.upper,sp.ali.lower), col = "salmon")+
  geom_smooth(aes(sum.inv.upper,sp.ali.lower), method = "lm", col = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Pollen [%]", y = "No. taxa")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))

plot_cor.lower

plot_cor.upper <- ggplot(pollen_all) +
  geom_point(aes(sum.inv.lower,sp.ali.lower), col = "salmon")+
  geom_smooth(aes(sum.inv.lower,sp.ali.lower), method = "lm", col = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Pollen [%]", y = "No. taxa")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))

plot_cor.upper




# crop pollen in time
pollen_all_crop <- pollen_all[which(pollen_all$sum.crops > 0 & !is.na(pollen_all$sum.crops)),]

# build GLMM explain sum of crop pollen with time and with islands as random effect
# extreme values in Gran Canaria. Add the following term to exclude
# pollen_all_crop$island != "Gran Canaria"
# plot pollen abundance 5000 years (affinity 300*300)

plot_crop <- ggplot(pollen_all_crop[which(pollen_all_crop$time <= 2000),])+
  geom_point(aes(time,sum.crops), col = "salmon", alpha = 0.5) +
  geom_smooth(aes(time,sum.crops), span = 1,
              se = T, colour = "salmon")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Time [Cal. years BP]", y = "Cultivar pollen [%]")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  xlim(2000,0)

plot_crop

# endemic pollen in time
pollen_all_endem <- pollen_all[which(pollen_all$sum.endem > 0 & !is.na(pollen_all$sum.crops)),]

# build GLMM explain sum of crop pollen with time and with islands as random effect
# extreme values in Gran Canaria. Add the following term to exclude
# pollen_all_crop$island != "Gran Canaria"
# plot pollen abundance 5000 years (affinity 300*300)
plot_endem <- ggplot(pollen_all_endem)+
  geom_point(aes(time,sum.endem), col = "darkblue", alpha = 0.5) +
  geom_smooth(aes(time,sum.endem), span = 1,
              se = T, colour = "darkblue")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  labs(x = "Time [Cal. years BP]", y = "Endemic pollen [%]")+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title = element_text(size = 12))+
  xlim(5000,0)

plot_endem




# number of unique taxa looked at in this analysis
length(unique(pollen_tax$tax)) # 1158
length(unique(pollen_tax_filtered$tax_filtered)) # 652

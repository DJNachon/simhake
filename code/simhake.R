
# Forecasting changes in the size-at-maturity of a commercially valuable demersal species under global warming
# 
# David José Nachón1*, Eduardo Ramírez-Romero2, A. Paz1, Marta Cousido-Rocha1, Francisco Izquierdo1, M. Grazia Pennino3, & Santiago Cerviño1
# 
# 1Instituto Español de Oceanografía (IEO, CSIC). Centro Oceanográfico de Vigo. Subida a Radio Faro, 50-52, 36390, Vigo (Pontevedra), Spain.
# 2Instituto de Ciencias Marinas de Andalucía (ICMAN-CSIC)
# 3Instituto Español de Oceanografía (IEO-CSIC). Sede Central. Calle del Corazón de María, 8, 28002 Madrid, Spain.
# *Corresponding author: David José Nachón, *david.nachon@ieo.csic.es


rm(list=ls())

# Packages  -----------------------------------------

library(data.table)
library(DT)
library(fields)
library(fitdistrplus)
library(gam)
library(ggplot2)
library(GGally)
library(gratia)
library(gridExtra)
library(Hmisc)
library(mgcv)
library(pander)
library(patchwork)
library(PerformanceAnalytics)
library(plyr)
library(ragg)
library(stats)
library(tidyverse)
library(tseries)
library(viridis)

source( 'code/functions.R')


directory <- getwd()


# L50 -----------------------------------------

L50A <- data.frame( readxl::read_excel( './data/L50.xlsx', sheet='Atlantic'))
L50M <- data.frame( readxl::read_excel( './data/L50.xlsx', sheet='Mediterranean'))

L50A$Ocean <- 'Atlantic Ocean'
L50A$Stock <- 'North European Atlantic Ocean'

L50A$Stock[ which( L50A$Area %in% c('ICES.8.c', 'ICES.8.c.9.a', 'ICES.9.a', 'Morocco'))] <- 
  'South European Atlantic Ocean'

L50M$Ocean <- L50M$Stock <- 'Mediterranean Sea'

cols <- c( 'Year', 'Year_min', 'Year_max', 'L50', 'L50_min', 'L50_max', 
           'Sex', 'Area', 'Ocean', 'Stock', 'Observations')

L50 <- rbind( L50A[ , cols], L50M[ , cols])


#### Areas & stock ----------------------------------

L50 <- L50 %>% arrange( Ocean, Area)

areas <- unique( L50$Area); areas


#### Decisions ---------------------------------

for(i in 1:nrow(L50)){
  if( is.na( L50$L50[i]) & !is.na( L50$L50_min[i]) & !is.na( L50$L50_max[i])){
    L50$L50[i] <- (L50$L50_min[i] + L50$L50_max[i]) /2}
  
  if( is.na( L50$Year[i]) & !is.na( L50$Year_min[i]) & !is.na( L50$Year_max[i])){
    L50$Year[i] <- (L50$Year_min[i] + L50$Year_max[i]) /2}
}

# conflict macro - micro
L50 <- L50[ -which( str_detect( L50$Observations, 'micro') == T), ]

# no L50
L50 <- L50[ -which( is.na(L50$L50)), ]


L50$Observations <- L50$Year_min <- L50$Year_max <- L50$L50_min <- L50$L50_max <- NULL



# Environmental variables ---------------------------------------

SST <- data.frame( readxl::read_excel( './data/Temp.xlsx', sheet='SST'))
SBT <- data.frame( readxl::read_excel( './data/Temp.xlsx', sheet='SBT'))

latlon <- data.frame( readxl::read_excel( './data/LatLon.xlsx'))


### Temperature ---------------------

temp <- list( SST = SST, SBT = SBT)

for (j in names(temp)){
  
  temp[[j]] <- temp[[j]] %>% select( where( ~!all(. == "NaN")))
  
  temp[[j]]$'ICES.3.a.4.a' <- tmeans( temp[[j]], latlon, c('ICES.3.a','ICES.4.a'))
  temp[[j]]$'ICES.6.7' <- tmeans( temp[[j]], latlon, c('ICES.7.j','ICES.7.h','ICES.7.g','ICES.7.k','ICES.7.c','ICES.7.e','ICES.7.b','ICES.6.a','ICES.6.b'))
  temp[[j]]$'ICES.6.7.8' <- tmeans( temp[[j]], latlon, c('ICES.7.j','ICES.7.h','ICES.7.g','ICES.7.k','ICES.7.c','ICES.7.e','ICES.7.f','ICES.7.b','ICES.8.c','ICES.8.b','ICES.8.d','ICES.8.a','ICES.6.a','ICES.6.b'))
  temp[[j]]$'ICES.7' <- tmeans( temp[[j]], latlon, c('ICES.7.j','ICES.7.h','ICES.7.f','ICES.7.g','ICES.7.k','ICES.7.c','ICES.7.e','ICES.7.b'))
  temp[[j]]$'ICES.7.8.a' <- tmeans( temp[[j]], latlon, c('ICES.7.j','ICES.7.h','ICES.7.g','ICES.7.k','ICES.7.c','ICES.7.e','ICES.7.b','ICES.7.f','ICES.8.a'))
  temp[[j]]$'ICES.7.b.c.k.j.g' <- tmeans( temp[[j]], latlon, c('ICES.7.j','ICES.7.g','ICES.7.k','ICES.7.c','ICES.7.b'))
  temp[[j]]$'ICES.8' <- tmeans( temp[[j]], latlon, c('ICES.8.c','ICES.8.b','ICES.8.d','ICES.8.a'))
  temp[[j]]$'ICES.8.a.b' <- tmeans( temp[[j]], latlon, c('ICES.8.b','ICES.8.a'))
  temp[[j]]$'ICES.8.a.b.c' <- tmeans( temp[[j]], latlon, c('ICES.8.c','ICES.8.b','ICES.8.a'))
  temp[[j]]$'ICES.8.a.b.d' <- tmeans( temp[[j]], latlon, c('ICES.8.b','ICES.8.d','ICES.8.a'))
  temp[[j]]$'ICES.8.c.9.a' <- tmeans( temp[[j]], latlon, c('ICES.8.c','ICES.9.a'))
  
  temp[[j]]$'GSA1.10' <- tmeans( temp[[j]], latlon, c('GSA1','GSA10'))
  temp[[j]]$'GSA6.7' <- tmeans( temp[[j]], latlon, c('GSA7','GSA6'))
  temp[[j]]$'GSA12.13' <- tmeans( temp[[j]], latlon, c('GSA12','GSA13'))
  temp[[j]]$'GSA12.13.14' <- tmeans( temp[[j]], latlon, c('GSA12','GSA13','GSA14'))
  temp[[j]]$'GSA12.14' <- tmeans( temp[[j]], latlon, c('GSA12','GSA14'))
  temp[[j]]$'GSA12.16' <- tmeans( temp[[j]], latlon, c('GSA12','GSA16'))
  temp[[j]]$'GSA15.16' <- tmeans( temp[[j]], latlon, c('GSA15','GSA16'))
  temp[[j]]$'GSA17.18' <- tmeans( temp[[j]], latlon, c('GSA17','GSA18'))
  temp[[j]]$'GSA22.27' <- tmeans( temp[[j]], latlon, c('GSA22','GSA27'))
  
}


#### Mean by year ------------------

temp_mean <- list()

for( j in names(temp)){
  temp_mean[[j]] <- temp[[j]][,-2]      # subset(temp[[j]], Month<=6)[,-2]
  temp_mean[[j]] <- aggregate(. ~ Year, data = temp_mean[[j]], FUN = mean, na.rm = TRUE)
}


### Lat & Lon -------------------------------------

latlon <- rbind( latlon,
                 c( Area = 'ICES.3.a.4.a', ameans( latlon, c('ICES.3.a', 'ICES.4.a'))),
                 c( Area = 'ICES.4', ameans( latlon, c('ICES.4.a', 'ICES.4.b', 'ICES.4.c'))),
                 c( Area = 'ICES.6', ameans( latlon, c('ICES.6.a', 'ICES.6.b'))),
                 c( Area = 'ICES.6.7', ameans( latlon, c('ICES.7.j','ICES.7.h','ICES.7.g','ICES.7.k','ICES.7.c','ICES.7.e','ICES.7.b','ICES.6.a','ICES.6.b'))),
                 c( Area = 'ICES.6.7.8', ameans( latlon, c('ICES.7.j','ICES.7.h','ICES.7.g','ICES.7.k','ICES.7.c','ICES.7.e','ICES.7.f','ICES.7.b','ICES.8.c','ICES.8.b','ICES.8.d','ICES.8.a','ICES.6.a','ICES.6.b'))),
                 c( Area = 'ICES.7', ameans( latlon, c('ICES.7.a', 'ICES.7.b', 'ICES.7.c', 'ICES.7.d', 'ICES.7.e', 'ICES.7.f', 'ICES.7.g', 'ICES.7.h', 'ICES.7.j', 'ICES.7.k'))),
                 c( Area = 'ICES.7.b.c.k.j.g', ameans( latlon, c('ICES.7.b', 'ICES.7.c', 'ICES.7.g','ICES.7.j', 'ICES.7.k'))),
                 c( Area = 'ICES.7.8.a', ameans( latlon, c('ICES.7.j','ICES.7.h','ICES.7.g','ICES.7.k','ICES.7.c','ICES.7.e','ICES.7.b','ICES.7.f','ICES.8.a'))),
                 c( Area = 'ICES.8', ameans( latlon, c('ICES.8.a', 'ICES.8.b', 'ICES.8.c', 'ICES.8.d', 'ICES.8.e'))),
                 c( Area = 'ICES.8.a.b', ameans( latlon, c('ICES.8.a', 'ICES.8.b'))),
                 c( Area = 'ICES.8.a.b.c', ameans( latlon, c('ICES.8.a', 'ICES.8.b','ICES.8.c'))),
                 c( Area = 'ICES.8.a.b.d', ameans( latlon, c('ICES.8.a', 'ICES.8.b','ICES.8.d'))),
                 c( Area = 'ICES.8.c.9.a', ameans( latlon, c('ICES.9.a', 'ICES.8.c'))),
                 c( Area = 'GSA1.10', ameans( latlon, c('GSA1', 'GSA10'))),
                 c( Area = 'GSA.5.6', ameans( latlon, c('GSA5', 'GSA6'))),
                 c( Area = 'GSA6.7', ameans( latlon, c('GSA6', 'GSA7'))),
                 c( Area = 'GSA12.13', ameans( latlon, c('GSA12', 'GSA13'))),
                 c( Area = 'GSA12.13.14', ameans( latlon, c('GSA12', 'GSA13', 'GSA14'))),
                 c( Area = 'GSA12.14', ameans( latlon, c('GSA12', 'GSA14'))),
                 c( Area = 'GSA12.16', ameans( latlon, c('GSA12', 'GSA16'))),
                 c( Area = 'GSA15.16', ameans( latlon, c('GSA15', 'GSA16'))),
                 c( Area = 'GSA17.18', ameans( latlon, c('GSA17', 'GSA18'))),
                 c( Area = 'GSA22.27', ameans( latlon, c('GSA22', 'GSA27'))))



## Match -------------------------------------------

L50$Lat <- L50$Lon <- NA

for( i in 1:nrow(latlon)){
  
  iar <- latlon$Area[i]
  ind <- which( L50$Area == iar)
  
  L50$Lat[ind] <- latlon$Lat[i]    
  L50$Lon[ind] <- latlon$Lon[i]
  
}


## Save for map

save( L50, file = './map/Data/L50/L50.RData')


# No year remove

fullL50 <- L50

L50 <- L50[ -which( is.na( L50$Year)), ]
L50$Year <- round( L50$Year, 0)



for (j in names(temp_mean)){
  
  L50[ ,paste0(j)] <- L50[ ,paste0(j,'_1')] <- L50[ ,paste0(j,'_2')] <- NA
  
  for( i in 1:nrow(L50)){
    
    i_year <- L50$Year[i]
    i_area <- L50$Area[i]
    
    ind <- ifelse( any( temp_mean[[j]]$Year == i_year) == TRUE,
                   which( temp_mean[[j]]$Year == i_year), 0)
    
    if( ind >= 2){ ind_1 <- ind-1} else { ind_1 <- 0}
    if( ind >= 3){ ind_2 <- ind-2} else { ind_2 <- 0}
    
    if( is.null( temp_mean[[j]][ ind, i_area])){
      
      L50[ i, paste0(j)] <- L50[ i, paste0(j,'_1')] <-L50[ i, paste0(j,'_2')] <- NA
      
    } else {
      
      if( ind >= 1){ L50[ i, paste0(j)] <- temp_mean[[j]][ ind, i_area]}
      if( ind_1 >= 1){ L50[ i, paste0(j,'_1')] <- temp_mean[[j]][ ind_1, i_area]}
      if( ind_2 >= 1){ L50[ i, paste0(j,'_2')] <- temp_mean[[j]][ ind_2, i_area]}
      
    }
  }
}


## Save ----------------------------------------------

save( L50, SST, SBT, latlon, fullL50, file = './data/data.RData')



# Exploratory analysis -----------------------------------------

mL50 <- subset( L50, Sex == 1)
fL50 <- subset( L50, Sex == 2)

L50$Sex <- factor( L50$Sex); levels( L50$Sex) <- c( 'Males', 'Females', 'Combined')
L50$Stock <- factor( L50$Stock, levels = c('North European Atlantic Ocean',                
                                           'South European Atlantic Ocean', 'Mediterranean Sea'))
L50$Area <- as.factor(L50$Area)
L50$Ocean <- as.factor(L50$Ocean)



## Density plot -----------------------------------------

max_Value <- max( L50$L50)

fullL50$Sex <- factor(fullL50$Sex); levels(fullL50$Sex) <- c('Males', 'Females', 'Combined')
fullL50$Stock <- factor( fullL50$Stock, levels = c('North European Atlantic Ocean',                
                                                   'South European Atlantic Ocean', 'Mediterranean Sea'))

DensityPlot <- ggplot( fullL50, aes( L50, fill = Sex)) +
  geom_density(alpha = 0.8, color = NA) + theme_light() +
  theme(axis.text.y = element_text(face = 'plain', colour = 1, size = rel(1)),
        legend.position = 'bottom', legend.justification = 'center',
        legend.box.margin = margin(0, 0, 0, 0), legend.background = element_rect(fill = 'transparent', color = NA),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +  
  scale_fill_manual(values = c('#FF6666', '#66CC99', '#6699CC')) +
  labs(x = 'Size-at-maturity (cm)', y = 'Density') + xlim(0, max_Value) +
  facet_wrap(~ Stock, nrow = 3, scales = 'fixed')

DensityPlot

plotdir <- paste0(directory, '/exploratory_plots/')
dir.create( path = plotdir, showWarnings = TRUE, recursive = TRUE)

ggsave( paste0(plotdir, 'DensityPlot_L50byArea.tif'),
        plot = DensityPlot, width = 8, height = 6, dpi = 300) 


## Descriptive statistics ---------------

### by sex and stock

summary_stats <- dplyr::group_by(fullL50, Sex, Stock) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(L50, na.rm = TRUE),
    sd = sd(L50, na.rm = TRUE),
    min = min(L50, na.rm = TRUE),
    max = max(L50, na.rm = TRUE)
  )

print(summary_stats)

summary_stats <- dplyr::group_by(fullL50, Stock) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(L50, na.rm = TRUE),
    sd = sd(L50, na.rm = TRUE),
    min = min(L50, na.rm = TRUE),
    max = max(L50, na.rm = TRUE)
  )

print(summary_stats)



## Statistics ---------------------------------

shapiro <- array( NA, dim = c( length(unique(fullL50$Sex)), length(unique(fullL50$Stock))),
                  dimnames = list( unique(fullL50$Sex), unique(fullL50$Stock)))

for( j in unique(fullL50$Sex)) for( i in unique( fullL50$Stock)) 
  shapiro[j,i] <- shapiro.test( subset( fullL50, Stock == i & Sex == j)$L50)$p.value

shapiro

kruskal_stock <- rep( NA, length(unique(fullL50$Stock)))
names(kruskal_stock) <- unique(fullL50$Stock)

kruskal_sex <- rep( NA, length(unique(fullL50$Sex)))
names(kruskal_sex) <- unique(fullL50$Sex)

for(i in unique(fullL50$Stock))
  kruskal_stock[i] <- kruskal.test(L50 ~ Sex, data = subset( fullL50, Stock == i))$p.value

for(i in unique(fullL50$Sex))
  kruskal_sex[i] <- kruskal.test(L50 ~ Stock, data = subset( fullL50, Sex == i))$p.value

kruskal_stock
kruskal_sex


pairwise_stock <- pairwise_sex <- list()

for(i in unique(fullL50$Stock)){
  il50 <- subset( fullL50, Stock == i)
  pairwise_stock[[i]] <- pairwise.wilcox.test(il50$L50, il50$Sex, p.adjust.method = "bonferroni")$p.value
}

pairwise_stock

for(i in unique(fullL50$Sex)){
  il50 <- subset( fullL50, Sex == i)
  pairwise_sex[[i]] <- pairwise.wilcox.test(il50$L50, il50$Stock, p.adjust.method = "BH")$p.value
}

pairwise_sex



## Count plot ------------------------------------

TotalYearPlotbyStock <- ggplot( L50, aes( x = Year, color = Sex, fill = Sex)) +
  geom_histogram(position = 'stack', alpha = 0.8) + labs(x = 'Year', y = 'Records number') + 
  theme_light() + theme( legend.position = 'bottom', legend.title = element_text(size = 12),
                         axis.text = element_text(size = 12), axis.title = element_text(size = 12), strip.text = element_text(size = 12)) +
  scale_color_manual(values = c('#FF6666', '#66CC99', '#6699CC')) +
  scale_fill_manual(values = c('#FF6666', '#66CC99', '#6699CC')) +  
  facet_wrap(~ Stock, nrow = 3, scales = 'fixed')

TotalYearPlotbyStock

ggsave( paste0( plotdir, filename = 'Histogram_L50bySex.tif'),
        plot = TotalYearPlotbyStock, width = 8, height = 6, dpi = 300)


## L50 vs SST, SBT and Year ------------------------------------------

color_palette <- viridis(3)

### SST
SSTStock <- ggplot( L50, aes( x = SST, y = L50)) +
  geom_point(aes(color = Stock), shape = 19, size = 3, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE, option = 'viridis') +
  labs(x = 'SST', y = 'Size-at-maturity (cm)', color = 'Geographic Area') +
  theme_light()

### SBT
SBTStock <- ggplot( L50, aes( x = SBT, y = L50)) +
  geom_point(aes(color = Stock), shape = 19, size = 3, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE, option = 'viridis') +
  labs(x = 'SBT', y='', color = 'Geographic Area') +
  theme_light()


combined_plot <- (SSTStock + SBTStock) +
  plot_layout(ncol = 2, guides = 'collect') &
  theme(legend.position = 'bottom')

print(combined_plot)

ggsave( paste0( plotdir, 'L50by_SST_SBT.tif'),
        plot = combined_plot, width = 8, height = 6, dpi = 300)


## Correlation --------------------

shapiro.test(L50$L50)
shapiro.test(L50$SST)
shapiro.test(L50$SBT)

cor.test(L50$L50, L50$SST, method = "spearman")

L50_clean <- na.omit(L50)
cor.test(L50_clean$L50, L50_clean$SBT, method = "spearman")


# GAMs analysis  -----------------------------------------

### Subset Year > 1950

allL50 <- L50
L50 <- subset( L50, Year >= 1958)


## Step 0 -----------------------------------------

gam0 <- gam( L50 ~ s( SBT, k = 4), method = 'REML', data = L50) 
gam_res( gam0)
gam_check( gam0)


### Non-retained combinations

gam0_1 <- gam( L50 ~ s( SST,k=4), method='REML', data=L50); gam_res( gam0_1)
gam0_2 <- gam( L50 ~ SST, method='REML', data=L50); gam_res( gam0_2)
gam0_3 <- gam( L50 ~ SBT, method='REML', data=L50); gam_res( gam0_3)
gam0_4 <- gam( L50 ~ Sex, method='REML', data=L50); gam_res( gam0_4)
gam0_5 <- gam( L50 ~ Ocean, method='REML', data=L50); gam_res( gam0_5)
gam0_6 <- gam( L50 ~ Stock, method='REML', data=L50); gam_res( gam0_6)
gam0_7 <- gam( L50 ~ te( Lon,Lat), method='REML', data=L50); gam_res( gam0_7)

LLvec <- c( min(L50$Lon, na.rm=T), max(L50$Lon, na.rm=T), 
            min(L50$Lat, na.rm=T), max(L50$Lat, na.rm=T))

gam_map( gam0_7, LLvec)


## Step 1 -----------------------------------------

gam1 <- gam( L50 ~ s( SBT, k = 4) + Sex, method = 'REML', data = L50) 
gam_res( gam1)
gam_check( gam1)

### Non-retained combinations
gam1_1 <- gam( L50 ~ s(SBT,k=4) + Ocean, method='REML', data=L50); gam_res( gam1_1)
gam1_2 <- gam( L50 ~ s(SBT,k=4) + Stock, method='REML', data=L50); gam_res( gam1_2)
gam1_3 <- gam( L50 ~ s(SBT,k=4) + te(Lon,Lat), method='REML', data=L50); gam_res( gam1_1)

gam_map( gam1_3, LLvec)

## Step 2 -----------------------------------------

gam2 <- gam( L50 ~ s( SBT, k = 4) + Sex + te(Lon,Lat), method='REML' , data = L50) 
gam_res( gam2)
gam_check( gam2)
gam_map( gam2, LLvec)


### Step 2. Non-retained combinations 
gam2_1 <- gam( L50 ~ s(SBT,k=4) + Sex + Ocean, method='REML', data=L50); gam_res( gam2_1)
gam2_2 <- gam( L50 ~ s(SBT,k=4) + Sex + Stock, method='REML', data=L50); gam_res( gam2_2)


# Final GAM -----------------------------------------

plot( gam2)
o
## Lineal SBT

## Decision ----------------------------------

gamf <- gam( L50 ~ SBT + Sex + te(Lon,Lat), method='REML', data = L50) 
gam_res( gamf)
gam_check( gamf)
gam_map( gamf, LLvec)

winnerT <- names(gamf$coefficients)[2]; winnerT   # 'SBT'


## With SST --------------------------------------------------

sstgam0 <- gam( L50 ~ s( SST, k = 4), method = 'REML', data = L50) 
gam_res( sstgam0)
gam_check( sstgam0)

sstgam1 <- gam( L50 ~ s( SST, k = 4) + Sex, method = 'REML', data = L50) 
gam_res( sstgam1)
gam_check( sstgam1)

sstgam2 <- gam( L50 ~ s( SST, k = 4) + Sex + te(Lon,Lat), method='REML' , data = L50) 
gam_res( sstgam2)
gam_check( sstgam2)
gam_map( sstgam2, LLvec)

plot( sstgam2)
o

sstgamf <- gam( L50 ~ SST + Sex + te(Lon,Lat), method='REML', data = L50) 
gam_res( sstgamf)
gam_check( sstgamf)
gam_map( sstgamf, LLvec)

## SBT fit plot -------------------------

colors <- c('Males' = '#FF6666', 'Combined' = '#6699CC', 'Females' = '#66CC99')

coef <- gamf$coefficients

intercept_males <- coef[1]
intercept_combined <- coef[1] + coef['SexCombined']
intercept_females <- coef[1] + coef['SexFemales']

GAMSBTplot <- ggplot( L50, aes(x = SBT, y = L50, color = Sex)) +
  geom_point(shape = 19, size = 3, alpha = 0.7) +
  scale_color_manual(values = colors) +   theme_light() +
  labs(x = 'SBT', y = 'Size-at-maturity (cm)', color = 'Sex') +
  theme( legend.position = 'bottom',
         axis.text = element_text(size = 12), axis.title = element_text(size = 14),
         legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  # Agregar líneas de ajuste
  geom_abline(intercept = intercept_males, slope = coef[2], linetype = 'solid', color = colors['Males']) +
  geom_abline(intercept = intercept_combined, slope = coef[2], linetype = 'solid', color = colors['Combined']) +
  geom_abline(intercept = intercept_females, slope = coef[2], linetype = 'solid', color = colors['Females'])

GAMSBTplot

ggsave( paste0( plotdir, filename = 'GAM_SBTfit.tif'), plot = GAMSBTplot,
        width = 8, height = 6, dpi = 300)



# L50 proyections ---------------------------------------

tempnames <- c( 'SST', 'SBT')

tproj <- list()
tproj_list <- list()
tPTs <- list()


for( tn in tempnames){
  
  ## Choose T ----------------------

  proj_data <- ifelse( tn == 'SST',
                       paste0(directory, '/data/SST.xlsx'), 
                       paste0(directory, '/data/SBT.xlsx'))
  
  Historical <- data.frame(readxl::read_excel( proj_data, sheet = 'historico'))
  RCP4.5 <- data.frame(readxl::read_excel( proj_data, sheet = 'rcp45'))
  RCP8.5 <- data.frame(readxl::read_excel( proj_data, sheet = 'rcp85'))
  
  PTs <- list( Historical = Historical, RCP4.5 = RCP4.5, RCP8.5 = RCP8.5)
  
  areas <- c( 'ICES.4', 'ICES.7', 'ICES.8.c.9.a', 'GSA.5.6')
  
  scenarios <- names( PTs)
  
  
  # Scenarios ----------------------
  
  oldPTs <- PTs
  
  for( sc in scenarios){ 
    
    PTs[[sc]] <- oldPTs[[sc]][,-2]      # subset( oldPTs[[sc]], Month <= 6)[,-2]
    
    PTs[[sc]]$'ICES.8.c.9.a' <- tmeans( PTs[[sc]], latlon, c('ICES.8.c','ICES.9.a'))
    PTs[[sc]]$'GSA.5.6' <- tmeans( PTs[[sc]], latlon, c('GSA5','GSA6'))
    
    PTs[[sc]]$'ICES.8.c.9.a_std' <- tmeans( PTs[[sc]], latlon, c('ICES.8.c_std','ICES.9.a_std'), c('ICES.8.c','ICES.9.a'))
    PTs[[sc]]$'GSA.5.6_std' <- tmeans( PTs[[sc]], latlon, c('GSA5_std','GSA6_std'), c('GSA5','GSA6'))
    
    if( sc == 'Historical'){
      
      PTs[[sc]]$'GSA.5.6_std' <- PTs[[sc]]$'ICES.4_std' <- 
        PTs[[sc]]$'ICES.7_std' <- PTs[[sc]]$'ICES.8.c.9.a_std' <- 0
      
    }
    
    PTs[[sc]] <- aggregate(. ~ Year, data = PTs[[sc]], FUN = mean, na.rm = TRUE)
    PTs[[sc]] <- PTs[[sc]][, c('Year', areas, 'ICES.4_std', 'ICES.7_std', 'ICES.8.c.9.a_std',
                               'GSA.5.6_std')]
    
  }
  
  
  # Loop ----------------------------
  
  dummyy <- PTs[['Historical']][which(PTs[['Historical']]$Year==2022),]
  
  for( j in names(PTs)[ -which(names(PTs)=='Historical')]){ 
    PTs[[j]] <- rbind( dummyy, PTs[[j]])}
  
  
  allscenarios <- scenarios    # c( scenarios, scenarios2, scenarios3)
  
  proj <- data.frame()
  proj_list <- list()
  
  
  for( sc in allscenarios){
    for( ar in areas){
      
      Tdf <- PTs[[sc]][,c('Year', ar, paste0(ar,'_std'))]
      colnames(Tdf) <- c('Year', 'Temp', 'sd') 
      
      if( Tdf$Year[length(Tdf$Year)] >= 2050){
        fyi <- which(Tdf$Year == 2022); eyi <- which(Tdf$Year == 2099)
      } else {
        fyi <- 1; eyi <- which(Tdf$Year == 2022)}
      
      lpy <- nrow(Tdf)
      
      vec <- unique(L50$Sex)
      
      proj_list[[sc]] <- list()
      
      for( i in 1:length(vec)){
        
        newdf <- data.frame( Year = Tdf$Year[fyi:eyi],
                             Lat = latlon$Lat[which(latlon$Area == ar)],
                             Lon = latlon$Lon[which(latlon$Area == ar)],
                             Sex = vec[i], SBT = Tdf$Temp[fyi:eyi], 
                             SBT_sd = Tdf$sd[fyi:eyi],
                             Scenario = sc, Area = ar)
        
        if( tn == 'SST'){ 
          newdf$SST <- newdf$SBT
          newdf$SST_sd <- newdf$SBT_sd
          tmod <- sstgamf} else {tmod <- gamf}
        
        predictions <- predict( tmod, type = 'response', newdata = newdf, se.fit = T)
        newdf$L50 <- predictions$fit
        if( sc == 'Historical'){ newdf$L50_sd <- 0} else { newdf$L50_sd <- predictions$se.fit}
        
        proj <- rbind( proj, newdf)
        proj_list[[sc]][[vec[i]]] <- newdf
        
      }
    }
  }
  
  
  
  # Plots -----------------------------------------------------------
  
  
  if( winnerT == tn){
    plotd <- paste0( directory, '/proj_plots/')} else { 
      plotd <- paste0( directory, '/proj_plots/', tn, '/')}
  
  combo_fem_ann_sd <- comboplot( proj, path = plotd, name = 'females_ann_sd', 
    text = T, dif_mf = round( c( gamf$coefficients[3], gamf$coefficients[3]-gamf$coefficients[4]), 1), sd = T)
  
  if( tn == 'SST'){ proj$SBT <- proj$SBT_sd <- NULL}
  
  tproj[[tn]] <- proj
  tproj_list[[tn]] <- proj_list
  tPTs[[tn]] <- PTs
  
}


proj <- tproj[['SBT']]


# Save -----------------------------------

resdir <- paste0( directory, '/results')
dir.create( path = resdir, showWarnings = TRUE, recursive = TRUE)

save.image( paste0( resdir, '/results.RData'))

save( proj, file = paste0( resdir, '/L50.RData'))



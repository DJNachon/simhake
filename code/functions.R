tmeans <- function( temp, ll, areas, llareas = areas){ 
  
  mareas <- areas[ areas %in% colnames(temp)]
  mllareas <- llareas[ llareas %in% colnames(temp)]
  ill <- ll %>% filter( Area %in% mllareas)
  ill <- ill[ order(ill$Area),]
  ill <- matrix( ill$n / sum(ill$n), ncol = 1)
  
  it <- as.matrix( temp[ ,mareas])
  it <- it[, order(colnames(it))]
  
  wmean <- it%*%ill
  
  return( as.vector(wmean))
  
}


weighted_or_mean <- function(x, w) {
  if (any(is.na(w))) {
    return(mean(x, na.rm = TRUE))
  } else {
    return(weighted.mean(x, w, na.rm = TRUE))
  }
}



ameans <- function( df, areas) {
  
  idf <- df %>% filter( Area %in% areas) 
  
  total_n <- sum(idf$n)
  
  lat <- sum(idf$Lat * idf$n) / total_n
  lon <- sum(idf$Lon * idf$n) / total_n
  
  return( data.frame( Lat = lat, Lon = lon, n = total_n))
}


gam_check <- function( model){
  
  res <- c( jarquebera = jarque.bera.test( resid(model))$p.value,
            shapiro = shapiro.test( resid(model))$p.value,
            ttest = t.test( resid(model),mu=0)$p.value)
  
  p <- 1:20; for(i in 1:20){ 
    p[i] = Box.test (resid(model), lag = i, type = 'Ljung')$p.value}
  
  par(mfrow=c(3,1))
  plot( p,ylim=c(0,1),col='black',ylab='P-valor',xlab='Lag', main='Box test')
  abline( h=0.05,lty=2,col='blue')
  acf( resid(model),lag.max=38, main='ACF')
  pacf( resid(model),lag.max =38, main='PACF')
  par(mfrow=c(1,1))
  
  return(res)
}


gam_res <- function( model){
  res <- c( deviance = summary(model)$dev.expl, AIC = model$aic)
  return( res)}


gam_map <- function( model, LLvec){
  sm <- smooth_estimates( model, smooth = 'te(Lon,Lat)')
  draw(sm) +  theme( panel.background = element_blank(), axis.line = element_blank()) +
    geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),fill = 'white', 
                 color = 'black') + coord_fixed(ratio=1.6, xlim = c(LLvec[1],LLvec[2]),
                                                ylim=c(LLvec[3], LLvec[4]))}


comboplot <- function( proj, path = paste0( getwd(),'/'), name = 'plot', 
                       sex = 'Females', plotcolor = 'Scenario', 
                       text = FALSE, dif_mf = NULL,
                       sd = FALSE, alpha = 0.05){
  
  tname <- 'SST (ºC)'

  df <- proj %>%
      pivot_longer( cols = c(SST, L50), names_to = "Type", values_to = "Value") %>%
      mutate( sd = case_when( Type == "SST" ~ SST_sd, Type == "L50" ~ L50_sd),
              Type = recode( Type, SST = tname, L50 = "Size-at-maturity (cm)")) %>%
      mutate( Type = factor( Type, levels = c( tname, "Size-at-maturity (cm)")))
   
  df$Scenario <- factor( df$Scenario)
  df$Area <- factor( df$Area, levels = c( 'ICES.4', 'ICES.7', 'ICES.8.c.9.a', 'GSA.5.6'))
    
  colors <- c( 1,4,6)
  
  dir.create( path = path, showWarnings = FALSE, recursive = TRUE)
  
  arl <- subset( df, Sex %in% sex)
  
  if( length(sex) == 1){ 
    
    plL50 <- ggplot(data = arl, aes(x = Year, y = Value, color = Scenario, fill = Scenario)) +
      geom_line(size = 0.5) +
      geom_vline( xintercept = 2022, linetype = 'dashed', size = .2) +
      facet_grid( Type ~ Area, scales = 'free') +
      labs(x = 'Year', y = '') +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme( legend.position = 'bottom',
             panel.border = element_rect(color = 'gray', fill = NA, size = 0.5),
             panel.grid.major = element_line(color = 'lightgray', size = 0.3), 
             panel.grid.minor = element_line(color = 'lightgray', size = 0.15),
             panel.background = element_rect(fill = 'white', color = NA))
    
    if( sd == T){
      z <- qnorm( 1 - alpha/2)
      plL50 <- plL50 + 
        geom_ribbon( aes( ymin = Value - z * sd, ymax = Value + z * sd, 
                          fill = Scenario), linetype = 0, alpha = 0.15)}
    
    if( text == TRUE){
      plL50 <- plL50 + 
        geom_text( data = subset(arl, Type == "Size-at-maturity (cm)" & Area == "GSA.5.6"), 
            aes(x = Inf, y = Inf, label = 
            paste0( 'Only females represented \n(Males = Females - ', dif_mf[1],' cm; \n Combined = Females - ', dif_mf[2],' cm)')),
            color = "black", size = 2, hjust = 1.1, vjust = 1.2)}
    
  } else {
    
    if( plotcolor == 'Scenario'){
      
      plL50 <- 
        ggplot(data = arl, aes(x = Year, y = Value, 
                               color = Scenario, fill = Scenario, linetype = Sex)) +
        geom_line(size = 0.5) +
        geom_vline( xintercept = 2022, linetype = 'dashed', size = .2) +
        facet_grid( Type ~ Area, scales = 'free') +
        labs(x = 'Year', y = '') +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        theme( legend.position = 'bottom',
               panel.border = element_rect(color = 'gray', fill = NA, size = 0.5),
               panel.grid.major = element_line(color = 'lightgray', size = 0.3), 
               panel.grid.minor = element_line(color = 'lightgray', size = 0.15),
               panel.background = element_rect(fill = 'white', color = NA))
      
      if( sd == T){
        z <- qnorm( 1 - alpha/2)
        plL50 <- plL50 + 
          geom_ribbon( aes( ymin = Value - z * sd, ymax = Value + z * sd, 
                            fill = Scenario, linetype = Sex), colour = NA, alpha = 0.15)}
      
    } else {
      
      plL50 <- ggplot(data = arl, aes(x = Year, y = Value, linetype = Scenario)) +
        geom_line(data = subset(arl, Type == tname), 
                  aes(color = NA), size = 0.5, color = "black", show.legend = FALSE) +
        geom_line(data = subset(arl, Type != tname), 
                  aes(color = Sex), size = 0.5) +
        geom_vline(xintercept = 2022, linetype = 'dashed', size = .2) +
        facet_grid(Type ~ Area, scales = 'free') +
        labs(x = 'Year', y = '') +
        scale_color_manual( name = "Sex",
              values = c('Males' = '#FF6666', 'Combined' = '#6699CC', 'Females' = '#66CC99')) +
        theme(legend.position = 'bottom',
              panel.border = element_rect(color = 'gray', fill = NA, size = 0.5),
              panel.grid.major = element_line(color = 'lightgray', size = 0.3), 
              panel.grid.minor = element_line(color = 'lightgray', size = 0.15),
              panel.background = element_rect(fill = 'white', color = NA))
      
      if(sd == TRUE) {
        z <- qnorm(1 - alpha/2)
        alphasc <- c( 0, .1, .2)
        plL50 <- plL50 + 
          geom_ribbon(data = subset(arl, Type == tname), 
                      aes(ymin = Value - z * sd, ymax = Value + z * sd, 
                          fill = NA, alpha = Scenario, linetype = Scenario), size=.1, fill = "black", 
                      show.legend = FALSE) +
          geom_ribbon(data = subset(arl, Type != tname), 
                      aes(ymin = Value - z * sd, ymax = Value + z * sd, color = Sex,
                          fill = Sex, alpha = Scenario, linetype = Scenario), size = .1) +
          scale_fill_manual(name = "Sex",
                            values = c('Males' = '#FF6666', 'Combined' = '#6699CC', 'Females' = '#66CC99')) +
          scale_alpha_manual( values = alphasc)}
      
      
    }
  }
  
  
  print(plL50)
  
  ggsave( paste0( path, name, '.tif'), plot = plL50, 
          width = 8, height = 6, dpi = 300) 
  
  
  return( plL50)
  
}


comboplot2 <- function( proj, path = paste0( getwd(),'/'), name = 'plot', 
                       sex = 'Females', plotcolor = 'Scenario', 
                       text = FALSE, dif_mf = NULL,
                       sd = FALSE, alpha = 0.05){
  
  tname <- 'SBT (ºC)'
  
  df <- proj %>%
    pivot_longer( cols = c(SBT, L50), names_to = "Type", values_to = "Value") %>%
    mutate( sd = case_when( Type == "SBT" ~ SBT_sd, Type == "L50" ~ L50_sd),
            Type = recode( Type, SBT = tname, L50 = "Size-at-maturity (cm)")) %>%
    mutate( Type = factor( Type, levels = c( tname, "Size-at-maturity (cm)")))
  
  df$Scenario <- factor( df$Scenario)
  df$Area <- factor( df$Area, levels = c( 'ICES.4', 'ICES.7', 'ICES.8.c.9.a', 'GSA.5.6'))
  
  colors <- c( 1,4,6)
  
  dir.create( path = path, showWarnings = FALSE, recursive = TRUE)
  
  arl <- subset( df, Sex %in% sex)
  
  if( length(sex) == 1){ 
    
    plL50 <- ggplot(data = arl, aes(x = Year, y = Value, color = Scenario, fill = Scenario)) +
      geom_line(size = 0.5) +
      geom_vline( xintercept = 2022, linetype = 'dashed', size = .2) +
      facet_grid( Type ~ Area, scales = 'free') +
      labs(x = 'Year', y = '') +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme( legend.position = 'bottom',
             panel.border = element_rect(color = 'gray', fill = NA, size = 0.5),
             panel.grid.major = element_line(color = 'lightgray', size = 0.3), 
             panel.grid.minor = element_line(color = 'lightgray', size = 0.15),
             panel.background = element_rect(fill = 'white', color = NA))
    
    if( sd == T){
      z <- qnorm( 1 - alpha/2)
      plL50 <- plL50 + 
        geom_ribbon( aes( ymin = Value - z * sd, ymax = Value + z * sd, 
                          fill = Scenario), linetype = 0, alpha = 0.15)}
    
    if( text == TRUE){
      plL50 <- plL50 + 
        geom_text( data = subset(arl, Type == "Size-at-maturity (cm)" & Area == "GSA.5.6"), 
                   aes(x = Inf, y = Inf, label = 
                         paste0( 'Only females represented \n(Males = Females - ', dif_mf[1],' cm; \n Combined = Females - ', dif_mf[2],' cm)')),
                   color = "black", size = 2, hjust = 1.1, vjust = 1.2)}
    
  } else {
    
    if( plotcolor == 'Scenario'){
      
      plL50 <- 
        ggplot(data = arl, aes(x = Year, y = Value, 
                               color = Scenario, fill = Scenario, linetype = Sex)) +
        geom_line(size = 0.5) +
        geom_vline( xintercept = 2022, linetype = 'dashed', size = .2) +
        facet_grid( Type ~ Area, scales = 'free') +
        labs(x = 'Year', y = '') +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        theme( legend.position = 'bottom',
               panel.border = element_rect(color = 'gray', fill = NA, size = 0.5),
               panel.grid.major = element_line(color = 'lightgray', size = 0.3), 
               panel.grid.minor = element_line(color = 'lightgray', size = 0.15),
               panel.background = element_rect(fill = 'white', color = NA))
      
      if( sd == T){
        z <- qnorm( 1 - alpha/2)
        plL50 <- plL50 + 
          geom_ribbon( aes( ymin = Value - z * sd, ymax = Value + z * sd, 
                            fill = Scenario, linetype = Sex), colour = NA, alpha = 0.15)}
      
    } else {
      
      plL50 <- ggplot(data = arl, aes(x = Year, y = Value, linetype = Scenario)) +
        geom_line(data = subset(arl, Type == tname), 
                  aes(color = NA), size = 0.5, color = "black", show.legend = FALSE) +
        geom_line(data = subset(arl, Type != tname), 
                  aes(color = Sex), size = 0.5) +
        geom_vline(xintercept = 2022, linetype = 'dashed', size = .2) +
        facet_grid(Type ~ Area, scales = 'free') +
        labs(x = 'Year', y = '') +
        scale_color_manual( name = "Sex",
                            values = c('Males' = '#FF6666', 'Combined' = '#6699CC', 'Females' = '#66CC99')) +
        theme(legend.position = 'bottom',
              panel.border = element_rect(color = 'gray', fill = NA, size = 0.5),
              panel.grid.major = element_line(color = 'lightgray', size = 0.3), 
              panel.grid.minor = element_line(color = 'lightgray', size = 0.15),
              panel.background = element_rect(fill = 'white', color = NA))
      
      if(sd == TRUE) {
        z <- qnorm(1 - alpha/2)
        alphasc <- c( 0, .1, .2)
        plL50 <- plL50 + 
          geom_ribbon(data = subset(arl, Type == tname), 
                      aes(ymin = Value - z * sd, ymax = Value + z * sd, 
                          fill = NA, alpha = Scenario, linetype = Scenario), size=.1, fill = "black", 
                      show.legend = FALSE) +
          geom_ribbon(data = subset(arl, Type != tname), 
                      aes(ymin = Value - z * sd, ymax = Value + z * sd, color = Sex,
                          fill = Sex, alpha = Scenario, linetype = Scenario), size = .1) +
          scale_fill_manual(name = "Sex",
                            values = c('Males' = '#FF6666', 'Combined' = '#6699CC', 'Females' = '#66CC99')) +
          scale_alpha_manual( values = alphasc)}
      
      
    }
  }
  
  
  print(plL50)
  
  ggsave( paste0( path, name, '.tif'), plot = plL50, 
          width = 8, height = 6, dpi = 300) 
  
  
  return( plL50)
  
}



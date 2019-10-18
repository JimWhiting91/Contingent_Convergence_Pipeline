# General functions for slim analysis

# Extract values per measure
fill_matrix<-function(vec){
  
  # First divide vec
  vec<-vec/(0.05*length(unique(slim_dd$Gene_ID)))
  
  # Empty mat
  mat<-matrix(ncol=length(levels(slim_dd$Treatment_Demo)),
              nrow=length(levels(slim_dd$Treatment_Demo)))
  # Add diagonal
  for(i in 1:length(levels(slim_dd$Treatment_Demo))){
    mat[i,i]<-1
  }
  # Add values from vec
  for(i in 1:length(comparisons_to_make[1,])){
    mat[comparisons_to_make[1,i],comparisons_to_make[2,i]]<-vec[i]
  }
  return(mat)
}

# Functions
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

plot_cor_mat<-function(input_matrix,cor=FALSE){
  if(cor==TRUE){
    cormat <- round(cor((input_matrix)),2)
  } else {
    cormat <- round((input_matrix),2)
  }
  
  library(reshape2)
  melted_cormat <- melt(cormat)
  library(ggplot2)
  ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  upper_tri <- get_upper_tri(cormat)
  upper_tri
  
  # Melt the correlation matrix
  library(reshape2)
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Heatmap
  library(ggplot2)
  ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  # Reorder the correlation matrix
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  # Print the heatmap
  #print(ggheatmap)
  
  ggheatmap + 
    #geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.title = element_text(size=24),
      legend.text = element_text(size=12),
      legend.position = c(0.4, 0.7),
      axis.text = element_text(size=16),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
}

reorder_cormat <- function(cormat){
  # Mirror matrix
  cormat[lower.tri(cormat)]<-t(cormat)[lower.tri(cormat)]
  dd <- dist(cormat)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


# Proportion Matrices
plot_prop_matrix<-function(mat){
  
  reorder_cormat <- function(cormat){
    # Mirror matrix
    cormat[lower.tri(cormat)]<-t(cormat)[lower.tri(cormat)]
    dd <- dist(cormat)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  mat[mat==1]<-NA
  
plot1<-plot_cor_mat(reorder_cormat(mat))
plot2<-plot1 +
  scale_fill_gradient2(low = "white", high = "red",
                       limit=c(0,max(mat)), space = "Lab", 
                       name="Proportion\nOverlap")+
  theme(legend.text = element_text(angle = 45,vjust=1))

return(plot2)
}

# ------------------------------------------------------ 
# Downsample the dataset for random top selection genes
# ------------------------------------------------------

random_downsample<-function(data){
  
  # Take only one treatment and one generation to remove repetition
  data2<-data[data$Treatment_Demo == data$Treatment_Demo[1] & data$Generation == "Gen_1",]
  
  # Randomise
  tmp<- data2[sample(1:nrow(data2)), ]
  
  # For the strongest selection coefficients, downsample relative to strength
  tmp1<-tmp[tmp$Selection_Coef %in% c(0.1,0.2,1.0),]
  
  tmp2<-tmp[tmp$Selection_Coef == 2.0,]
  tmp2<-tmp2[1:(length(tmp2$Gene_ID)/2),] # We half the size of this
  
  tmp3<-tmp[tmp$Selection_Coef == 10.0,]
  tmp3<-tmp2[1:(length(tmp3$Gene_ID)/10),] # We take 10$ the size of this
  
  # Rbind back together and sort
  to_keep<-rbind(tmp1,tmp2,tmp3)
  to_keep<-unique(to_keep$Gene_ID)
  tmp_out<-data[data$Gene_ID %in% to_keep,]
  
  # Return
  return(tmp_out)
}

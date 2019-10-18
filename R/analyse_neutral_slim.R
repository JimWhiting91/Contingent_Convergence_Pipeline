##########
# This script takes the output files from slim (1 per gene window, 20 iterations per 16 treatments per 12 sampling points) and analyses
# SPECIFICALLY THIS SCRIPT REPLACES THE PHENOdiv DATA WITH NEUTRAL AND REPEATS ALL MAJOR ANALYSES
# STEP 1 - Correlations

# Data = Spatial model, fixed DXY, fixed evolv_gen
##########

# Packages
lib=c("ggedit","tidyr","dplyr","pbapply","randomForest","VennDiagram","ghibli","ggpubr","data.table","ggplot2","parallel","doParallel","effects","knitr","lme4","MASS","lmerTest","dplyr","rpart","randomForest","MuMIn")
lapply(lib,library,character.only=T)

# Source functions
source("R/Slim_General_Functions.R")

##### STEP 1 #####

# Define inputs
pop_size<-1000

# Run over mut rates
data_files<-c("MutRate-6")
neutral_files<-c("Neutral_mut-6")


# Make output/figs directories
dir.create(file.path(paste0("outputs/",data_files,"_Neutral")))
dir.create(file.path(paste0("figs/",data_files,"_Neutral")))

# Better
fig_out<-paste0("figs/",data_files,"_Neutral")
output_out<-paste0("outputs/",data_files,"_Neutral")

# Read in all the data in the dir
to_read<-list.files(paste0("data/",data_files))
to_read_divergent<-to_read[grep("SLIM.out",to_read)]
to_read_null<-to_read[grep("NULL.out",to_read)]
to_read_neutral<-list.files(paste0("data/",neutral_files))

# Read in null and divergent separately
divergent_dd<-data.frame(rbindlist(lapply(to_read_divergent,function(x){return(fread(paste0("data/",data_files,"/",x)))})))
divergent_dd$run_type<-"True"
null_dd<-data.frame(rbindlist(lapply(to_read_null,function(x){return(fread(paste0("data/",data_files,"/",x)))})))
null_dd$run_type<-"Null"
neutral_dd<-data.frame(rbindlist(lapply(to_read_neutral,function(x){return(fread(paste0("data/",neutral_files,"/",x)))})))
neutral_dd$run_type<-"Neutral"

# Add in dummy columns for neutral
neutral_dd$Pop1_Mean<-NA
neutral_dd$Pop2_Mean<-NA
neutral_dd$Pop1_sd<-NA
neutral_dd$Pop2_sd<-NA
neutral_dd$Evolv_Gen<-NA

# Reorder all the columns
divergent_dd2<-divergent_dd[,order(colnames(divergent_dd))]
null_dd2<-null_dd[,order(colnames(null_dd))]
neutral_dd2<-neutral_dd[,order(colnames(divergent_dd))]

# Combine
slim_dd<-rbind(divergent_dd2,null_dd2,neutral_dd2)

# No Selection
slim_dd$Treatment_Demo<-paste0("Bot=",slim_dd$Bottleneck,"/Pop2=",slim_dd$Pop2_Size,"/Mig=",slim_dd$Migration)
slim_dd$Treatment_Demo<-as.factor(slim_dd$Treatment_Demo)

# Calculate deltaPI which is log(10)-transformed
# Add the bare minimum increase to PI, which is 1 more site in 1 more individual ie. 1/pop_size
slim_dd$deltaPI<-log10((slim_dd$Pop1_PI+(1/pop_size))/(slim_dd$Pop2_PI+(1/pop_size)))

# Transform selection_coef so that it scales positively with selection
slim_dd$Selection_Coef<--1*(log10(slim_dd$Selection_Coef))

# Add in gene info detail
gene_info<-read.table("data/slim_gene_info.txt")
colnames(gene_info)<-c("ID","Length","Exon_N","Exon_Length")
slim_dd$gene_length<-NA
slim_dd$exon_N<-NA
slim_dd$Target_sel<-NA
# loop to fill
for(i in 1:100){
  slim_dd[slim_dd$Gene_ID == i,"gene_length"]<-gene_info[gene_info$ID == i, "Length"]  
  slim_dd[slim_dd$Gene_ID == i,"exon_N"]<-gene_info[gene_info$ID == i, "Exon_N"]  
  slim_dd[slim_dd$Gene_ID == i,"Target_sel"]<-gene_info[gene_info$ID == i, "Exon_Length"]/25000
}

# Remove old, useless column
slim_dd<-slim_dd[,!(colnames(slim_dd) %in% c("Exon_Length_Total"))]

# Transform gen times to numbers rather than factors
sample_times<-c(100,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
slim_dd$Gen_N<-NA
for(i in 1:12){
  slim_dd[slim_dd$Generation == paste0("Gen_",i),"Gen_N"]<-sample_times[i]
}

# Separate out Nulls - HERE WE CHANGE THE DATA THATS GIVEN TO ALL MAJOR ANALYES TO PHENOnull
slim_neutral<-slim_dd[slim_dd$run_type == "Neutral",]
slim_dd<-slim_dd[slim_dd$run_type == "Neutral",]


#### STEP 2 ####
# What is the correlation of FST and DXY across different treatments?
cor_dd<-data.frame(Treatment = levels(slim_dd$Treatment_Demo))

# Calculate correlations over each generation
cor_dd_list<-list()
for (i in 1:length(cor_dd$Treatment)){
  # Treatments
  cor_dd$Mig[i]<-mean(slim_dd[slim_dd$Iteration==1 & slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"Migration"])
  cor_dd$Bot[i]<-mean(slim_dd[slim_dd$Iteration==1 & slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"Bottleneck"])
  cor_dd$Pop2[i]<-mean(slim_dd[slim_dd$Iteration==1 & slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"Pop2_Size"])
}


#-------------------------------------------
# EFFECT OF TREATMENTS ON MEASURES
#-------------------------------------------

# Plot relationships across generations
# Figures to examine effects of demographic variables across generation
treatment_vec<-c("Bottleneck","Pop2_Size","Migration")

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  ###################################################################################
  # This gets run over each iteration and averaged over
  
  iteration_list<-lapply(1:20,function(iter){
    # Edit the data such that Fst and dXY are same column with identifier column
    # Summarise
    tmp1<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","FST",treatment_vec[x],"Treatment_Demo")],Measure=rep("FST"))
    # Remove any NA values from calcs
    tmp1<-na.omit(tmp1)
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise(mean(FST),(sd(FST)/sqrt(length(FST)))))
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Mean","SE")
    
    
    # Average again...
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]]) %>%
                          dplyr::summarise(Mean=mean(Mean),
                                           SE=mean(SE)))
    
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Mean","SE")
    tmp1$Metric<-rep("FST")
    
    # Summarise
    tmp2<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","dXY",treatment_vec[x],"Treatment_Demo")],Measure=rep("dXY"))
    # Remove any NA values from calcs
    tmp2<-na.omit(tmp2)
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise(mean(dXY),(sd(dXY)/sqrt(length(dXY)))))
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Mean","SE")
    
    
    # Average again...
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]]) %>%
                          dplyr::summarise(Mean=mean(Mean),
                                           SE=mean(SE)))
    
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Mean","SE")
    tmp2$Metric<-rep("dXY")
    
    # Summarise
    tmp3<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","deltaPI",treatment_vec[x],"Treatment_Demo")],Measure=rep("deltaPI"))
    # Remove any NA values from calcs
    tmp3<-na.omit(tmp3)
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise(mean(deltaPI),(sd(deltaPI)/sqrt(length(deltaPI)))))
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Mean","SE")
    
    # Average again...
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]]) %>%
                          dplyr::summarise(Mean=mean(Mean),
                                           SE=mean(SE)))
    
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Mean","SE")
    tmp3$Metric<-rep("deltaPI")
    
    
    # Rbind
    tmp<-rbind(tmp1,tmp2,tmp3)
    tmp[,2]<-as.factor(tmp[,2])
    
    tmp$Metric<-factor(tmp$Metric,levels=c("FST","dXY","deltaPI"))
    levels(tmp$Metric)<-c("F[ST]","D[XY]","Delta*pi")
    return(tmp)
  })
  
  # Average over the list
  tmp<-Reduce("+", iteration_list) / length(iteration_list)
  tmp[,2]<-iteration_list[[1]][,2]
  tmp[,5]<-iteration_list[[1]][,5]
  
  ###############################################################
  
  # Plot
  pd <- position_dodge(0.5)
  
  g1<-ggplot(tmp,aes(x=Gen_N,y=Mean,shape=tmp[,2],colour=tmp[,2]))+
    #geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=pd,show.legend = FALSE) +
    geom_point(position=pd,size=4)+
    geom_line(show.legend = F)+
    facet_wrap(~Metric,ncol=1,scales = "free",labeller = label_parsed)+
    theme_bw()+
    theme(axis.title= element_text(size=18),
          axis.text=element_text(size=16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          strip.text = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "top")+
    scale_x_continuous(labels=c(0.1,2,4,6,8,10),breaks=c(100,2000,4000,6000,8000,10000))+
    xlab(expression(Generation~(10^3)))+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium"))) +
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab("Divergence")
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="DP Pop Size",shape="DP Pop Size")
  }
  
  return(g1)
})

pdf(paste0(fig_out,"/Spatial_Model_Divergence_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height=10)
print(ggarrange(plot_list[[1]],
                plot_list[[2]],
                plot_list[[3]],ncol = 3,labels = "AUTO"))
dev.off()

#-------------------------------------------
# CORRELATIONS WITH SELECTION
#-------------------------------------------

# Figures to examine effects of demographic variables across Gen_Ns on correlation coefficients with selection
treatment_vec<-c("Bottleneck","Pop2_Size","Migration")

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  ###################################################################################
  # This gets run over each iteration separately and averaged over
  
  iteration_list<-lapply(1:20,function(iter){
    
    # Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column
    # Summarise FST
    tmp1<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","FST","Selection_Coef",treatment_vec[x],"Treatment_Demo")],Measure=rep("FST"))
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(FST,Selection_Coef)$estimate),
                                           (cor.test(FST,Selection_Coef)$conf.int[1]),
                                           (cor.test(FST,Selection_Coef)$conf.int[2])))
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    tmp1$Metric<-"FST"
    
    # Summarise 
    tmp2<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","dXY","Selection_Coef",treatment_vec[x],"Treatment_Demo")],Measure=rep("dXY"))
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(dXY,Selection_Coef)$estimate),
                                           (cor.test(dXY,Selection_Coef)$conf.int[1]),
                                           (cor.test(dXY,Selection_Coef)$conf.int[2])))
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    tmp2$Metric<-"dXY"
    
    # Summarise 
    tmp3<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","deltaPI","Selection_Coef",treatment_vec[x],"Treatment_Demo")],Measure=rep("deltaPI"))
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(deltaPI,Selection_Coef)$estimate),
                                           (cor.test(deltaPI,Selection_Coef)$conf.int[1]),
                                           (cor.test(deltaPI,Selection_Coef)$conf.int[2])))
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    tmp3$Metric<-"deltaPI"
    
    
    # Rbind
    tmp<-rbind(tmp1,tmp2,tmp3)
    tmp[,2]<-as.factor(tmp[,2])
    
    tmp$Metric<-factor(tmp$Metric,levels=c("FST","dXY","deltaPI"))
    levels(tmp$Metric)<-c("F[ST]","D[XY]","Delta*pi")
    
    return(tmp)
  })
  
  # Average over the list
  tmp<-Reduce("+", iteration_list) / length(iteration_list)
  tmp[,2]<-iteration_list[[1]][,2]
  tmp[,6]<-iteration_list[[1]][,6]
  
  # Plot
  pd <- position_dodge(0.5)
  
  ###############################################################
  
  
  g1<-ggplot(tmp,aes(x=Gen_N,y=Cor.Coef,shape=tmp[,2],colour=tmp[,2]))+
    #geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, position=pd,show.legend = FALSE) +
    geom_line(show.legend = F)+
    geom_point(position=pd,size=3)+
    facet_wrap(~Metric,ncol=1,scales = "free",labeller = label_parsed)+
    theme_bw()+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium"))) +
    theme(axis.title= element_text(size=18),
          axis.text=element_text(size=16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          strip.text = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "top")+
    scale_x_continuous(labels=c(0.1,2,4,6,8,10),breaks=c(100,2000,4000,6000,8000,10000))+
    xlab(expression(Generation~(10^3)))+
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab("Correlation with Selection")
  
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="DP Pop Size",shape="DP Pop Size")
  }
  
  
  return(list(g1,tmp))
})

pdf(paste0(fig_out,"/Spatial_Model_DivergenceSelectionCoef_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height=10)
print(ggarrange(plot_list[[1]][[1]],
                plot_list[[2]][[1]],
                plot_list[[3]][[1]],ncol = 3,labels = "AUTO"))
dev.off()

# We also want to save the correlation coefficients

# Bottlenecks
write.table(plot_list[[1]][[2]],
            paste0(output_out,"/Spatial_Model_DivergenceSelectionCoef_Bottleneck.txt"),
            row.names = F,quote = F,sep="\t")
# DP Pop Sizes
write.table(plot_list[[2]][[2]],
            paste0(output_out,"/Spatial_Model_DivergenceSelectionCoef_DP Pop Sizes.txt"),
            row.names = F,quote = F,sep="\t")
# Migration
write.table(plot_list[[3]][[2]],
            paste0(output_out,"/Spatial_Model_DivergenceSelectionCoef_Migration.txt"),
            row.names = F,quote = F,sep="\t")

# Plot selection coefficients for all 16 treatments in 1-1 comparison to aid visibility

# Reformat cor_dd so its by generation
# Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column

###################################################################################
# This gets run over each iteration separately and averaged over

iteration_list<-lapply(1:20,function(iter){
  
  cor_dd2<-cbind(slim_dd[slim_dd$Iteration==iter,c("Treatment_Demo","Gen_N","FST","Selection_Coef","dXY","deltaPI")])
  cor_dd2<-as.data.frame(cor_dd2 %>% 
                           group_by(Gen_N,Treatment_Demo) %>%
                           dplyr::summarise((cor.test(FST,Selection_Coef)$estimate),
                                            (cor.test(dXY,Selection_Coef)$estimate),
                                            (cor.test(deltaPI,Selection_Coef)$estimate)))
  colnames(cor_dd2)<-c("Generation","Treatment","FST.cor","Dxy.cor","deltaPI.cor")
  
  # Subset for 4 timepoints
  cor_dd2<-cor_dd2[cor_dd2$Generation %in% c(100,500,3000,10000),]
  
  # Reformat
  cor_dd3<-data.frame(Generation=rep(cor_dd2$Generation,3),
                      Treatment=rep(cor_dd2$Treatment,3),
                      Correlation1=c(cor_dd2$FST.cor,cor_dd2$FST.cor,cor_dd2$Dxy.cor),
                      Correlation2=c(cor_dd2$Dxy.cor,cor_dd2$deltaPI.cor,cor_dd2$deltaPI.cor),
                      Mig=rep(cor_dd$Mig,12),
                      Pop2=rep(cor_dd$Pop2,12),
                      Bot=rep(cor_dd$Bot,12),
                      Measure=rep(c("F[ST]~x~D[XY]","F[ST]~x~Delta*pi","D[XY]~x~Delta*pi"),each=64))
  cor_dd3$Generation <- as.character(cor_dd3$Generation)
  
  
  # Plot
  cor_dd3$Measure<-factor(cor_dd3$Measure,levels=c("F[ST]~x~D[XY]","F[ST]~x~Delta*pi","D[XY]~x~Delta*pi"))
  
  # Plot with colour DP Pop Size and migration shape
  cor_dd3$Mig<-as.factor(cor_dd3$Mig)
  cor_dd3$Pop2<-as.factor(cor_dd3$Pop2)
  
  # Re-order facets
  cor_dd3$Generation<-factor(cor_dd3$Generation,levels=c(100,500,3000,10000))
  
  return(cor_dd3)
})

# Average over the list
cor_dd3<-Reduce("+", iteration_list) / length(iteration_list)
cor_dd3[,c(1:2,5:8)]<-iteration_list[[1]][,c(1:2,5:8)]

###############################################################

g1<-ggplot(data=cor_dd3,aes(x=Correlation1,y=Correlation2,colour=Pop2,shape=Mig))+
  geom_point(size=4)+
  geom_abline(slope = 1,intercept = 0)+
  facet_grid(Measure~Generation,labeller = label_parsed,scales = "free")+
  theme_bw()+
  theme(axis.title= element_text(size=18),
        axis.text=element_text(size=14),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        strip.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "top")+
  labs(colour="DP Pop Size",x="Correlation 1",y="Correlation 2",shape="Migration")+
  scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium")))

pdf(paste0(fig_out,"/Spatial_Model_Comparing_Correlations_with_Selection.pdf"),width=13,height=10)
print(g1)
dev.off()

#-------------------------------------------
# CORRELATIONS WITH PI
#-------------------------------------------


# Figures to examine effects of demographic variables across generations on correlation coefficients with selection
treatment_vec<-c("Bottleneck","Pop2_Size","Migration")

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  ###################################################################################
  # This gets run over each iteration separately and averaged over
  
  iteration_list<-lapply(1:20,function(iter){
    
    # Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column
    # Summarise FST
    tmp1<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","FST",treatment_vec[x],"Treatment_Demo")],Measure=rep("FST"))
    tmp1$Pop1_PI<-rep(slim_dd[slim_dd$Iteration==iter & slim_dd$Gen_N == 100,"Pop1_PI"],each=12)
    
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(FST,Pop1_PI)$estimate),
                                           (cor.test(FST,Pop1_PI)$conf.int[1]),
                                           (cor.test(FST,Pop1_PI)$conf.int[2])))
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    tmp1$Metric<-"FST"
    
    # Summarise 
    tmp2<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","dXY",treatment_vec[x],"Treatment_Demo")],Measure=rep("dXY"))
    tmp2$Pop1_PI<-rep(slim_dd[slim_dd$Iteration==iter & slim_dd$Gen_N == 100,"Pop1_PI"],each=12)
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(dXY,Pop1_PI)$estimate),
                                           (cor.test(dXY,Pop1_PI)$conf.int[1]),
                                           (cor.test(dXY,Pop1_PI)$conf.int[2])))
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    tmp2$Metric<-"dXY"
    
    # Summarise 
    tmp3<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","deltaPI",treatment_vec[x],"Treatment_Demo")],Measure=rep("deltaPI"))
    tmp3$Pop1_PI<-rep(slim_dd[slim_dd$Iteration==iter & slim_dd$Gen_N == 100,"Pop1_PI"],each=12)
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(deltaPI,Pop1_PI)$estimate),
                                           (cor.test(deltaPI,Pop1_PI)$conf.int[1]),
                                           (cor.test(deltaPI,Pop1_PI)$conf.int[2])))
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    tmp3$Metric<-"deltaPI"
    
    # Rbind
    tmp<-rbind(tmp1,tmp2,tmp3)
    tmp[,2]<-as.factor(tmp[,2])
    
    # Plot
    tmp$Metric<-factor(tmp$Metric,levels=c("FST","dXY","deltaPI"))
    levels(tmp$Metric)<-c("F[ST]","D[XY]","Delta*pi")
    
    return(tmp)
  })
  
  # Average over the list
  tmp<-Reduce("+", iteration_list) / length(iteration_list)
  tmp[,2]<-iteration_list[[1]][,2]
  tmp[,6]<-iteration_list[[1]][,6]
  
  ###############################################################
  
  pd <- position_dodge(0.5)
  
  g1<-ggplot(tmp,aes(x=Gen_N,y=Cor.Coef,shape=tmp[,2],colour=tmp[,2]))+
    # geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, position=pd,show.legend = FALSE) +
    geom_point(position=pd,size=3)+
    geom_line()+
    facet_wrap(~Metric,ncol=1,scales = "free",labeller = label_parsed)+
    theme_bw()+
    theme(axis.title= element_text(size=18),
          axis.text=element_text(size=16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          strip.text = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "top")+
    scale_x_continuous(labels=c(0.1,2,4,6,8,10),breaks=c(100,2000,4000,6000,8000,10000))+
    xlab(expression(Generation~(10^3)))+
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab(expression("Correlation with Founding "~pi))+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium")))
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="DP Pop Size",shape="DP Pop Size")
  }
  
  return(list(g1,tmp))
})

pdf(paste0(fig_out,"/Spatial_Model_DivergencePiCoef_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height=10)
print(ggarrange(plot_list[[1]][[1]],
                plot_list[[2]][[1]],
                plot_list[[3]][[1]],ncol = 3,labels = "AUTO"))
dev.off()

# We also want to save the correlation coefficients

# Bottlenecks
write.table(plot_list[[1]][[2]],
            paste0(output_out,"/Spatial_Model_DivergencePiCoef_Bottleneck.txt"),
            row.names = F,quote = F,sep="\t")
# DP Pop Sizes
write.table(plot_list[[2]][[2]],
            paste0(output_out,"/Spatial_Model_DivergencePiCoef_DP Pop Sizes.txt"),
            row.names = F,quote = F,sep="\t")
# Migration
write.table(plot_list[[3]][[2]],
            paste0(output_out,"/Spatial_Model_DivergencePiCoef_Migration.txt"),
            row.names = F,quote = F,sep="\t")

#------------------------------
# Correlation of measures
#------------------------------

# Figures to examine effects of demographic variables across generations on correlation coefficients between measures
treatment_vec<-c("Bottleneck","Pop2_Size","Migration")

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  ###################################################################################
  # This gets run over each iteration separately and averaged over
  
  iteration_list<-lapply(1:20,function(iter){
    
    # Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column
    # Summarise FST
    tmp1<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","FST","dXY",treatment_vec[x],"Treatment_Demo")],Measure=rep("FST x DXY"))
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(FST,dXY)$estimate),
                                           (cor.test(FST,dXY)$conf.int[1]),
                                           (cor.test(FST,dXY)$conf.int[2])))
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    
    
    # Summarise 
    tmp2<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","FST","deltaPI",treatment_vec[x],"Treatment_Demo")],Measure=rep("FST x deltaPI"))
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(FST,deltaPI)$estimate),
                                           (cor.test(FST,deltaPI)$conf.int[1]),
                                           (cor.test(FST,deltaPI)$conf.int[2])))
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    
    # Summarise 
    tmp3<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","dXY","deltaPI",treatment_vec[x],"Treatment_Demo")],Measure=rep("dXY x deltaPI"))
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(dXY,deltaPI)$estimate),
                                           (cor.test(dXY,deltaPI)$conf.int[1]),
                                           (cor.test(dXY,deltaPI)$conf.int[2])))
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    
    
    # Rbind
    tmp<-rbind(tmp1,tmp2,tmp3)
    
    return(tmp)
  })
  
  # Average over the list
  tmp<-Reduce("+", iteration_list) / length(iteration_list)
  tmp[,2]<-as.factor(iteration_list[[1]][,2])
  tmp$Metric<-as.factor(rep(c("FST x DXY","FST x deltaPI","dXY x deltaPI"),each=nrow(tmp)/3))
  
  ###############################################################
  
  
  # Plot
  pd <- position_dodge(0.5)
  
  tmp$Metric<-factor(tmp$Metric,levels=c("FST x DXY","FST x deltaPI","dXY x deltaPI"))
  levels(tmp$Metric)<-c("F[ST]~x~D[XY]","F[ST]~x~Delta*pi","D[XY]~x~Delta*pi")
  
  
  g1<-ggplot(tmp,aes(x=Gen_N,y=Cor.Coef,shape=tmp[,2],colour=tmp[,2]))+
    #geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, position=pd,show.legend = FALSE) +
    geom_point(position=pd,size=3)+
    geom_line()+
    facet_wrap(~Metric,ncol=1,scales = "free",labeller = label_parsed)+
    theme_bw()+
    theme(axis.title= element_text(size=18),
          axis.text=element_text(size=16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          strip.text = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "top")+
    scale_x_continuous(labels=c(0.1,2,4,6,8,10),breaks=c(100,2000,4000,6000,8000,10000))+
    xlab(expression(Generation~(10^3)))+
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab("Correlation between measures")+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium")))
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="DP Pop Size",shape="DP Pop Size")
  }
  
  return(list(g1,tmp))
})

pdf(paste0(fig_out,"/Spatial_model_DivergenceCoefs_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height = 10)
print(ggarrange(plot_list[[1]][[1]],
                plot_list[[2]][[1]],
                plot_list[[3]][[1]],ncol = 3,labels = "AUTO"))
dev.off()

# We also want to save the correlation coefficients

# Bottlenecks
write.table(plot_list[[1]][[2]],
            paste0(output_out,"/Spatial_Model_DivergenceCoefs_Bottleneck.txt"),
            row.names = F,quote = F,sep="\t")
# DP Pop Sizes
write.table(plot_list[[2]][[2]],
            paste0(output_out,"/Spatial_Model_DivergenceCoefs_DP Pop Sizes.txt"),
            row.names = F,quote = F,sep="\t")
# Migration
write.table(plot_list[[3]][[2]],
            paste0(output_out,"/Spatial_Model_DivergenceCoefs_Migration.txt"),
            row.names = F,quote = F,sep="\t")

#-----------------------------------
# Correlations with dynamic PI in diverging pop
#-----------------------------------

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  ###################################################################################
  # This gets run over each iteration separately and averaged over
  
  iteration_list<-lapply(1:20,function(iter){
    
    
    # Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column
    # Summarise FST
    tmp1<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","FST","Pop2_PI",treatment_vec[x],"Treatment_Demo")],Measure=rep("FST"))
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(FST,Pop2_PI)$estimate),
                                           (cor.test(FST,Pop2_PI)$conf.int[1]),
                                           (cor.test(FST,Pop2_PI)$conf.int[2])))
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp1<-as.data.frame(tmp1 %>% 
                          group_by(Gen_N,tmp1[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp1)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    tmp1$Metric<-"FST"
    
    # Summarise 
    tmp2<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","dXY","Pop2_PI",treatment_vec[x],"Treatment_Demo")],Measure=rep("dXY"))
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(dXY,Pop2_PI)$estimate),
                                           (cor.test(dXY,Pop2_PI)$conf.int[1]),
                                           (cor.test(dXY,Pop2_PI)$conf.int[2])))
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp2<-as.data.frame(tmp2 %>% 
                          group_by(Gen_N,tmp2[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp2)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    tmp2$Metric<-"dXY"
    
    # Summarise 
    tmp3<-cbind(slim_dd[slim_dd$Iteration==iter,c("Gen_N","deltaPI","Pop2_PI",treatment_vec[x],"Treatment_Demo")],Measure=rep("deltaPI"))
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]],Treatment_Demo) %>%
                          dplyr::summarise((cor.test(deltaPI,Pop2_PI)$estimate),
                                           (cor.test(deltaPI,Pop2_PI)$conf.int[1]),
                                           (cor.test(deltaPI,Pop2_PI)$conf.int[2])))
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Treatment_Demo","Cor.Coef","LL","UL")
    
    # Average again...
    tmp3<-as.data.frame(tmp3 %>% 
                          group_by(Gen_N,tmp3[,treatment_vec[x]]) %>%
                          dplyr::summarise(Cor.Coef=mean(Cor.Coef),
                                           LL=mean(LL),
                                           UL=mean(UL)))
    
    colnames(tmp3)<-c("Gen_N",treatment_vec[x],"Cor.Coef","LL","UL")
    tmp3$Metric<-"deltaPI"
    
    
    # Rbind
    tmp<-rbind(tmp1,tmp2,tmp3)
    
    return(tmp)
  })
  
  # Average over the list
  tmp<-Reduce("+", iteration_list) / length(iteration_list)
  tmp[,2]<-as.factor(iteration_list[[1]][,2])
  tmp$Metric<-rep(c("FST","dXY","deltaPI"),each=nrow(tmp)/3)
  
  ###############################################################
  
  
  # Plot
  pd <- position_dodge(0.5)
  tmp$Metric<-factor(tmp$Metric,levels=c("FST","dXY","deltaPI"))
  levels(tmp$Metric)<-c("F[ST]","D[XY]","Delta*pi")
  
  
  g1<-ggplot(tmp,aes(x=Gen_N,y=Cor.Coef,shape=tmp[,2],colour=tmp[,2]))+
    #geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, position=pd,show.legend = FALSE) +
    geom_point(position=pd,size=3)+
    geom_line()+
    facet_wrap(~Metric,ncol=1,scales = "free",labeller = label_parsed)+
    theme_bw()+
    theme(axis.title= element_text(size=18),
          axis.text=element_text(size=16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          strip.text = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "top")+
    scale_x_continuous(labels=c(0.1,2,4,6,8,10),breaks=c(100,2000,4000,6000,8000,10000))+
    xlab(expression(Generation~(10^3)))+
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab(expression("Correlation with Diverging "~pi))+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium")))
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="DP Pop Size",shape="DP Pop Size")
  }
  
  
  return(list(g1,tmp))
})

pdf(paste0(fig_out,"/Spatial_Model_DivergencePi2Coef_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height=10)
print(ggarrange(plot_list[[1]][[1]],
                plot_list[[2]][[1]],
                plot_list[[3]][[1]],ncol = 3,labels = "AUTO"))
dev.off()

# We also want to save the correlation coefficients

# Bottlenecks
write.table(plot_list[[1]][[2]],
            paste0(output_out,"/Spatial_Model_DivergencePi2Coef_Bottleneck.txt"),
            row.names = F,quote = F,sep="\t")
# DP Pop Sizes
write.table(plot_list[[2]][[2]],
            paste0(output_out,"/Spatial_Model_DivergencePi2Coef_DP Pop Sizes.txt"),
            row.names = F,quote = F,sep="\t")
# Migration
write.table(plot_list[[3]][[2]],
            paste0(output_out,"/Spatial_Model_DivergencePi2Coef_Migration.txt"),
            row.names = F,quote = F,sep="\t")


##### Correlation Matrices #####
# Run over all 16 treatments
treatments_to_keep<-unique(slim_dd$Treatment_Demo)
generations_to_plot<-c(100,500,3000,10000)

###################################################################################
# This gets permuted over different combinations of iterations

iteration_list<-mclapply(1:100,function(iter){
  
  # Draw the permutation iterations
  perm_iters<-sample(1:20,16)
  
  comps_to_make<-combn(1:16,2)
  
  # Get the permutation iterations
  iters_to_perm<-comps_to_make
  for(i in 1:ncol(iters_to_perm)){
    iters_to_perm[1,i]<-perm_iters[comps_to_make[1,i]]
    iters_to_perm[2,i]<-perm_iters[comps_to_make[2,i]]
  }
  
  #for (level_gen in generations_to_plot) {
  per_gen_list<-lapply(generations_to_plot,function(level_gen){
    slim_dd_GEN4<-slim_dd[slim_dd$Gen_N == level_gen,!(colnames(slim_dd) %in% c("Exon_Length_Total"))]
    #slim_dd_GEN4<-na.omit(slim_dd_GEN4)
    
    FST_cov<-vector()
    for(i in 1:length(treatments_to_keep)) {
      tmp<-as.vector(slim_dd_GEN4[slim_dd_GEN4$Iteration == perm_iters[i] & slim_dd_GEN4$Treatment_Demo == treatments_to_keep[i],"FST"])
      FST_cov<-cbind(FST_cov,tmp)
    }
    colnames(FST_cov)<-treatments_to_keep
    FST_cor<-cor(FST_cov,method = "spearman")


    dXY_cov<-vector()
    for(i in 1:length(treatments_to_keep)) {
      tmp<-as.vector(slim_dd_GEN4[slim_dd_GEN4$Iteration == perm_iters[i] & slim_dd_GEN4$Treatment_Demo == treatments_to_keep[i],"dXY"])
      dXY_cov<-cbind(dXY_cov,tmp)
    }
    colnames(dXY_cov)<-treatments_to_keep
    dxy_cor<-cor(dXY_cov,method = "spearman")
    
    
    deltaPI_cov<-vector()
    for(i in 1:length(treatments_to_keep)) {
      tmp<-as.vector(slim_dd_GEN4[slim_dd_GEN4$Iteration == perm_iters[i] & slim_dd_GEN4$Treatment_Demo == treatments_to_keep[i],"deltaPI"])
      deltaPI_cov<-cbind(deltaPI_cov,tmp)
    }
    colnames(deltaPI_cov)<-treatments_to_keep
    deltaPI_cor<-cor(deltaPI_cov,method = "spearman")
    
    return(list(FST_cor,dxy_cor,deltaPI_cor))  
  })
  
  # Return the iteration list
  return(per_gen_list)
},mc.cores=3)

# Need to average Per Gen and plot, Per Measure
lapply(1:length(generations_to_plot),function(gen_X){
  
  # Set this
  level_gen<-generations_to_plot[gen_X]
  
  # Retrieve all matrices for each
  FST_list<-lapply(iteration_list,function(list_x){return(list_x[[gen_X]][[1]])})
  dXY_list<-lapply(iteration_list,function(list_x){return(list_x[[gen_X]][[2]])})
  deltaPI_list<-lapply(iteration_list,function(list_x){return(list_x[[gen_X]][[3]])})
  
  # Average each
  FST_cor<-Reduce("+", FST_list) / length(FST_list)
  dxy_cor<-Reduce("+", dXY_list) / length(dXY_list)
  deltaPI_cor<-Reduce("+", deltaPI_list) / length(deltaPI_list)
  
  # Plot matrices
  pdf(paste0(fig_out,"/Spatial_model_",level_gen,"_FST_dXY_deltaPI_Correlation_matrices.pdf"),width=10,height=10)
  print(plot_cor_mat(as.matrix(FST_cor),cor = F)+
          theme(legend.text = element_text(angle=45))+
          scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                               midpoint = 0, limit = c(-1,1), space = "Lab", 
                               name=expression(Spearman~rho))) 
  print(plot_cor_mat(as.matrix(dxy_cor),cor = F)+
          theme(legend.text = element_text(angle=45))+
          scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                               midpoint = 0, limit = c(-1,1), space = "Lab", 
                               name=expression(Spearman~rho))) 
  print(plot_cor_mat(as.matrix(deltaPI_cor),cor = F)+
          theme(legend.text = element_text(angle=45))+
          scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                               midpoint = 0, limit = c(-1,1), space = "Lab", 
                               name=expression(Spearman~rho))) 
  dev.off()
})

##### False-positive rate across time
# This is new for the revised manuscript, compares the distribution of True results to Null results across time and across treatments...
library(ggridges)
FPR_data<-rbind(slim_dd,slim_null,slim_neutral)

# We will plot over the previously-used plotting generations
FPR_analysis<-lapply(generations_to_plot,function(gen){
  
  # Subset for generation
  FPR_sub<-FPR_data[FPR_data$Gen_N==gen,]
  
  # For each treatment group, get the 0.95 quantile of each measure
  quantile_list<-data.frame(rbindlist(lapply(unique(FPR_sub$run_type),function(run){
    tmp<-FPR_sub[FPR_sub$run_type == run,]
    quantiles_to_plot<-matrix(ncol = 4,nrow = 16)
    quantiles_to_plot[,1]<-as.character(unique(tmp$Treatment_Demo))
    for(i in 1:nrow(quantiles_to_plot)){
      quantiles_to_plot[i,2]<-quantile(tmp[tmp$Treatment_Demo == quantiles_to_plot[i,1],"FST"],0.95)
      quantiles_to_plot[i,3]<-quantile(tmp[tmp$Treatment_Demo == quantiles_to_plot[i,1],"dXY"],0.95)
      quantiles_to_plot[i,4]<-quantile(tmp[tmp$Treatment_Demo == quantiles_to_plot[i,1],"deltaPI"],0.95)
    }
    
    out<-as.data.frame(quantiles_to_plot)
    colnames(out)<-c("Treatment_Demo","FST","dXY","deltaPI")
    out$run_type<-run
    return(out)
  })))
  
  # Change y-axes
  #FPR_sub[FPR_sub$run_type=="Neutral","run_type"]<-"Neutrak"
  #quantile_list[quantile_list$run_type=="Neutral","run_type"]<-"Divergent"
  
  # Plot for each measure
  measures<-c("FST","dXY","deltaPI")
  labs<-c(expression(F[ST]),expression(D[XY]),expression(Delta*pi))
  
  # Change the order of treatments
  FPR_sub$Treatment_Demo<-factor(FPR_sub$Treatment_Demo,levels=c("Bot=100/Pop2=0.01/Mig=0","Bot=100/Pop2=0.1/Mig=0","Bot=100/Pop2=0.01/Mig=0.002","Bot=100/Pop2=0.1/Mig=0.002",
                                                                 "Bot=100/Pop2=0.5/Mig=0","Bot=100/Pop2=1/Mig=0","Bot=100/Pop2=0.5/Mig=0.002","Bot=100/Pop2=1/Mig=0.002",
                                                                 "Bot=1000/Pop2=0.01/Mig=0","Bot=1000/Pop2=0.1/Mig=0","Bot=1000/Pop2=0.01/Mig=0.002","Bot=1000/Pop2=0.1/Mig=0.002",
                                                                 "Bot=1000/Pop2=0.5/Mig=0","Bot=1000/Pop2=1/Mig=0","Bot=1000/Pop2=0.5/Mig=0.002","Bot=1000/Pop2=1/Mig=0.002")
  )
  

  # We are also interested in aspects of the distributions
  # What proportion of the divergent selected genes lie above the neutral 0.95 cut-off, ie. false-negative
  neutral_quantiles<-quantile_list[quantile_list$run_type=="Neutral",]
  false_negatives<-data.frame(Treatment=neutral_quantiles$Treatment_Demo,
                              FNR_FST=sapply(1:nrow(neutral_quantiles),function(x){
                                tmp<-FPR_sub[as.character(FPR_sub$Treatment_Demo) == as.character(neutral_quantiles$Treatment_Demo[x]) & FPR_sub$run_type=="Divergent",]
                                return(nrow(tmp[tmp$FST < as.numeric(as.character(neutral_quantiles$FST[x])),])/nrow(tmp))}),
                              FNR_dXY=sapply(1:nrow(neutral_quantiles),function(x){
                                tmp<-FPR_sub[as.character(FPR_sub$Treatment_Demo) == as.character(neutral_quantiles$Treatment_Demo[x]) & FPR_sub$run_type=="Divergent",]
                                return(nrow(tmp[tmp$dXY < as.numeric(as.character(neutral_quantiles$dXY[x])),])/nrow(tmp))}),
                              FNR_deltaPI=sapply(1:nrow(neutral_quantiles),function(x){
                                tmp<-FPR_sub[as.character(FPR_sub$Treatment_Demo) == as.character(neutral_quantiles$Treatment_Demo[x]) & FPR_sub$run_type=="Divergent",]
                                return(nrow(tmp[tmp$deltaPI < as.numeric(as.character(neutral_quantiles$deltaPI[x])),])/nrow(tmp))}))
  
  false_negatives$Generation<-gen
  
  # What proportion of the upper 1% quantile are truly selected genes
  # We can permute this...
  FPR_perms<-mclapply(1:100,function(perm){
    set.seed(perm)
    
    downsampled_divergents<-sample(1:20,100,replace=T)
    
    # Run separately for each treatment group
    treats<-lapply(1:nrow(neutral_quantiles),function(y){
      
      subset_N<-nrow(FPR_sub[as.character(FPR_sub$Treatment_Demo) == as.character(neutral_quantiles$Treatment_Demo[y]) & 
                               FPR_sub$run_type=="Neutral",])
      
      # Subset for the right treatment group and replace 
      neutral_tmp<-FPR_sub[as.character(FPR_sub$Treatment_Demo) == as.character(neutral_quantiles$Treatment_Demo[y]) & 
                             FPR_sub$run_type=="Neutral",]
      neutral_tmp$index<-paste0(neutral_tmp$Gene_ID,"-",neutral_tmp$Iteration,"-",neutral_tmp$run_type)
      neutral_tmp<-neutral_tmp[neutral_tmp$Iteration != sample(1:20,1),]
      divergent_tmp<-FPR_sub[as.character(FPR_sub$Treatment_Demo) == as.character(neutral_quantiles$Treatment_Demo[y]) & FPR_sub$run_type=="Neutral",]
      
      # Get the right iteration for each gene
      divergent_tmp<-data.frame(rbindlist(lapply(1:100,function(gene){return(divergent_tmp[divergent_tmp$Gene_ID == gene & divergent_tmp$Iteration == downsampled_divergents[gene],])})))
      divergent_tmp$index<-paste0(divergent_tmp$Gene_ID,"-",divergent_tmp$run_type)
      # Combine and index
      perm_tmp<-rbind(neutral_tmp,divergent_tmp)
      
      # Get upper 5%
      fst_top<-perm_tmp[perm_tmp$FST > quantile(perm_tmp$FST,0.95),"index"]
      dxy_top<-perm_tmp[perm_tmp$dXY > quantile(perm_tmp$dXY,0.95),"index"]
      deltaPI_top<-perm_tmp[perm_tmp$deltaPI > quantile(perm_tmp$deltaPI,0.95),"index"]
      
      # Get combined upper 5%'s
      fst_dxy_top<-perm_tmp[perm_tmp$FST > quantile(perm_tmp$FST,0.95) &
                              perm_tmp$dXY > quantile(perm_tmp$dXY,0.95),"index"]
      fst_deltaPI_top<-perm_tmp[perm_tmp$FST > quantile(perm_tmp$FST,0.95) &
                                  perm_tmp$deltaPI > quantile(perm_tmp$deltaPI,0.95),"index"]
      dxy_deltaPI_top<-perm_tmp[perm_tmp$dXY > quantile(perm_tmp$dXY,0.95) &
                                  perm_tmp$deltaPI > quantile(perm_tmp$deltaPI,0.95),"index"]
      all_top<-Reduce(intersect,list(fst_dxy_top,fst_deltaPI_top,dxy_deltaPI_top))
      
      # Sample the upper quantile of each measure, bearing in mind that divergent indices are always > 2000
      false_pos<-data.frame(Treatment=neutral_quantiles$Treatment_Demo[y],
                            FPR_FST=length(fst_top[grep("Neutral",fst_top)])/length(fst_top),
                            FPR_dXY=length(dxy_top[grep("Neutral",dxy_top)])/length(dxy_top),
                            FPR_deltaPI=length(deltaPI_top[grep("Neutral",deltaPI_top)])/length(deltaPI_top),
                            
                            Outliers_FST_dXY=length(fst_dxy_top),
                            FPR_FST_dXY=length(fst_dxy_top[grep("Neutral",fst_dxy_top)])/length(fst_dxy_top),
                            FNR_FST_dXY=(nrow(perm_tmp[perm_tmp$run_type=="Divergent",])-length(fst_dxy_top[grep("Divergent",fst_dxy_top)]))/nrow(perm_tmp[perm_tmp$run_type=="Divergent",]),
                            
                            Outliers_FST_deltaPI=length(fst_deltaPI_top),
                            FPR_FST_deltaPI=length(fst_deltaPI_top[grep("Neutral",fst_deltaPI_top)])/length(fst_deltaPI_top),
                            FNR_FST_deltaPI=(nrow(perm_tmp[perm_tmp$run_type=="Divergent",])-length(fst_deltaPI_top[grep("Divergent",fst_deltaPI_top)]))/nrow(perm_tmp[perm_tmp$run_type=="Divergent",]),
                            
                            Outliers_dXY_deltaPI=length(dxy_deltaPI_top),
                            FPR_dXY_deltaPI=length(dxy_deltaPI_top[grep("Neutral",dxy_deltaPI_top)])/length(dxy_deltaPI_top),
                            FNR_dXY_deltaPI=(nrow(perm_tmp[perm_tmp$run_type=="Divergent",])-length(dxy_deltaPI_top[grep("Divergent",dxy_deltaPI_top)]))/nrow(perm_tmp[perm_tmp$run_type=="Divergent",]),
                            
                            Outliers_all=length(all_top),
                            FPR_all=length(all_top[grep("Neutral",all_top)])/length(all_top),
                            FNR_all=(nrow(perm_tmp[perm_tmp$run_type=="Divergent",])-length(all_top[grep("Divergent",all_top)]))/nrow(perm_tmp[perm_tmp$run_type=="Divergent",]))
      
      
      # Remove empty FPR, assuming empty means no overlapping outliers at all, so false positive of 0
      false_pos[is.na(false_pos)]<-0
      
      # Also return a list with all outliers
      outlier_list<-list(fst_top,dxy_top,deltaPI_top,fst_dxy_top,fst_deltaPI_top,dxy_deltaPI_top,all_top)
      
      return(list(false_pos,outlier_list))                       
    })
    
    # Retrieve treats_out
    treats_out<-data.frame(rbindlist(lapply(1:length(treats),function(treat){return(treats[[treat]][[1]])})))
    
    # Retrieve outlier lists
    outlier_lists<-lapply(1:length(treats),function(treat){return(treats[[treat]][[2]])})
    
    # End the permutations
    return(list(treats_out,outlier_lists))  
  },mc.cores=3)
  
  # Retrieve FPR/FNR
  FPR_perms_out<-lapply(1:length(FPR_perms),function(perm){return(FPR_perms[[perm]][[1]])})
  perm_outlier_lists<-lapply(1:length(FPR_perms),function(perm){return(FPR_perms[[perm]][[2]])})
  
  # Average the results
  perm_out<-Reduce("+", FPR_perms_out) / length(FPR_perms_out)
  perm_out[,1]<-FPR_perms_out[[1]][,1]
  
  # Set the generation
  perm_out$Generation<-gen
  
  # Re-order the tables
  treatment_order<-c("Bot=100/Pop2=0.01/Mig=0","Bot=1000/Pop2=0.01/Mig=0","Bot=100/Pop2=0.1/Mig=0","Bot=1000/Pop2=0.1/Mig=0","Bot=100/Pop2=0.5/Mig=0","Bot=1000/Pop2=0.5/Mig=0","Bot=100/Pop2=1/Mig=0","Bot=1000/Pop2=1/Mig=0","Bot=1000/Pop2=0.01/Mig=0.002","Bot=100/Pop2=0.01/Mig=0.002","Bot=1000/Pop2=0.1/Mig=0.002","Bot=100/Pop2=0.1/Mig=0.002","Bot=100/Pop2=0.5/Mig=0.002","Bot=1000/Pop2=0.5/Mig=0.002","Bot=100/Pop2=1/Mig=0.002","Bot=1000/Pop2=1/Mig=0.002")
  perm_out<-perm_out[match(treatment_order,perm_out$Treatment),]
  
  
  # Before returning, we need to compare outlier_lists across treatments
  avg_heat_list<-mclapply(1:length(perm_outlier_lists),function(perm){
    
    # Subset
    perm_tmp<-perm_outlier_lists[[perm]]
    
    # We have 7 sets of outliers, each needs to compared in all pairwise comparisons of the 16 treatments
    comparisons_to_make<-combn(1:16,2)
    
    # Output to a matrix, 1 per set of outliers
    matrix_overlap<-lapply(1:7,function(set){
      
      # Fill a matrix
      mat_tmp<-matrix(ncol=16,nrow=16)
      for(comp in 1:ncol(comparisons_to_make)){
        row_N<-comparisons_to_make[1,comp]
        col_N<-comparisons_to_make[2,comp]
        overlap<-Reduce(intersect,list(perm_tmp[[row_N]][[set]],perm_tmp[[col_N]][[set]]))
        mat_tmp[row_N,col_N]<-length(overlap)/min(c(length(perm_tmp[[row_N]][[set]]),length(perm_tmp[[col_N]][[set]])))
      }
      
      mat_tmp[is.na(mat_tmp)]<-0
      
      rownames(mat_tmp)<-unique(slim_dd$Treatment_Demo)
      colnames(mat_tmp)<-unique(slim_dd$Treatment_Demo)
      
      return(mat_tmp)
    })
    
    # So now we have the matrices for the 7 outlier sets, across the 16 comparisons for 1 permutation  
    return(matrix_overlap)    
  },mc.cores=3)
  
  # We now need to average over all the permuted heatmaps
  avg_out<-lapply(1:7,function(set){
    tmp_list<-lapply(avg_heat_list,function(heat){return(heat[[set]])})
    tmp_out<-Reduce("+", tmp_list) / length(tmp_list)
    for(i in 1:nrow(tmp_out)){
      tmp_out[i,i]<-1
    }
    return(tmp_out)
  })
  
  # And plot with titles
  heat_titles<-c(expression(F[ST]),expression(D[XY]),expression(Delta~pi),expression(F[ST]~+~D[XY]),expression(F[ST]~+~Delta~pi),expression(D[XY]~+~Delta~pi),expression(All))
  prop_mat_plots<-lapply(1:length(avg_out),function(plot){
    return(plot_prop_matrix(avg_out[[plot]])+
             ggtitle(heat_titles[plot])+
             theme(title = element_text(size=24),
                   legend.title = element_text(size=15))+
             scale_fill_gradient2(low = "white", high = "red2",
                                  limit=c(0,1), space = "Lab", 
                                  name="Proportion\nOverlap"))
  })
  
  # Plot them all separately
  pdf(paste0(fig_out,"/Downsampled_outlier_overlap_separated_gen_",gen,".pdf"),width=8,height=8)
  for(i in 1:length(prop_mat_plots)){
    print(prop_mat_plots[[i]])
  }
  dev.off()
  
  # Plot the first 6 together without heat legends
  prop_mat_plots<-lapply(1:length(avg_out),function(plot){
    return(plot_prop_matrix(avg_out[[plot]])+
             ggtitle(heat_titles[plot])+
             theme(title = element_text(size=24),
                   legend.position = "none")+
             scale_fill_gradient2(low = "white", high = "red2",
                                  limit=c(0,1), space = "Lab", 
                                  name="Proportion\nOverlap"))
  })
  pdf(paste0(fig_out,"/Downsampled_outlier_overlap_together_gen_",gen,".pdf"),width=14,height=18)
  print(ggarrange(prop_mat_plots[[1]],prop_mat_plots[[4]],prop_mat_plots[[2]],prop_mat_plots[[5]],prop_mat_plots[[3]],prop_mat_plots[[6]],
                  ncol=2,nrow=3))
  dev.off()
  
  # Return false negs and false pos
  return(list(false_negatives,perm_out))
})


# Retrieve the outputs and save
for(i in 1:4){
  write.table(FPR_analysis[[i]][[1]],
              paste0(output_out,"FNR_95q_generation_",generations_to_plot[i],".txt"),quote = F,row.names = F,sep="\t")
  write.table(FPR_analysis[[i]][[2]],
              paste0(output_out,"FPR_95q_generation_",generations_to_plot[i],".txt"),quote = F,row.names = F,sep="\t")
}


##### What explains measures? #####

###########################################################################
# Explain FST at final stage

# Subset for final generation
slim_dd_GEN4<-slim_dd[slim_dd$Gen_N==10000,]

######################################################
# FST first for treatments WITH migration
lmm1<-lmer(FST~Selection_Coef+log(exon_N)+log(Target_sel)+log(gene_length)+log(Evolv_Gen)+(1|Gene_ID)+(1|Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0.002 &
                                                                                                                                        slim_dd_GEN4$Evolv_Gen > 0,])
drop1(lmm1)
model.step1<-step(lmm1,scope=~.^2,direction="both")
lmm_fixed<-lmer(FST ~ Selection_Coef +log(exon_N) + log(Target_sel) + log(Evolv_Gen) + (1 | Gene_ID) + (1 | Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0.002&
                                                                                                                                           slim_dd_GEN4$Evolv_Gen > 0,])
# Evaluate
drop1(lmm_fixed)

# plot
plot(allEffects(lmm_fixed))
summary(lmm_fixed)
# Selection has strongest effect

# Plot
Fst_selection_plot<-slim_dd_GEN4 %>% group_by(Gene_ID,Migration) %>% summarise_all(.funs=mean)
Fst_selection_plot2<-slim_dd_GEN4 %>% group_by(Gene_ID,Migration,Treatment_Demo) %>% summarise_all(.funs=mean)
Fst_selection_fig<-ggplot(Fst_selection_plot,aes(Selection_Coef,FST))+
  geom_point(size=4,shape=17)+
  geom_smooth(method="lm")+
  geom_point(data=Fst_selection_plot2,alpha=0.2)+
  facet_wrap(~Migration)+
  theme_bw()+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        strip.text = element_text(size=18))+
  labs(y=expression(F[ST]),x=expression(Selection~Coefficient~(log[10])))

######################################################
# FST first for treatments WITHOUT migration
lmm2<-lmer(FST~Selection_Coef+log(exon_N)+log(Target_sel)+log(gene_length)+log(Evolv_Gen)+(1|Gene_ID)+(1|Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0&
                                                                                                                                        slim_dd_GEN4$Evolv_Gen > 0,])
drop1(lmm2)
model.step2<-step(lmm2,scope=~.^2,direction="both")
lmm2_fixed<-lmer(FST ~ Selection_Coef + log(Target_sel) + (1 | Gene_ID) + (1 | Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0 &
                                                                                                              slim_dd_GEN4$Evolv_Gen > 0,])
# Evaluate
drop1(lmm2_fixed)

# plot
plot(allEffects(lmm2_fixed))
summary(lmm2_fixed)
# Selection has strongest effect, alongside positive relationships of exon N and target selection

# Plot
Fst_target_plot<-slim_dd_GEN4 %>% group_by(Gene_ID,Migration) %>% summarise_all(.funs=mean)
Fst_target_plot2<-slim_dd_GEN4 %>% group_by(Gene_ID,Migration,Treatment_Demo) %>% summarise_all(.funs=mean)

Fst_target_fig<-ggplot(Fst_target_plot,aes(log(Target_sel*100),FST))+
  geom_point(size=4,shape=17)+
  geom_smooth(method="lm")+
  geom_point(data=Fst_target_plot2,alpha=0.2)+
  facet_wrap(~Migration)+
  theme_bw()+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        strip.text = element_text(size=18))+
  labs(y=expression(F[ST]),x="Selection Target (log[%])")

# Plot these together
FST_plots<-ggarrange(Fst_selection_fig,Fst_target_fig,ncol=1,nrow=2)

pdf(paste0(fig_out,"/Fst_effects_demography.pdf"),width=12,height=12)
print(FST_plots)
dev.off()

######################################################
# DXY first for treatments WITH migration
lmm3<-lmer(dXY~Selection_Coef+log(exon_N)+log(Target_sel)+log(gene_length)+log(Evolv_Gen)+(1|Gene_ID)+(1|Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0.002 &
                                                                                                                                        slim_dd_GEN4$Evolv_Gen > 0,])

drop1(lmm3)
model.step3<-step(lmm3,scope=~.^2,direction="both")
lmm3_fixed<-lmer(dXY ~ Selection_Coef + log(Target_sel) + log(Evolv_Gen) + (1 | Gene_ID) + (1 | Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0.002 &
                                                                                                                               slim_dd_GEN4$Evolv_Gen > 0,])

# Evaluate
drop1(lmm3_fixed)

# plot
plot(allEffects(lmm3_fixed))
summary(lmm3_fixed)
# Selection has strongest effect but also strong effect of target of selection

######################################################
# DXY now for treatments WITHOUT migration
lmm4<-lmer(dXY~Selection_Coef+log(exon_N)+log(Target_sel)+log(gene_length)+log(Evolv_Gen)+(1|Gene_ID/Iteration)+(1|Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0 &
                                                                                                                                                  slim_dd_GEN4$Evolv_Gen > 0,])
drop1(lmm4)
model.step4<-step(lmm4,scope=~.^2,direction="both")
lmm4_fixed<-lmer(dXY ~ log(Target_sel) + (1 | Gene_ID) + (1 | Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0 &
                                                                                             slim_dd_GEN4$Evolv_Gen > 0,])

# Evaluate
drop1(lmm4_fixed)

# plot
plot(allEffects(lmm4_fixed))
summary(lmm4_fixed)
# Target has strongest effect but also strong effect of target of selection

# Plot
dXY_selection_plot<-slim_dd_GEN4 %>% group_by(Gene_ID,Migration) %>% summarise_all(.funs=mean)
dXY_selection_plot2<-slim_dd_GEN4 %>% group_by(Gene_ID,Migration,Treatment_Demo) %>% summarise_all(.funs=mean)
dXY_selection_fig<-ggplot(dXY_selection_plot,aes(Selection_Coef,dXY))+
  geom_point(size=4,shape=17)+
  geom_smooth(method="lm")+
  geom_point(data=dXY_selection_plot2,alpha=0.2)+
  facet_wrap(~Migration)+
  theme_bw()+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        strip.text = element_text(size=18))+
  labs(y=expression(D[XY]),x=expression(Selection~Coefficient~(log[10])))

# Repeat the dXY figures for DXY
dXY_target_fig<-ggplot(dXY_selection_plot,aes(log(Target_sel*100),dXY))+
  geom_point(size=4,shape=17)+
  geom_smooth(method="lm")+
  geom_point(data=dXY_selection_plot2,alpha=0.2)+
  facet_wrap(~Migration)+
  theme_bw()+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        strip.text = element_text(size=18))+
  labs(y=expression(D[XY]),x="Selection Target (log[%])")

# Plot these together
dXY_plots<-ggarrange(dXY_selection_fig,dXY_target_fig,ncol=1,nrow=2)


######################################################
# deltaPI first for treatments WITH migration
lmm5<-lmer(deltaPI~Selection_Coef+log(exon_N)+log(Target_sel)+log(gene_length)+log(Evolv_Gen)+(1|Gene_ID)+(1|Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0.002 &
                                                                                                                                            slim_dd_GEN4$Evolv_Gen > 0,])

drop1(lmm5)
model.step5<-step(lmm5,scope=~.^2,direction="both")
lmm5_fixed<-lmer(deltaPI ~ Selection_Coef + log(Target_sel) + (1 | Gene_ID) + (1 | Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0.002 & slim_dd_GEN4$Evolv_Gen > 0,])

# Evaluate
drop1(lmm5_fixed)

# plot
plot(allEffects(lmm5_fixed))
summary(lmm5_fixed)
# Selection has strongest effect but also strong effect of target of selection

######################################################
# deltaPI now for treatments WITHOUT migration
lmm6<-lmer(deltaPI~Selection_Coef+log(exon_N)+log(Target_sel)+log(gene_length)+log(Evolv_Gen)+(1|Gene_ID)+(1|Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0 &
                                                                                                                                            slim_dd_GEN4$Evolv_Gen > 0,])
drop1(lmm6)
model.step6<-step(lmm6,scope=~.^2,direction="both")
lmm6_fixed<-lmer(deltaPI ~ Selection_Coef + log(Evolv_Gen) + (1 | Gene_ID) + (1 | Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0 &
                                                                                                                 slim_dd_GEN4$Evolv_Gen > 0,])
# Evaluate
drop1(lmm6_fixed)

# plot
plot(allEffects(lmm6_fixed))
summary(lmm6_fixed)
# Selection has strongest effect but also strong effect of target of selection

# Plot
deltaPI_selection_plot<-slim_dd_GEN4 %>% group_by(Gene_ID,Migration) %>% summarise_all(.funs=mean)
deltaPI_selection_plot2<-slim_dd_GEN4 %>% group_by(Gene_ID,Migration,Treatment_Demo) %>% summarise_all(.funs=mean)
deltaPI_selection_fig<-ggplot(deltaPI_selection_plot,aes(Selection_Coef,deltaPI))+
  geom_point(size=4,shape=17)+
  geom_smooth(method="lm")+
  geom_point(data=deltaPI_selection_plot2,alpha=0.2)+
  facet_wrap(~Migration)+
  theme_bw()+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        strip.text = element_text(size=18))+
  labs(y=expression(Delta~pi),x=expression(Selection~Coefficient~(log[10])))


# Plot these together
#deltaPI_plots<-ggarrange(deltaPI_selection_fig,deltaPI_target_fig,ncol=1,nrow=2)

# Can we get a better model by squaring selection
# deltaPI now for treatments WITHOUT migration
lmm7<-lmer(deltaPI~Selection_Coef^2+log(exon_N)+log(Target_sel)+log(gene_length)+log(Evolv_Gen)+(1|Gene_ID)+(1|Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0 &
                                                                                                                                              slim_dd_GEN4$Evolv_Gen > 0,])
drop1(lmm7)
model.step7<-step(lmm7,scope=~.^2,direction="both")
lmm7_fixed<-lmer(deltaPI ~ Selection_Coef^2 + log(Evolv_Gen)  + (1 | Gene_ID) + (1 | Treatment_Demo),slim_dd_GEN4[slim_dd_GEN4$Migration == 0 &
                                                                                                                    slim_dd_GEN4$Evolv_Gen > 0,])

# plot
plot(allEffects(lmm7_fixed))
summary(lmm7_fixed)

deltaPI_selection_fig<-ggplot(deltaPI_selection_plot,aes(Selection_Coef^2,deltaPI))+
  geom_point(size=4,shape=17)+
  geom_smooth(method="lm")+
  geom_point(data=deltaPI_selection_plot2,alpha=0.2)+
  facet_wrap(~Migration)+
  theme_bw()+
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        strip.text = element_text(size=18))+
  labs(y=expression(Delta~pi),x=expression((Selection~Coefficient~(log[10])^2)))


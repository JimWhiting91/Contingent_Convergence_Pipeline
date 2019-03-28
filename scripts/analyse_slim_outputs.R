#######################
# THIS SCRIPT TAKES SLIM OUTPUTS FOR LOOKING AT FST-DXY-delPI UNDER DIFFERENT DEMOGRAPHIC MODELS

# STEP 1: READ AND CLEAN
# STEP 2: ANALYSIS AND FIG PRODUCTION

#######################

# Packages
lib=c("tidyr","dplyr","pbapply","randomForest","VennDiagram","ghibli","ggpubr","data.table","ggplot2","ggedit","parallel","doParallel","effects","knitr","lme4","MASS","lmerTest","dplyr","rpart","partykit","rattle","randomForest")
lapply(lib,library,character.only=T)

# Source functions
source("R/Slim_General_Functions.R")

# Set up parallel back-end
cl <- makeCluster(detectCores())
registerDoParallel(cl)

####### STEP 1 - Set variables first ########

# Define inputs
pop_size<-1000

# Run over all mutation rates. These values don't affect the analysis, but name future directories and outputs
# The analysis is run over all mutation rates, producing figures and outputs for all steps

rates<-c(-5,-6,-7)

# Where is the data? Data should be in the demographics_model dir
data_dirs<-c("demographics_model/path_1",
			 "demographics_model/path_2",
			 "demographics_model/path_3")
			 
#######################             


# Loop over all mutation rates, comment out to do individual analyses
for(i in 1:length(rates)){
  mut<-i

# Make output/figs directories
dir.create(file.path(paste0("outputs/slim_output_",mut_rate)))
dir.create(file.path(paste0("figs/slim_output_",mut_rate)))

# Read in the data
N_Outputs<-list.files(paste0("data/",data_dirs[mut]))

# Remove the gene 2_Chunks
to_remove<-grep("2_Chunk",N_Outputs)
N_Outputs<-N_Outputs[(max(to_remove)+1):length(N_Outputs)]

# Read in
slim_dd<-data.frame(rbindlist(mclapply(1:length(N_Outputs),function(x){
  tmp_in<-fread(paste0("data/",data_dirs[mut],"/",N_Outputs[x]),fill=T)
  return(tmp_in)
}),fill = T))

# Remove incorrect runs
slim_dd<-slim_dd[slim_dd$Bottleneck != 10000 & slim_dd$Selection_Coef %in% c(0.1,0.5,1.0,5.0,10.0),]

# Make a treatment column
slim_dd$Treatment_All<-paste0("Sel=",slim_dd$Selection_Coef,"/Bot=",slim_dd$Bottleneck,"/Pop2=",slim_dd$Pop2_Size,"/Mig=",slim_dd$Migration)
slim_dd$Treatment_All<-as.factor(slim_dd$Treatment_All)

# No Selection
slim_dd$Treatment_Demo<-paste0("Bot=",slim_dd$Bottleneck,"/Pop2=",slim_dd$Pop2_Size,"/Mig=",slim_dd$Migration)
slim_dd$Treatment_Demo<-as.factor(slim_dd$Treatment_Demo)

# Calculate deltaPI which is log(10)-transformed
# Add the bare minimum increase to PI, which is 1 more site in 1 more individual ie. 1/pop_size
slim_dd$deltaPI<-log10((slim_dd$Pop1_PI+(1/pop_size))/(slim_dd$Pop2_PI+(1/pop_size)))

# Transform selection_coef so that it scales positively with selection
slim_dd$Selection_Coef<-1/(slim_dd$Selection_Coef)

# Targets of selection
slim_dd$Target_sel<-slim_dd$Exon_Length_Total/slim_dd$Gene_Size

# Find the Exon N from the gene files 
# Read in the data
N_Outputs<-list.files(paste0("data/",data_dirs[mut]))

# Remove the gene 2_Chunks
to_keep<-grep("2_Chunk",N_Outputs)
N_Outputs<-N_Outputs[1:max(to_keep)]

tidying<-as.data.frame(N_Outputs) %>% separate(N_Outputs, into = paste("V", 1:5, sep = "_"))

# Remove genes with no data
tidying2<-tidying[tidying$V_5 %in% unique(slim_dd$Gene_ID),]
N_Outputs<-paste(tidying2[,1],tidying2[,2],tidying2[,3],tidying2[,4],tidying2[,5],sep="_")

# Read in
genes_list<-lapply(1:length(N_Outputs),function(x){
  tmp_in<-fread(paste0("data/",data_dirs[mut],"/",N_Outputs[x]),fill=T)
  exon_N<-(length(tmp_in$V1)-1)/2
  return(exon_N)
})

exons_dd<-data.frame(exons=unlist(genes_list),
                     id=unique(slim_dd$Gene_ID))

# Add to slim_dd 32*4 for all treatments and all timepoints
for(i in 1:nrow(exons_dd)){
  slim_dd$Exon_N[i]<-exons_dd[exons_dd$id == slim_dd$Gene_ID[i],"exons"]
}


#######################			 STEP 2					####################################


# ------------------------------------------
#### Various Correlation Analyses #####
# ------------------------------------------

#-------------------------------------------
# Calculate correlation matrices
#-------------------------------------------

cor_dd<-data.frame(Treatment = levels(slim_dd$Treatment_Demo))

# Calculate correlations over each generation
cor_dd_list<-list()
for (i in 1:length(cor_dd$Treatment)){
  # Treatments
  cor_dd$Mig[i]<-mean(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"Migration"])
  cor_dd$Bot[i]<-mean(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"Bottleneck"])
  cor_dd$Pop2[i]<-mean(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"Pop2_Size"])
  
  
  # Correlations within generation
  for (j in 1:4){
    
    # Between measures
    cor_dd$FST.dXY[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"FST"])),
                           as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"dXY"])),method = c("pearson"))
    cor_dd$FST.delPI[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"FST"])),
                             as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"deltaPI"])),method = c("pearson"))
    cor_dd$dXY.delPI[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"dXY"])),
                             as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"deltaPI"])),method = c("pearson"))
    
    # With Selection
    cor_dd$FST.SEL[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"FST"])),
                           as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"Selection_Coef"])),method = c("pearson"))
    cor_dd$dXY.SEL[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"dXY"])),
                           as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"Selection_Coef"])),method = c("pearson"))
    cor_dd$delPI.SEL[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"deltaPI"])),
                             as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"Selection_Coef"])),method = c("pearson"))
    
    # With founder pi
    cor_dd$FST.SEL[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"FST"])),
                           as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == "Gen_1","Pop1_PI"])),method = c("pearson"))
    cor_dd$dXY.SEL[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"dXY"])),
                           as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == "Gen_1","Pop1_PI"])),method = c("pearson"))
    cor_dd$delPI.SEL[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == paste0("Gen_",j),"deltaPI"])),
                             as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i] & slim_dd$Generation == "Gen_1","Pop1_PI"])),method = c("pearson"))
    
    cor_dd_list[[j]]<-cor_dd
    
  }
  
  # Correlations across all generations
  cor_dd$FST.dXY[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"FST"])),
                         as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"dXY"])),method = c("pearson"))
  cor_dd$FST.SEL[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"FST"])),
                         as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"Selection_Coef"])),method = c("pearson"))
  cor_dd$dXY.SEL[i]<-cor(as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"dXY"])),
                         as.numeric(as.character(slim_dd[slim_dd$Treatment_Demo == levels(slim_dd$Treatment_Demo)[i],"Selection_Coef"])),method = c("pearson"))
  cor_dd_list[[5]]<-cor_dd
  
}

#-------------------------------------------
# EFFECT OF TREATMENTS ON MEASURES
#-------------------------------------------

# Plot relationships across generations
# Figures to examine effects of demographic variables across generation
treatment_vec<-c("Bottleneck","Pop2_Size","Migration")

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  # Edit the data such that Fst and dXY are same column with identifier column
  # Summarise
  tmp1<-cbind(slim_dd[,c("Generation","FST",treatment_vec[x])],Measure=rep("FST"))
  # Remove any NA values from calcs
  tmp1<-na.omit(tmp1)
  tmp1<-as.data.frame(tmp1 %>% 
                        group_by(Generation,tmp1[,treatment_vec[x]]) %>%
                        dplyr::summarise(mean(FST),(sd(FST)/sqrt(length(FST)))))
  colnames(tmp1)<-c("Generation",treatment_vec[x],"Mean","SE")
  tmp1$Metric<-rep("FST")
  
  # Summarise dXY
  tmp2<-cbind(slim_dd[,c("Generation","dXY",treatment_vec[x])],Measure=rep("dXY"))
  tmp2<-as.data.frame(tmp2 %>% 
                        group_by(Generation,slim_dd[,treatment_vec[x]]) %>%
                        dplyr::summarise(mean(dXY),(sd(dXY)/sqrt(length(dXY)))))
  colnames(tmp2)<-c("Generation",treatment_vec[x],"Mean","SE")
  tmp2$Metric<-rep("dXY")
  
  # Summarise deltaPI
  tmp3<-cbind(slim_dd[,c("Generation","deltaPI",treatment_vec[x])],Measure=rep("deltaPI"))
  tmp3<-as.data.frame(tmp3 %>% 
                        group_by(Generation,slim_dd[,treatment_vec[x]]) %>%
                        dplyr::summarise(mean(deltaPI),(sd(deltaPI)/sqrt(length(deltaPI)))))
  colnames(tmp3)<-c("Generation",treatment_vec[x],"Mean","SE")
  tmp3$Metric<-rep("deltaPI")
  
  
  # Rbind
  tmp<-rbind(tmp1,tmp2,tmp3)
  tmp[,2]<-as.factor(tmp[,2])
  
  tmp$Metric<-factor(tmp$Metric,levels=c("FST","dXY","deltaPI"))
  levels(tmp$Metric)<-c("F[ST]","D[XY]","Delta*pi")
  
  
  # Plot
  pd <- position_dodge(0.5)
  
  g1<-ggplot(tmp,aes(x=Generation,y=Mean,shape=tmp[,2],colour=tmp[,2]))+
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=pd,show.legend = FALSE) +
    geom_point(position=pd,size=4)+
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
    scale_x_discrete(labels=c(0.05*2*pop_size,
                              0.158*2*pop_size,
                              1.581*2*pop_size,
                              5*2*pop_size))+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium"))) +
    xlab("Generation")+
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab("Divergence")
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="Contraction",shape="Contraction")
  }
  
  return(g1)
})

pdf(paste0("figs/slim_output_",mut_rate,"/Spatial_Model_Divergence_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height=10)
print(ggarrange(plot_list[[1]],
          plot_list[[2]],
          plot_list[[3]],ncol = 3,labels = "AUTO"))
dev.off()

#-------------------------------------------
# CORRELATIONS WITH SELECTION
#-------------------------------------------

# Figures to examine effects of demographic variables across generations on correlation coefficients with selection
treatment_vec<-c("Bottleneck","Pop2_Size","Migration")

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  # Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column
  # Summarise FST
  tmp1<-cbind(slim_dd[,c("Generation","FST","Selection_Coef",treatment_vec[x])],Measure=rep("FST"))
  tmp1<-as.data.frame(tmp1 %>% 
                        group_by(Generation,tmp1[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(FST,Selection_Coef)$estimate),
                                         (cor.test(FST,Selection_Coef)$conf.int[1]),
                                         (cor.test(FST,Selection_Coef)$conf.int[2])))
  colnames(tmp1)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp1$Metric<-rep("FST")
  
  # Summarise dXY
  tmp2<-cbind(slim_dd[,c("Generation","dXY","Selection_Coef",treatment_vec[x])],Measure=rep("dXY"))
  tmp2<-as.data.frame(tmp2 %>% 
                        group_by(Generation,tmp2[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(dXY,Selection_Coef)$estimate),
                                         (cor.test(dXY,Selection_Coef)$conf.int[1]),
                                         (cor.test(dXY,Selection_Coef)$conf.int[2])))
  colnames(tmp2)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp2$Metric<-rep("dXY")
  
  # Summarise deltaPI
  tmp3<-cbind(slim_dd[,c("Generation","deltaPI","Selection_Coef",treatment_vec[x])],Measure=rep("deltaPI"))
  tmp3<-as.data.frame(tmp3 %>% 
                        group_by(Generation,tmp3[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(deltaPI,Selection_Coef)$estimate),
                                         (cor.test(deltaPI,Selection_Coef)$conf.int[1]),
                                         (cor.test(deltaPI,Selection_Coef)$conf.int[2])))
  colnames(tmp3)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp3$Metric<-rep("deltaPI")
  
  
  # Rbind
  tmp<-rbind(tmp1,tmp2,tmp3)
  tmp[,2]<-as.factor(tmp[,2])
  
  # Plot
  pd <- position_dodge(0.5)
  
  tmp$Metric<-factor(tmp$Metric,levels=c("FST","dXY","deltaPI"))
  levels(tmp$Metric)<-c("F[ST]","D[XY]","Delta*pi")
  
  
  g1<-ggplot(tmp,aes(x=Generation,y=Cor.Coef,shape=tmp[,2],colour=tmp[,2]))+
    geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, position=pd,show.legend = FALSE) +
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
    scale_x_discrete(labels=c(0.05*2*pop_size,
                              0.158*2*pop_size,
                              1.581*2*pop_size,
                              5*2*pop_size))+
    xlab("Generation")+
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab("Correlation with Selection")
  
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="Contraction",shape="Contraction")
  }
  
  
  return(list(g1,tmp))
})

pdf(paste0("figs/slim_output_",mut_rate,"/Spatial_Model_DivergenceSelectionCoef_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height=10)
print(ggarrange(plot_list[[1]][[1]],
          plot_list[[2]][[1]],
          plot_list[[3]][[1]],ncol = 3,labels = "AUTO"))
dev.off()

# We also want to save the correlation coefficients

# Bottlenecks
write.table(plot_list[[1]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergenceSelectionCoef_Bottleneck.txt"),
            row.names = F,quote = F,sep="\t")
# Contractions
write.table(plot_list[[2]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergenceSelectionCoef_Contractions.txt"),
            row.names = F,quote = F,sep="\t")
# Migration
write.table(plot_list[[3]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergenceSelectionCoef_Migration.txt"),
            row.names = F,quote = F,sep="\t")

# Plot selection coefficients for all 32 treatments in 1-1 comparison to aid visibility

# Reformat cor_dd so its by generation
# Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column

cor_dd2<-cbind(slim_dd[,c("Treatment_Demo","Generation","FST","Selection_Coef","dXY","deltaPI")])
cor_dd2<-as.data.frame(cor_dd2 %>% 
                         group_by(Generation,Treatment_Demo) %>%
                         dplyr::summarise((cor.test(FST,Selection_Coef)$estimate),
                                          (cor.test(dXY,Selection_Coef)$estimate),
                                          (cor.test(deltaPI,Selection_Coef)$estimate)))
colnames(cor_dd2)<-c("Generation","Treatment","FST.cor","Dxy.cor","deltaPI.cor")

# Reformat
cor_dd3<-data.frame(Generation=rep(cor_dd2$Generation,3),
                    Treatment=rep(cor_dd2$Treatment,3),
                    Correlation1=c(cor_dd2$FST.cor,cor_dd2$FST.cor,cor_dd2$Dxy.cor),
                    Correlation2=c(cor_dd2$Dxy.cor,cor_dd2$deltaPI.cor,cor_dd2$deltaPI.cor),
                    Mig=rep(cor_dd$Mig,12),
                    Pop2=rep(cor_dd$Pop2,12),
                    Bot=rep(cor_dd$Bot,12),
                    Measure=rep(c("F[ST]~x~D[XY]","F[ST]~x~Delta*pi","D[XY]~x~Delta*pi"),each=128))
cor_dd3$Generation <- as.character(cor_dd3$Generation)

# Add in years
cor_dd3[cor_dd3$Generation=="Gen_1","Generation"]<-paste0(0.05*2*pop_size,"~Gens")
cor_dd3[cor_dd3$Generation=="Gen_2","Generation"]<-paste0(0.158*2*pop_size,"~Gens")
cor_dd3[cor_dd3$Generation=="Gen_3","Generation"]<-paste0(1.581*2*pop_size,"~Gens")
cor_dd3[cor_dd3$Generation=="Gen_4","Generation"]<-paste0(5*2*pop_size,"~Gens")

# Plot
cor_dd3$Measure<-factor(cor_dd3$Measure,levels=c("F[ST]~x~D[XY]","F[ST]~x~Delta*pi","D[XY]~x~Delta*pi"))
cor_dd3$Generation<-factor(cor_dd3$Generation,levels=c(paste0(0.05*2*pop_size,"~Gens"),
                                                       paste0(0.158*2*pop_size,"~Gens"),
                                                       paste0(1.581*2*pop_size,"~Gens"),
                                                       paste0(5*2*pop_size,"~Gens")))

# Plot with colour Contraction and migration shape
cor_dd3$Mig<-as.factor(cor_dd3$Mig)
cor_dd3$Pop2<-as.factor(cor_dd3$Pop2)

g1<-ggplot(data=cor_dd3,aes(x=Correlation1,y=Correlation2,colour=Pop2,shape=Mig))+
  geom_point(size=2.5)+
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
  labs(colour="Contraction",x="Correlation 1",y="Correlation 2",shape="Migration")+
  scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium")))

pdf(paste0("figs/slim_output_",mut_rate,"/Spatial_Model_Comparing_Correlations_with_Selection.pdf"),width=13,height=10)
print(g1)
dev.off()


#-------------------------------------------
# CORRELATIONS WITH PI
#-------------------------------------------


# Figures to examine effects of demographic variables across generations on correlation coefficients with selection
treatment_vec<-c("Bottleneck","Pop2_Size","Migration")

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  # Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column
  # Summarise FST
  tmp1<-cbind(slim_dd[,c("Generation","FST",treatment_vec[x])],Measure=rep("FST"))
  tmp1$Pop1_PI<-rep(slim_dd[slim_dd$Generation == "Gen_1","Pop1_PI"],each=4)
  tmp1<-as.data.frame(tmp1 %>% 
                        group_by(Generation,tmp1[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(FST,Pop1_PI)$estimate),
                                         (cor.test(FST,Pop1_PI)$conf.int[1]),
                                         (cor.test(FST,Pop1_PI)$conf.int[2])))
  colnames(tmp1)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp1$Metric<-rep("FST")
  
  # Summarise dXY
  tmp2<-cbind(slim_dd[,c("Generation","dXY",treatment_vec[x])],Measure=rep("dXY"))
  tmp2$Pop1_PI<-rep(slim_dd[slim_dd$Generation == "Gen_1","Pop1_PI"],each=4)
  tmp2<-as.data.frame(tmp2 %>% 
                        group_by(Generation,tmp2[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(dXY,Pop1_PI)$estimate),
                                         (cor.test(dXY,Pop1_PI)$conf.int[1]),
                                         (cor.test(dXY,Pop1_PI)$conf.int[2])))
  colnames(tmp2)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp2$Metric<-rep("dXY")
  
  # Summarise deltaPI
  tmp3<-cbind(slim_dd[,c("Generation","deltaPI",treatment_vec[x])],Measure=rep("deltaPI"))
  tmp3$Pop1_PI<-rep(slim_dd[slim_dd$Generation == "Gen_1","Pop1_PI"],each=4)
  tmp3<-as.data.frame(tmp3 %>% 
                        group_by(Generation,tmp3[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(deltaPI,Pop1_PI)$estimate),
                                         (cor.test(deltaPI,Pop1_PI)$conf.int[1]),
                                         (cor.test(deltaPI,Pop1_PI)$conf.int[2])))
  colnames(tmp3)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp3$Metric<-rep("deltaPI")
  
  
  # Rbind
  tmp<-rbind(tmp1,tmp2,tmp3)
  tmp[,2]<-as.factor(tmp[,2])
  
  # Plot
  pd <- position_dodge(0.5)
  tmp$Metric<-factor(tmp$Metric,levels=c("FST","dXY","deltaPI"))
  levels(tmp$Metric)<-c("F[ST]","D[XY]","Delta*pi")
  
  
  g1<-ggplot(tmp,aes(x=Generation,y=Cor.Coef,shape=tmp[,2],colour=tmp[,2]))+
    geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, position=pd,show.legend = FALSE) +
    geom_point(position=pd,size=3)+
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
    scale_x_discrete(labels=c(0.05*2*pop_size,
                              0.158*2*pop_size,
                              1.581*2*pop_size,
                              5*2*pop_size))+
    xlab("Generation")+
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab(expression("Correlation with Founding "~pi))+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium")))
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="Contraction",shape="Contraction")
  }
  
  return(list(g1,tmp))
})

pdf(paste0("figs/slim_output_",mut_rate,"/Spatial_Model_DivergencePiCoef_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height=10)
print(ggarrange(plot_list[[1]][[1]],
          plot_list[[2]][[1]],
          plot_list[[3]][[1]],ncol = 3,labels = "AUTO"))
dev.off()

# We also want to save the correlation coefficients

# Bottlenecks
write.table(plot_list[[1]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergencePiCoef_Bottleneck.txt"),
            row.names = F,quote = F,sep="\t")
# Contractions
write.table(plot_list[[2]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergencePiCoef_Contractions.txt"),
            row.names = F,quote = F,sep="\t")
# Migration
write.table(plot_list[[3]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergencePiCoef_Migration.txt"),
            row.names = F,quote = F,sep="\t")

#------------------------------
# Correlation of measures
#------------------------------

# Figures to examine effects of demographic variables across generations on correlation coefficients between measures
treatment_vec<-c("Bottleneck","Pop2_Size","Migration")

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  # Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column
  # Summarise FST
  tmp1<-cbind(slim_dd[,c("Generation","FST","dXY",treatment_vec[x])],Measure=rep("FST x DXY"))
  tmp1<-as.data.frame(tmp1 %>% 
                        group_by(Generation,slim_dd[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(FST,dXY)$estimate),
                                         (cor.test(FST,dXY)$conf.int[1]),
                                         (cor.test(FST,dXY)$conf.int[2])))
  colnames(tmp1)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp1$Metric<-rep("FST x DXY")
  
  # Summarise dXY
  tmp2<-cbind(slim_dd[,c("Generation","FST","deltaPI",treatment_vec[x])],Measure=rep("FST x deltaPI"))
  tmp2<-as.data.frame(tmp2 %>% 
                        group_by(Generation,slim_dd[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(FST,deltaPI)$estimate),
                                         (cor.test(FST,deltaPI)$conf.int[1]),
                                         (cor.test(FST,deltaPI)$conf.int[2])))
  colnames(tmp2)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp2$Metric<-rep("FST x deltaPI")
  
  # Summarise deltaPI
  tmp3<-cbind(slim_dd[,c("Generation","dXY","deltaPI",treatment_vec[x])],Measure=rep("dXY x deltaPI"))
  tmp3<-as.data.frame(tmp3 %>% 
                        group_by(Generation,slim_dd[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(deltaPI,dXY)$estimate),
                                         (cor.test(deltaPI,dXY)$conf.int[1]),
                                         (cor.test(deltaPI,dXY)$conf.int[2])))
  colnames(tmp3)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp3$Metric<-rep("dXY x deltaPI")
  
  
  # Rbind
  tmp<-rbind(tmp1,tmp2,tmp3)
  tmp[,2]<-as.factor(tmp[,2])
  
  # Plot
  pd <- position_dodge(0.5)
  
  tmp$Metric<-factor(tmp$Metric,levels=c("FST x DXY","FST x deltaPI","dXY x deltaPI"))
  levels(tmp$Metric)<-c("F[ST]~x~D[XY]","F[ST]~x~Delta*pi","D[XY]~x~Delta*pi")
  
  
  g1<-ggplot(tmp,aes(x=Generation,y=Cor.Coef,shape=tmp[,2],colour=tmp[,2]))+
    geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, position=pd,show.legend = FALSE) +
    geom_point(position=pd,size=3)+
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
    scale_x_discrete(labels=c(0.05*2*pop_size,
                              0.158*2*pop_size,
                              1.581*2*pop_size,
                              5*2*pop_size))+
    xlab("Generation")+
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab("Correlation between measures")+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium")))
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="Contraction",shape="Contraction")
  }
  
  return(list(g1,tmp))
})

pdf(paste0("figs/slim_output_",mut_rate,"/Spatial_model_DivergenceCoefs_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height = 10)
print(ggarrange(plot_list[[1]][[1]],
          plot_list[[2]][[1]],
          plot_list[[3]][[1]],ncol = 3,labels = "AUTO"))
dev.off()

# We also want to save the correlation coefficients

# Bottlenecks
write.table(plot_list[[1]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergenceCoefs_Bottleneck.txt"),
            row.names = F,quote = F,sep="\t")
# Contractions
write.table(plot_list[[2]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergenceCoefs_Contractions.txt"),
            row.names = F,quote = F,sep="\t")
# Migration
write.table(plot_list[[3]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergenceCoefs_Migration.txt"),
            row.names = F,quote = F,sep="\t")

#-----------------------------------
# Correlations with dynamic PI in diverging pop
#-----------------------------------

plot_list<-lapply(1:length(treatment_vec),function(x){
  
  # Edit the data such that Fst.sel and dXY.sel and delPI.sel are same column with identifier column
  # Summarise FST
  tmp1<-cbind(slim_dd[,c("Generation","FST","Pop2_PI",treatment_vec[x])],Measure=rep("FST"))
  tmp1<-as.data.frame(tmp1 %>% 
                        group_by(Generation,tmp1[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(FST,Pop2_PI)$estimate),
                                         (cor.test(FST,Pop2_PI)$conf.int[1]),
                                         (cor.test(FST,Pop2_PI)$conf.int[2])))
  colnames(tmp1)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp1$Metric<-rep("FST")
  
  # Summarise dXY
  tmp2<-cbind(slim_dd[,c("Generation","dXY","Pop2_PI",treatment_vec[x])],Measure=rep("dXY"))
  tmp2$Pop1_PI<-rep(slim_dd[slim_dd$Generation == "Gen_1","Pop1_PI"],each=4)
  tmp2<-as.data.frame(tmp2 %>% 
                        group_by(Generation,tmp2[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(dXY,Pop2_PI)$estimate),
                                         (cor.test(dXY,Pop2_PI)$conf.int[1]),
                                         (cor.test(dXY,Pop2_PI)$conf.int[2])))
  colnames(tmp2)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp2$Metric<-rep("dXY")
  
  # Summarise deltaPI
  tmp3<-cbind(slim_dd[,c("Generation","deltaPI","Pop2_PI",treatment_vec[x])],Measure=rep("deltaPI"))
  tmp3$Pop1_PI<-rep(slim_dd[slim_dd$Generation == "Gen_1","Pop1_PI"],each=4)
  tmp3<-as.data.frame(tmp3 %>% 
                        group_by(Generation,tmp3[,treatment_vec[x]]) %>%
                        dplyr::summarise((cor.test(deltaPI,Pop2_PI)$estimate),
                                         (cor.test(deltaPI,Pop2_PI)$conf.int[1]),
                                         (cor.test(deltaPI,Pop2_PI)$conf.int[2])))
  colnames(tmp3)<-c("Generation",treatment_vec[x],"Cor.Coef","LL","UL")
  tmp3$Metric<-rep("deltaPI")
  
  
  # Rbind
  tmp<-rbind(tmp1,tmp2,tmp3)
  tmp[,2]<-as.factor(tmp[,2])
  
  # Plot
  pd <- position_dodge(0.5)
  tmp$Metric<-factor(tmp$Metric,levels=c("FST","dXY","deltaPI"))
  levels(tmp$Metric)<-c("F[ST]","D[XY]","Delta*pi")
  
  
  g1<-ggplot(tmp,aes(x=Generation,y=Cor.Coef,shape=tmp[,2],colour=tmp[,2]))+
    geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, position=pd,show.legend = FALSE) +
    geom_point(position=pd,size=3)+
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
    scale_x_discrete(labels=c(0.05*2*pop_size,
                              0.158*2*pop_size,
                              1.581*2*pop_size,
                              5*2*pop_size))+
    xlab("Generation")+
    labs(colour=treatment_vec[x],shape=treatment_vec[x])+
    ylab(expression("Correlation with Diverging "~pi))+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium")))
  
  if(treatment_vec[x]=="Pop2_Size"){
    g1<-g1+labs(colour="Contraction",shape="Contraction")
  }
  
  
  return(list(g1,tmp))
})

pdf(paste0("figs/slim_output_",mut_rate,"/Spatial_Model_DivergencePi2Coef_by_Generation_for_treatments_FST_DXY_DELTAPI.pdf"),width=16,height=10)
print(ggarrange(plot_list[[1]][[1]],
          plot_list[[2]][[1]],
          plot_list[[3]][[1]],ncol = 3,labels = "AUTO"))
dev.off()

# We also want to save the correlation coefficients

# Bottlenecks
write.table(plot_list[[1]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergencePi2Coef_Bottleneck.txt"),
            row.names = F,quote = F,sep="\t")
# Contractions
write.table(plot_list[[2]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergencePi2Coef_Contractions.txt"),
            row.names = F,quote = F,sep="\t")
# Migration
write.table(plot_list[[3]][[2]],
            paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_DivergencePi2Coef_Migration.txt"),
            row.names = F,quote = F,sep="\t")

##### Correlating PI, Selection, and contractions #####

# What are the correlation coefs for demo treatments on pi
contraction_pi_dd<-as.data.frame(slim_dd %>% 
                                   group_by(Generation,Selection_Coef) %>%
                                   dplyr::summarise((cor.test(log(Pop2_Size),Pop2_PI)$estimate),
                                                    (cor.test(log(Pop2_Size),Pop2_PI)$conf.int[1]),
                                                    (cor.test(log(Pop2_Size),Pop2_PI)$conf.int[2])))

# Add colours to figure
plot_dd<-data.frame(slim_dd,
                    cor_colour=NA)
for(i in 1:length(contraction_pi_dd$Generation)){
  plot_dd[plot_dd$Generation == contraction_pi_dd$Generation[i] & plot_dd$Selection_Coef == contraction_pi_dd$Selection_Coef[i],"cor_colour"]<-contraction_pi_dd$`(cor.test(log(Pop2_Size), Pop2_PI)$estimate)`[i]
}


plot_dd[plot_dd$Generation=="Gen_1","Generation"]<-paste0(0.05*2*pop_size," Gens")
plot_dd[plot_dd$Generation=="Gen_2","Generation"]<-paste0(0.158*2*pop_size," Gens")
plot_dd[plot_dd$Generation=="Gen_3","Generation"]<-paste0(1.581*2*pop_size," Gens")
plot_dd[plot_dd$Generation=="Gen_4","Generation"]<-paste0(5*2*pop_size," Gens")


contraction_pi<-ggplot(plot_dd,aes(y=Pop2_PI,x=log(Pop2_Size)))+
  geom_jitter(size=0.5,alpha=0.2,aes(colour=cor_colour))+
  scale_colour_gradient2(high="red2")+
  geom_smooth(method="lm",colour="black")+
  facet_grid(Generation~Selection_Coef,
             scales= "free_y")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        strip.text=element_text(size=14))+
  labs(x="Contraction (log)",
       y=expression("Diverging"~pi),
       colour="Corr Coef")

## Migration
migration_pi_dd<-as.data.frame(slim_dd %>% 
                                 group_by(Generation,Selection_Coef) %>%
                                 dplyr::summarise((cor.test(Migration,Pop2_PI)$estimate),
                                                  (cor.test(Migration,Pop2_PI)$conf.int[1]),
                                                  (cor.test(Migration,Pop2_PI)$conf.int[2])))

# Add colours to figure
plot_dd<-data.frame(slim_dd,
                    cor_colour=NA)
for(i in 1:length(migration_pi_dd$Generation)){
  plot_dd[plot_dd$Generation == migration_pi_dd$Generation[i] & plot_dd$Selection_Coef == migration_pi_dd$Selection_Coef[i],"cor_colour"]<-migration_pi_dd$`(cor.test(Migration, Pop2_PI)$estimate)`[i]
}


plot_dd[plot_dd$Generation=="Gen_1","Generation"]<-paste0(0.05*2*pop_size," Gens")
plot_dd[plot_dd$Generation=="Gen_2","Generation"]<-paste0(0.158*2*pop_size," Gens")
plot_dd[plot_dd$Generation=="Gen_3","Generation"]<-paste0(1.581*2*pop_size," Gens")
plot_dd[plot_dd$Generation=="Gen_4","Generation"]<-paste0(5*2*pop_size," Gens")

plot_dd$Generation_f = factor(plot_dd$Generation, levels=c(paste0(0.05*2*pop_size," Gens"),
                                                           paste0(0.158*2*pop_size," Gens"),
                                                           paste0(1.581*2*pop_size," Gens"),
                                                           paste0(5*2*pop_size," Gens")))

migration_pi<-ggplot(plot_dd,aes(y=Pop2_PI,x=factor(Migration)))+
  #geom_jitter(size=0.5,alpha=0.2)+
  #geom_smooth(method="lm",colour="black")+
  geom_violin(aes(fill=cor_colour),draw_quantiles = c(0.5))+
  #stat_summary(fun.data = median, geom = "bar")+
  facet_grid(Generation_f~Selection_Coef,
             scales= "free_y")+
  scale_fill_gradient2(high="red2")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        strip.text=element_text(size=14))+
  labs(x="Migration",
       y=expression("Diverging"~pi),
       colour="Corr Coef")

pdf(paste0("figs/slim_output_",mut_rate,"/Relationships_between_demography_and_pi.pdf"),width=10,height=14)
print(contraction_pi)
print(migration_pi)
dev.off()

##### Correlation Matrices #####
# Run over all 32 treatments
treatments_to_keep<-unique(slim_dd$Treatment_Demo)

for (level_gen in levels(as.factor(slim_dd$Generation))) {
  slim_dd_GEN4<-slim_dd[slim_dd$Generation == level_gen,]
  slim_dd_GEN4<-na.omit(slim_dd_GEN4)
  
  FST_cov<-vector()
  foreach(i=1:length(treatments_to_keep)) %do% {
    tmp<-as.vector(slim_dd_GEN4[slim_dd_GEN4$Treatment_Demo == treatments_to_keep[i],"FST"])
    FST_cov<-cbind(FST_cov,tmp)
    return(FST_cov)
  }
  colnames(FST_cov)<-treatments_to_keep
  
  dXY_cov<-vector()
  foreach(i=1:length(treatments_to_keep)) %do% {
    tmp<-as.vector(slim_dd_GEN4[slim_dd_GEN4$Treatment_Demo == treatments_to_keep[i],"dXY"])
    dXY_cov<-cbind(dXY_cov,tmp)
    return(dXY_cov)
  }
  colnames(dXY_cov)<-treatments_to_keep
  
  deltaPI_cov<-vector()
  foreach(i=1:length(treatments_to_keep)) %do% {
    tmp<-as.vector(slim_dd_GEN4[slim_dd_GEN4$Treatment_Demo == treatments_to_keep[i],"deltaPI"])
    deltaPI_cov<-cbind(deltaPI_cov,tmp)
    return(deltaPI_cov)
  }
  colnames(deltaPI_cov)<-treatments_to_keep
  
  # Plot matrices
  pdf(paste0("figs/slim_output_",mut_rate,"/Spatial_model_",level_gen,"_FST_dXY_deltaPI_Correlation_matrices.pdf"),width=16,height=16)
  print(plot_cor_mat(as.matrix(FST_cov),cor = TRUE))
  print(plot_cor_mat(as.matrix(dXY_cov),cor = TRUE))
  print(plot_cor_mat(as.matrix(deltaPI_cov),cor = TRUE))
  dev.off()
}

# ------------------------------------------
#### Proportion of Outliers #####
# ------------------------------------------

# Run over all generations
treatment_vec<-c("Bottleneck","Pop2_Size","Migration")
prop_figures<-lapply(1:4,function(gen){
  
  plot_dd<-data.frame(rbindlist(lapply(1:length(treatment_vec),function(x){
    
    # First split slim_dd into treatment levels
    treatment_slims<-data.frame(rbindlist(lapply(1:length(unique(slim_dd[,treatment_vec[x]])),function(y){
      
      tmp<-slim_dd[slim_dd$Generation==paste0("Gen_",gen),]
      tmp1<-data.frame(tmp,
                       treatment=tmp[,treatment_vec[x]])
      
      tmp<-tmp1[tmp1$treatment==unique(slim_dd[,treatment_vec[x]])[y],]
      
      # Mark 99% outliers
      FST_out<-tmp[tmp$FST > quantile(na.omit(tmp$FST),probs=0.99),]
      dXY_out<-tmp[tmp$dXY > quantile(tmp$dXY,probs=0.99),]
      pi_out<-tmp[tmp$deltaPI > quantile(tmp$deltaPI,probs=0.99),]
      
      outlier_list<-list(FST_out,dXY_out,pi_out)
      
      # Count proportion of outliers for each selection_coef
      measure_vec<-c("FST","DXY","deltaPI")
      outlier_props<-data.frame(rbindlist(lapply(1:length(outlier_list),function(m){
        
        sum <- data.frame( outlier_list[[m]] %>% 
                             group_by(Selection_Coef) %>%
                             summarise(no_rows = length(Selection_Coef)))
        
        sum$prop<-sum$no_rows/length(outlier_list[[m]][,1])
        sum$treatment<-unique(slim_dd[,treatment_vec[x]])[y]
        sum$measure<-measure_vec[m]
        
        return(sum)
      })))
      
      # Return
      return(outlier_props)
      
    })))
    
    treatment_slims$demo<-treatment_vec[x]
    
    return(treatment_slims)
  })))
  
  # Plot proportions
  # Edit plot_dd
  plot_dd[plot_dd$demo=="Pop2_Size","demo"]<-"Contraction"
  plot_dd$measure<-factor(plot_dd$measure,levels=c("FST","DXY","deltaPI"))
  levels(plot_dd$measure)<-c("F[ST]","D[XY]","Delta*pi")
  
  # Save Proportions to an RData file
  saveRDS(plot_dd,paste0("outputs/slim_output_",mut_rate,"/Spatial_Model_99Q_Proportion_by_Selection.rds"))
  
  # Calculate chi-square of non-randomness
  plot_dd$demo.treatment<-paste0(plot_dd$demo,".",plot_dd$treatment,".",plot_dd$measure)
  
  # Loop over each unique treatment.demo and perform chi
  plot_dd$chi<-NA
  treats<-unique(plot_dd$demo.treatment)
  for(i in 1:length(treats)){
    tmp<-plot_dd[plot_dd$demo.treatment == treats[i],"no_rows"]
    res<-chisq.test(tmp)
    
    if (res$p.value > 0.05){
      ast<-NA
    } else if (res$p.value > 0.01){
      ast<-"*"
    } else if (res$p.value > 0.001){
      ast<-"**"
    } else {
      ast<-"***"
    }
    plot_dd[plot_dd$demo.treatment == treats[i],"chi"]<-ast
  }
  
  # Plot with asterices
  
  prop_plot<-ggplot(plot_dd,aes(x=as.factor(treatment),y=prop,fill=as.factor(Selection_Coef)))+
    geom_col()+
    facet_grid(measure~demo,scales = "free_x",labeller = label_parsed)+
    scale_fill_manual(values = rev(ghibli_palette("PonyoMedium"))) +
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.title=element_text(size=32),
          axis.text=element_text(size=20),
          strip.text = element_text(size=24),
          legend.text = element_text(size=30),
          legend.title = element_text(size=32))+
    geom_text(aes(x=as.factor(treatment),y=1.1,label=chi),size=8)+
    scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1.0))+
    labs(y="Proportion of 99% Outliers",x="Treatment Level",fill="Selection")
  
  # Save proportions
  write.table(plot_dd,
              paste0("outputs/slim_output_",mut_rate,"/Prop_outliers_99q_Gen",gen,".txt"),
              quote = F,row.names = F,sep="\t")
  
  
  # Also plot for each treatment level (32 bars)
  
  plot_dd<-data.frame(rbindlist(lapply(1:length(unique(slim_dd$Treatment_Demo)),function(x){
    
    tmp<-slim_dd[slim_dd$Generation == paste0("Gen_",gen),]
    tmp<-tmp[tmp$Treatment_Demo == unique(tmp$Treatment_Demo)[x],]
    
    # Mark 99% outliers
    FST_out<-tmp[tmp$FST > quantile(na.omit(tmp$FST),probs=0.99),]
    dXY_out<-tmp[tmp$dXY > quantile(tmp$dXY,probs=0.99),]
    pi_out<-tmp[tmp$deltaPI > quantile(tmp$deltaPI,probs=0.99),]
    
    outlier_list<-list(FST_out,dXY_out,pi_out)
    
    # Count proportion of outliers for each selection_coef
    measure_vec<-c("FST","DXY","deltaPI")
    outlier_props<-data.frame(rbindlist(lapply(1:length(outlier_list),function(m){
      
      sum <- data.frame( outlier_list[[m]] %>% 
                           group_by(Selection_Coef) %>%
                           summarise(no_rows = length(Selection_Coef)))
      
      sum$prop<-sum$no_rows/length(outlier_list[[m]][,1])
      sum$treatment<-unique(slim_dd$Treatment_Demo)[x]
      sum$measure<-measure_vec[m]
      
      return(sum)
    })))
    return(outlier_props)
  })))
  
  # Plot
  # Edit plot_dd
  plot_dd$measure<-factor(plot_dd$measure,levels=c("FST","DXY","deltaPI"))
  levels(plot_dd$measure)<-c("F[ST]","D[XY]","Delta*pi")
  
  prop_plot2<-ggplot(plot_dd,aes(x=as.factor(treatment),y=prop,fill=as.factor(Selection_Coef)))+
    geom_col()+
    facet_wrap(~measure,ncol=1,strip.position = "right",labeller = label_parsed)+
    scale_fill_manual(values = rev(ghibli_palette("PonyoMedium"))) +
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.title=element_text(size=32),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=12,angle=45,hjust=1),
          strip.text = element_text(size=24),
          legend.text = element_text(size=30),
          legend.title = element_text(size=32))+
    labs(y="Proportion of 99% Outliers",x="Treatment",fill="Selection")
  
  # Save proportions
  write.table(plot_dd,
              paste0("outputs/slim_output_",mut_rate,"/Prop_outliers_99q_Gen",gen,"_allTreatments.txt"),
              quote = F,row.names = F,sep="\t")
  
  # Also plot for each treatment level without bottlenecks (16 bars)
  slim_dd$Treatment_noBot<-paste0("Pop2 = ",slim_dd$Pop2_Size,"/Mig = ",slim_dd$Migration)
  
  plot_list2<-lapply(1:length(unique(slim_dd$Treatment_noBot)),function(x){
    
    tmp<-slim_dd[slim_dd$Generation == paste0("Gen_",gen),]
    tmp<-tmp[tmp$Treatment_noBot == unique(tmp$Treatment_noBot)[x],]
    
    # Mark 99% outliers
    FST_out<-tmp[tmp$FST > quantile(na.omit(tmp$FST),probs=0.99),]
    dXY_out<-tmp[tmp$dXY > quantile(tmp$dXY,probs=0.99),]
    pi_out<-tmp[tmp$deltaPI > quantile(tmp$deltaPI,probs=0.99),]
    
    outlier_list<-list(FST_out,dXY_out,pi_out)
    
    # Count proportion of outliers for each selection_coef
    measure_vec<-c("FST","DXY","deltaPI")
    outlier_props<-data.frame(rbindlist(lapply(1:length(outlier_list),function(m){
      
      sum <- data.frame( outlier_list[[m]] %>% 
                           group_by(Selection_Coef) %>%
                           summarise(no_rows = length(Selection_Coef)))
      
      sum$prop<-sum$no_rows/length(outlier_list[[m]][,1])
      sum$treatment<-unique(slim_dd$Treatment_noBot)[x]
      sum$measure<-measure_vec[m]
      
      return(sum)
    })))
    return(list(outlier_props))
  })
  
  # Get each
  plot_dd<-data.frame(rbindlist(lapply(1:length(plot_list2),function(plot){return(plot_list2[[plot]][[1]])})))
  
  # Plot
  # Edit plot_dd
  plot_dd$measure<-factor(plot_dd$measure,levels=c("FST","DXY","deltaPI"))
  levels(plot_dd$measure)<-c("F[ST]","D[XY]","Delta*pi")
  
  prop_plot3<-ggplot(plot_dd,aes(x=as.factor(treatment),y=prop,fill=as.factor(Selection_Coef)))+
    geom_col()+
    facet_wrap(~measure,ncol=1,strip.position = "right",labeller = label_parsed)+
    scale_fill_manual(values = rev(ghibli_palette("PonyoMedium"))) +
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.title=element_text(size=32),
          axis.text.y=element_text(size=24),
          axis.text.x=element_text(size=22,angle=45,hjust=1),
          strip.text = element_text(size=24),
          legend.text = element_text(size=30),
          legend.title = element_text(size=32))+
    labs(y="Proportion of 99% Outliers",x="Treatment",fill="Selection")
  
  # Save proportions
  write.table(plot_dd,
              paste0("outputs/slim_output_",mut_rate,"/Prop_outliers_99q_Gen",gen,"_allTreatments_noBot.txt"),
              quote = F,row.names = F,sep="\t")
  
  # What are the characteristic features of outliers in each treatment
  
  # Save figs
  pdf(paste0("figs/slim_output_",mut_rate,"/Prop_outliers_bars_Gen",gen,".pdf"),width=18,height=16)
  print(prop_plot)
  print(prop_plot2)
  print(prop_plot3)
  dev.off()
  
  # End
})

##### Convergence of genes ######
# Get vector of comparisons
comparisons_to_make<-combn(1:32,2)

lapply(1:4,function(gen){
  converge_list<-lapply(1:length(levels(slim_dd$Treatment_Demo)),function(x){
    
    tmp<-slim_dd[slim_dd$Generation==paste0("Gen_",gen),]
    tmp<-tmp[tmp$Treatment_Demo == levels(tmp$Treatment_Demo)[x],]
    
    # Mark 95% outliers
    FST_out<-tmp[tmp$FST > quantile(na.omit(tmp$FST),probs=0.95),"Gene_ID"]
    dXY_out<-tmp[tmp$dXY > quantile(tmp$dXY,probs=0.95),"Gene_ID"]
    pi_out<-tmp[tmp$deltaPI > quantile(tmp$deltaPI,probs=0.95),"Gene_ID"]
    
    # Return
    return(list(FST_out,dXY_out,pi_out))
  })
  
  # Do all comparisons of overlap
  overlap<-lapply(1:length(comparisons_to_make[1,]),function(x){
    
    # Cols
    col1<-comparisons_to_make[1,x]
    col2<-comparisons_to_make[2,x]
    
    # Run comparisons over each measure for same measure eg. Fst - Fst and for all pairs eg Fst - Dxy
    tmp<-lapply(1:3,function(y){
      OL<-Reduce(intersect,list(converge_list[[col1]][[y]],
                                converge_list[[col2]][[y]]))
      return(length(OL))
    })
    
    return(tmp)  
  })
  
  FST_mat<-fill_matrix(unlist(lapply(overlap,function(x){return(x[[1]])})))
  DXY_mat<-fill_matrix(unlist(lapply(overlap,function(x){return(x[[2]])})))
  PI_mat<-fill_matrix(unlist(lapply(overlap,function(x){return(x[[3]])})))
 
  # Add Colnames and rownames
  colnames(FST_mat)<-levels(slim_dd$Treatment_Demo)
  rownames(FST_mat)<-levels(slim_dd$Treatment_Demo)
  colnames(DXY_mat)<-levels(slim_dd$Treatment_Demo)
  rownames(DXY_mat)<-levels(slim_dd$Treatment_Demo)
  colnames(PI_mat)<-levels(slim_dd$Treatment_Demo)
  rownames(PI_mat)<-levels(slim_dd$Treatment_Demo)
  
  # Plot
  prop_FST<-plot_prop_matrix(FST_mat)
  prop_DXY<-plot_prop_matrix(DXY_mat)
  prop_PI<-plot_prop_matrix(PI_mat)
  
  plots_to_colour<-list(prop_FST,prop_DXY,prop_PI)
  
  # Edit with cols
  colour_plots<-lapply(1:3,function(plot_in){
    
    # Set up colours for axes labels
    ghib_cols<-ghibli_palette("YesterdayMedium")[c(4,6,5,7)]
    prop_FST2<-plots_to_colour[[plot_in]]
    
    # Colour axes according to treatment
    ggbld <- ggplot_build(prop_FST2)
    xlabs<-ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.labels
    ylabs<-ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$y.labels
    
    pub_colours_y<-data.frame(labs=ylabs,
                                 mig_colour=as.character(rep(ghib_cols[1],31)),
                                 pop2_colour=rep(ghib_cols[1],31),
                                 bot_colour=rep(ghib_cols[1],31))
    
    pub_colours_x<-data.frame(labs=xlabs,
                                 mig_colour=rep(ghib_cols[1],31),
                                 pop2_colour=rep(ghib_cols[1],31),
                                 bot_colour=rep(ghib_cols[1],31))
    
    # Edit all
    # Mig
    pub_colours_y$mig_colour<-as.character(pub_colours_y$mig_colour)
    pub_colours_y$mig_colour[grep("Mig=0.002",ylabs)]<-ghib_cols[2]
    
    pub_colours_x$mig_colour<-as.character(pub_colours_x$mig_colour)
    pub_colours_x$mig_colour[grep("Mig=0.002",xlabs)]<-ghib_cols[2]
    
    # Contractions
    pub_colours_y$pop2_colour<-as.character(pub_colours_y$pop2_colour)
    pub_colours_y$pop2_colour[grep("Pop2=0.1",ylabs)]<-ghib_cols[2]
    pub_colours_y$pop2_colour[grep("Pop2=0.5",ylabs)]<-ghib_cols[3]
    pub_colours_y$pop2_colour[grep("Pop2=1",ylabs)]<-ghib_cols[4]
    
    pub_colours_x$pop2_colour<-as.character(pub_colours_x$pop2_colour)
    pub_colours_x$pop2_colour[grep("Pop2=0.1",xlabs)]<-ghib_cols[2]
    pub_colours_x$pop2_colour[grep("Pop2=0.5",xlabs)]<-ghib_cols[3]
    pub_colours_x$pop2_colour[grep("Pop2=1",xlabs)]<-ghib_cols[4]
    
    # Bots
    pub_colours_y$bot_colour<-as.character(pub_colours_y$bot_colour)
    pub_colours_y$bot_colour[grep("Bot=20",ylabs)]<-ghib_cols[2]
    pub_colours_y$bot_colour[grep("Bot=100",ylabs)]<-ghib_cols[3]
    pub_colours_y$bot_colour[grep("Bot=1000",ylabs)]<-ghib_cols[4]
    
    pub_colours_x$bot_colour<-as.character(pub_colours_x$bot_colour)
    pub_colours_x$bot_colour[grep("Bot=20",xlabs)]<-ghib_cols[2]
    pub_colours_x$bot_colour[grep("Bot=100",xlabs)]<-ghib_cols[3]
    pub_colours_x$bot_colour[grep("Bot=1000",xlabs)]<-ghib_cols[4]
    
    # Mig colours
    prop_FST2_mig<-prop_FST2+
      theme(axis.text.x = element_text(size=20,colour=pub_colours_x$mig_colour),
            axis.text.y = element_text(size=20,colour=pub_colours_y$mig_colour))
    
    return(prop_FST2_mig)
  })
  
  # Plot
  pdf(paste0("figs/slim_output_",mut_rate,"/Proportion_overlapping_genes_Gen",gen,".pdf"),width=14,height=14)
  print(prop_FST)
  print(prop_DXY)
  print(prop_PI)
  dev.off()
  
})

# ------------------------------------------
#### Downsampling Analysis #####
# ------------------------------------------


# Here we randomly downsample genes under strongest selection treatments to see if overlap is improved under more realistic assumptions
# For each iteration we return a matrix before averaging across matrices
iterations<-100
comparisons_to_make<-combn(1:32,2)


# Check if file exists
if(file.exists(paste0("outputs/slim_output_",mut_rate,"/downsampling_out",iterations,"_v2.rds"))=="FALSE"){
  
  downsampling_out<-lapply(1:iterations,function(iter){
    
    # Downsample slim randomly
    rand_slim<-random_downsample(slim_dd)
    
    converge_list<-lapply(1:length(levels(rand_slim$Treatment_Demo)),function(x){
      
      tmp<-rand_slim[rand_slim$Generation=="Gen_4",]
      tmp<-tmp[tmp$Treatment_Demo == levels(tmp$Treatment_Demo)[x],]
      
      # Mark 95% outliers
      FST_out<-tmp[tmp$FST > quantile(na.omit(tmp$FST),probs=0.95),"Gene_ID"]
      dXY_out<-tmp[tmp$dXY > quantile(tmp$dXY,probs=0.95),"Gene_ID"]
      pi_out<-tmp[tmp$deltaPI > quantile(tmp$deltaPI,probs=0.95),"Gene_ID"]
      
      # Return
      return(list(FST_out,dXY_out,pi_out))
    })
    
    # Return info on overlap across measures
    measure_overlap<-lapply(1:length(converge_list),function(treat){
      tmp<-calculate.overlap(converge_list[[treat]])
      return(c(length(tmp[[1]])/length(converge_list[[treat]][[1]]),
               length(tmp[[2]])/length(converge_list[[treat]][[1]]),
               length(tmp[[3]])/length(converge_list[[treat]][[1]]),
               length(tmp[[4]])/length(converge_list[[treat]][[1]])))
    })
    
    
    measure_overlap_mat<-matrix(ncol=4,nrow=32)
    for(i in 1:nrow(measure_overlap_mat)){
      measure_overlap_mat[i,]<-measure_overlap[[i]]
    }
    
    
    # Turn to dataframe
    measure_overlap_dd<-data.frame(overlap_all=measure_overlap_mat[,1],
                                   fst.dxy=measure_overlap_mat[,2],
                                   fst.pi=measure_overlap_mat[,3],
                                   dxy.pi=measure_overlap_mat[,4],
                                   treatment=levels(rand_slim$Treatment_Demo))
    
    # Also return add raw treatment groups
    for(i in 1:nrow(measure_overlap_dd)){
      measure_overlap_dd$Migration[i]<-slim_dd[slim_dd$Treatment_Demo == measure_overlap_dd$treatment[i],"Migration"][1]
      measure_overlap_dd$Contraction[i]<-slim_dd[slim_dd$Treatment_Demo == measure_overlap_dd$treatment[i],"Pop2_Size"][1]
      measure_overlap_dd$Bottleneck[i]<-slim_dd[slim_dd$Treatment_Demo == measure_overlap_dd$treatment[i],"Bottleneck"][1]
    }
    
    # Get vector of comparisons
    comparisons_to_make<-combn(1:32,2)
    
    # Do all comparisons of overlap
    overlap<-lapply(1:length(comparisons_to_make[1,]),function(x){
      
      # Cols
      col1<-comparisons_to_make[1,x]
      col2<-comparisons_to_make[2,x]
      
      comparisons_to_make2<-matrix(c(1,1,2,2,3,3,1,2,1,3,2,3),ncol=6)
      # Run comparisons over each measure for same measure eg. Fst - Fst and for all pairs eg Fst - Dxy
      tmp<-lapply(1:ncol(comparisons_to_make2),function(y){
        measure1<-comparisons_to_make2[1,y]
        measure2<-comparisons_to_make2[2,y]
        OL<-Reduce(intersect,list(converge_list[[col1]][[measure1]],
                                  converge_list[[col2]][[measure2]]))
        return(length(OL))
      })
      
      return(tmp)
    })
    
    FST_mat<-fill_matrix(unlist(lapply(overlap,function(x){return(x[[1]])})))
    DXY_mat<-fill_matrix(unlist(lapply(overlap,function(x){return(x[[2]])})))
    PI_mat<-fill_matrix(unlist(lapply(overlap,function(x){return(x[[3]])})))
    FST_DXY<-fill_matrix(unlist(lapply(overlap,function(x){return(x[[4]])})))
    FST_PI<-fill_matrix(unlist(lapply(overlap,function(x){return(x[[5]])})))
    DXY_PI<-fill_matrix(unlist(lapply(overlap,function(x){return(x[[6]])})))
    
    # Add Colnames and rownames
    colnames(FST_mat)<-levels(slim_dd$Treatment_Demo)
    rownames(FST_mat)<-levels(slim_dd$Treatment_Demo)
    colnames(DXY_mat)<-levels(slim_dd$Treatment_Demo)
    rownames(DXY_mat)<-levels(slim_dd$Treatment_Demo)
    colnames(PI_mat)<-levels(slim_dd$Treatment_Demo)
    rownames(PI_mat)<-levels(slim_dd$Treatment_Demo)
    
    colnames(FST_DXY)<-levels(slim_dd$Treatment_Demo)
    rownames(FST_DXY)<-levels(slim_dd$Treatment_Demo)
    colnames(FST_PI)<-levels(slim_dd$Treatment_Demo)
    rownames(FST_PI)<-levels(slim_dd$Treatment_Demo)
    colnames(DXY_PI)<-levels(slim_dd$Treatment_Demo)
    rownames(DXY_PI)<-levels(slim_dd$Treatment_Demo)
    
    
    return(list(FST_mat,DXY_mat,PI_mat,measure_overlap_dd,
                FST_DXY,FST_PI,DXY_PI))
  })
  
  # Save the outputs
  saveRDS(downsampling_out,paste0("outputs/slim_output_",mut_rate,"/downsampling_out",iterations,"_v2.rds"))
}

downsampling_out<-readRDS(paste0("outputs/slim_output_",mut_rate,"/downsampling_out",iterations,"_v2.rds"))

# Average across matrices
FST_list<-lapply(1:length(downsampling_out),function(x){return(downsampling_out[[x]][[1]])})
DXY_list<-lapply(1:length(downsampling_out),function(x){return(downsampling_out[[x]][[2]])})
PI_list<-lapply(1:length(downsampling_out),function(x){return(downsampling_out[[x]][[3]])})
measure_overlap_list<-lapply(1:length(downsampling_out),function(x){return(downsampling_out[[x]][[4]][,1:4])})

FSTxDXY_list<-lapply(1:length(downsampling_out),function(x){return(downsampling_out[[x]][[5]])})
FSTxPI_list<-lapply(1:length(downsampling_out),function(x){return(downsampling_out[[x]][[6]])})
DXYxPI_list<-lapply(1:length(downsampling_out),function(x){return(downsampling_out[[x]][[7]])})


# Reduce
FST_avg<-Reduce("+", FST_list) / length(FST_list)
DXY_avg<-Reduce("+", DXY_list) / length(DXY_list)
PI_avg<-Reduce("+", PI_list) / length(PI_list)
measure_overlap_avg<-Reduce("+", measure_overlap_list) / length(measure_overlap_list)

FSTxDXY_avg<-Reduce("+", FSTxDXY_list) / length(FSTxDXY_list)
FSTxPI_avg<-Reduce("+", FSTxPI_list) / length(FSTxPI_list)
DXYxPI_avg<-Reduce("+", DXYxPI_list) / length(DXYxPI_list)


# Plot as above
# Plot
prop_FST<-plot_prop_matrix(FST_avg)
prop_DXY<-plot_prop_matrix(DXY_avg)
prop_PI<-plot_prop_matrix(PI_avg)
prop_FSTxDXY<-plot_prop_matrix(FSTxDXY_avg)
prop_FSTxPI<-plot_prop_matrix(FSTxPI_avg)
prop_DXYxPI<-plot_prop_matrix(DXYxPI_avg)
plots_to_colour<-list(prop_FST,prop_DXY,prop_PI,
                      prop_FSTxDXY,prop_FSTxPI,prop_DXYxPI)


# Edit with cols
colour_plots<-lapply(1:6,function(plot_in){
  
  # Set up colours for axes labels
  ghib_cols<-ghibli_palette("YesterdayMedium")[c(4,6,5,7)]
  prop_FST2<-plots_to_colour[[plot_in]]
  
  # Colour axes according to treatment
  ggbld <- ggplot_build(prop_FST2)
  xlabs<-ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.labels
  ylabs<-ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$y.labels
  
  pub_colours_y<-data.frame(labs=ylabs,
                               mig_colour=as.character(rep(ghib_cols[1],31)),
                               pop2_colour=rep(ghib_cols[1],31),
                               bot_colour=rep(ghib_cols[1],31))
  
  pub_colours_x<-data.frame(labs=xlabs,
                               mig_colour=rep(ghib_cols[1],31),
                               pop2_colour=rep(ghib_cols[1],31),
                               bot_colour=rep(ghib_cols[1],31))
  
  # Edit all
  # Mig
  pub_colours_y$mig_colour<-as.character(pub_colours_y$mig_colour)
  pub_colours_y$mig_colour[grep("Mig=0.002",ylabs)]<-ghib_cols[2]
  
  pub_colours_x$mig_colour<-as.character(pub_colours_x$mig_colour)
  pub_colours_x$mig_colour[grep("Mig=0.002",xlabs)]<-ghib_cols[2]
  
  # Contractions
  pub_colours_y$pop2_colour<-as.character(pub_colours_y$pop2_colour)
  pub_colours_y$pop2_colour[grep("Pop2=0.1",ylabs)]<-ghib_cols[2]
  pub_colours_y$pop2_colour[grep("Pop2=0.5",ylabs)]<-ghib_cols[3]
  pub_colours_y$pop2_colour[grep("Pop2=1",ylabs)]<-ghib_cols[4]
  
  pub_colours_x$pop2_colour<-as.character(pub_colours_x$pop2_colour)
  pub_colours_x$pop2_colour[grep("Pop2=0.1",xlabs)]<-ghib_cols[2]
  pub_colours_x$pop2_colour[grep("Pop2=0.5",xlabs)]<-ghib_cols[3]
  pub_colours_x$pop2_colour[grep("Pop2=1",xlabs)]<-ghib_cols[4]
  
  # Bots
  pub_colours_y$bot_colour<-as.character(pub_colours_y$bot_colour)
  pub_colours_y$bot_colour[grep("Bot=20",ylabs)]<-ghib_cols[2]
  pub_colours_y$bot_colour[grep("Bot=100",ylabs)]<-ghib_cols[3]
  pub_colours_y$bot_colour[grep("Bot=1000",ylabs)]<-ghib_cols[4]
  
  pub_colours_x$bot_colour<-as.character(pub_colours_x$bot_colour)
  pub_colours_x$bot_colour[grep("Bot=20",xlabs)]<-ghib_cols[2]
  pub_colours_x$bot_colour[grep("Bot=100",xlabs)]<-ghib_cols[3]
  pub_colours_x$bot_colour[grep("Bot=1000",xlabs)]<-ghib_cols[4]
  
  # Mig colours
  prop_FST2_mig<-prop_FST2+
    theme(axis.text.x = element_text(size=20,colour=pub_colours_x$mig_colour),
          axis.text.y = element_text(size=20,colour=pub_colours_y$mig_colour))
  
  return(prop_FST2_mig)
})

for(i in 1:6){
colour_plots[[i]]<-remove_geom(colour_plots[[i]],"GeomText")
}

# Plot
pdf(paste0("figs/slim_output_",mut_rate,"/Downsampled_Proportion_overlapping_genes.pdf"),width=14,height=14)
for(i in 1:6){
print(colour_plots[[i]])
}
dev.off()

# Return measure_overlap to dd form
measure_overlap_final<-cbind(measure_overlap_avg,downsampling_out[[1]][[4]][,5:length(downsampling_out[[1]][[4]])])

# Plot figure for each treatment
treatment_vec2<-c("Bottleneck","Contraction","Migration")

plot_list<-lapply(1:length(treatment_vec2),function(x){
  
  tmp<-data.frame(overlap=unlist(list(measure_overlap_final[,1:4])),
                  treatment=rep(measure_overlap_final[,treatment_vec2[x]],4),
                  measures=rep(c("All","fst.dxy","fst.pi","dxy.pi"),each=32))
  
  tmp$measures<-factor(tmp$measures,levels=c("All","fst.dxy","fst.pi","dxy.pi"))
  levels(tmp$measures)<-c("All","F[ST]~x~D[XY]","F[ST]~x~Delta*pi","D[XY]~x~Delta*pi")
  
  tmp$treatment<-as.factor(tmp$treatment)
  
  g1<-ggplot(tmp,aes(x=treatment,y=overlap,fill=measures))+
    #geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, position=pd,show.legend = FALSE) +
    #geom_point(position=pd,size=3)+
    geom_boxplot()+
    #facet_wrap(~measures,ncol=1,scales = "free_x",labeller = label_parsed)+
    theme_bw()+
    theme(axis.title= element_text(size=18),
          axis.text=element_text(size=16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          strip.text = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "right")+
    xlab(treatment_vec2[x])+
    scale_fill_manual(breaks=levels(tmp$measures),
                      labels=c("All",
                               expression(F[ST]~x~D[XY]),
                               expression(F[ST]~x~Delta*pi),
                               expression(D[XY]~x~Delta*pi)),
                      values = rev(ghibli_palette("YesterdayMedium"))) +
    labs(fill="Measures")+
    ylab("Overlap Proportion")
  
  return(g1)
})

pdf(paste0("figs/slim_output_",mut_rate,"/Overlap_proportions_top_outliers_across_measures_Gen4.pdf"),height=10,width=6)
print(ggarrange(plot_list[[1]],
          plot_list[[2]],
          plot_list[[3]],
          ncol=1,nrow=3,common.legend = T,legend="right",labels = "AUTO"))
dev.off()

# --------------------------
##### Effect of Demography on Evolv Gen #####
# --------------------------

treatment_vec2<-c("Bottleneck","Contraction","Migration")

evolv_gen_list<-lapply(1:length(treatment_vec2),function(x){
  
  tmp<-as.data.frame(slim_dd %>% 
                       group_by(slim_dd[,treatment_vec[x]]) %>%
                       dplyr::summarise(mean=mean(Evolv_Gen[Evolv_Gen!= 0]),
                                        sd=sd(Evolv_Gen[Evolv_Gen!= 0]),
                                        se=(sd(Evolv_Gen[Evolv_Gen!= 0])/sqrt(length(Evolv_Gen[Evolv_Gen!= 0])))))
  
  
  tmp2<-as.data.frame(slim_dd %>% 
                        group_by(slim_dd[,treatment_vec[x]],Generation) %>%
                        dplyr::summarise(prop_survive=((length(Evolv_Gen[Evolv_Gen!= 0]))/length(Evolv_Gen))))
  
  plot_dd<-slim_dd[slim_dd$Evolv_Gen !=0,]
  
  g1<-ggplot(data = plot_dd, aes(x = as.factor(plot_dd[,treatment_vec[x]]), y = Evolv_Gen,fill=as.factor(plot_dd[,treatment_vec[x]]))) +
    geom_violin(draw_quantiles=0.5,width=1)+
    #geom_boxplot(outlier.colour=NA)+
    stat_summary(fun.y = "mean",geom="point")+
    theme_bw()+
    theme(axis.text.x=element_text(size=24,angle=45,hjust=1),
          axis.text.y=element_text(size=20),
          axis.title.y=element_text(size=22),
          panel.grid=element_blank(),
          legend.title=element_text(size=20),
          legend.text=element_text(size=18),
          legend.position = "none")+
    scale_fill_manual(values=rev(ghibli_palette("YesterdayMedium")))+
    labs(x=treatment_vec[x],y="Evolution Generation")
  
  return(list(tmp,tmp2,g1))
})

# Plot them and save
pdf(paste0("figs/slim_output_",mut_rate,"/Evolv_gen_by_demography.pdf"),height=4,width=14)
print(ggarrange(evolv_gen_list[[1]][[3]],
          evolv_gen_list[[2]][[3]],
          evolv_gen_list[[3]][[3]],ncol=3))
dev.off()


# ------------------------------------------
#### Effect of Demography on Phenotype #####
# ------------------------------------------

treatment_vec<-c("Bottleneck","Pop2_Size","Migration")
treatment_vec2<-c("Bottleneck","Contraction","Migration")

pheno_g_list<-lapply(1:4,function(gen){

# Only run for last gen
slim_dd_GEN4<-slim_dd[slim_dd$Generation==paste0("Gen_",gen),]

# Plot the distribution of phenotypes in a facet of demographic treatment and selection
pheno_dd<-data.frame(pheno=rep(slim_dd_GEN4$Pop2_Mean,3),
                     pheno_sd=rep(slim_dd_GEN4$Pop2_sd,3),
                     treatment=c(slim_dd_GEN4$Bottleneck,
                                 slim_dd_GEN4$Pop2_Size,
                                 slim_dd_GEN4$Migration),
                     treatment_facet=rep(treatment_vec2,each=length(slim_dd_GEN4$Pop2_Mean)),
                     selection=rep(slim_dd_GEN4$Selection_Coef,3))

pheno_dd$treatment<-factor(pheno_dd$treatment,
                           levels = levels(as.factor(pheno_dd$treatment)))
# Plot
pheno_g<-ggplot(pheno_dd,aes(x=treatment,y=pheno))+
  facet_grid(selection~treatment_facet,scales="free_x")+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  geom_hline(yintercept = 10,linetype="dashed")+
  theme_bw()+
  theme(axis.title= element_text(size=24),
        axis.text=element_text(size=20),
        legend.text = element_text(size=22),
        legend.title = element_text(size=22),
        strip.text = element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right")+
  xlab("Treatment")+
  ylab("Phenotype")

return(pheno_g)
})

pdf(paste0("figs/slim_output_",mut_rate,"/Phenotypes_by_demography.pdf"),height=12,width=12)
for(j in 1:4){
print(pheno_g_list[[j]])
}
dev.off()

# Also Plot Evolv_Gen as above
# Remove 0 Evolv_Gen
slim_dd_GEN4_EvolvGen<-slim_dd_GEN4[slim_dd_GEN4$Evolv_Gen != 0,]

evolvgen_dd<-data.frame(evolvgen=rep(slim_dd_GEN4_EvolvGen$Evolv_Gen,3),
                        treatment=c(slim_dd_GEN4_EvolvGen$Bottleneck,
                                    slim_dd_GEN4_EvolvGen$Pop2_Size,
                                    slim_dd_GEN4_EvolvGen$Migration),
                        treatment_facet=rep(treatment_vec2,each=length(slim_dd_GEN4_EvolvGen$Pop2_Mean)),
                        selection=rep(slim_dd_GEN4_EvolvGen$Selection_Coef,3))

evolvgen_dd$treatment<-factor(evolvgen_dd$treatment,
                              levels = levels(as.factor(evolvgen_dd$treatment)))
# Plot
evolvgen_g<-ggplot(evolvgen_dd,aes(x=treatment,y=evolvgen))+
  facet_grid(selection~treatment_facet,scales="free_x")+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  theme_bw()+
  theme(axis.title= element_text(size=18),
        axis.text=element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        strip.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right")+
  scale_fill_manual(values = rev(ghibli_palette("YesterdayMedium"))) +
  xlab("Treatment")+
  ylab("Generation Optimum Reached")

pdf(paste0("figs/slim_output_",mut_rate,"/Evolv_gen_by_demography_v2.pdf"),height=20,width=20)
print(evolvgen_g)
dev.off()

##### Generation Post Adaptation (GPA) #####

# Make a new dataset for modding
GPA_dd<-slim_dd
gens<-c(0.05*2*pop_size,
        0.158*2*pop_size,
        1.581*2*pop_size,
        5*2*pop_size)

# Fix generation to number
for(i in 1:4){
  GPA_dd[GPA_dd$Generation == paste0("Gen_",i),"Generation"]<-gens[i]
}

# Remove all incidences where phenotype was not reached
GPA_dd<-GPA_dd[GPA_dd$Evolv_Gen != 0,]
# Remove all incidences where measure was taken before phenotype was reached (+ 50 because of avg)
GPA_dd<-GPA_dd[GPA_dd$Evolv_Gen < (as.numeric(GPA_dd$Generation)+49),]

# Calculate GPA
GPA_dd$GPA<-as.numeric(GPA_dd$Generation)-GPA_dd$Evolv_Gen

# How does relationship over time vary demographically?
# Run over all treatments
plot_list<-lapply(1:length(treatment_vec),function(x){
  
  plot_dd<-data.frame(treatment=rep(GPA_dd[,treatment_vec[x]],3),
                      divergence=c(GPA_dd$FST,GPA_dd$dXY,GPA_dd$deltaPI),
                      selection=rep(GPA_dd$Selection_Coef,3),
                      GPA=rep(GPA_dd$GPA,3),
                      measure=rep(c("F[ST]","D[XY]","Delta*pi"),each=length(GPA_dd$FST)))
  
  plot_dd$selection<-as.factor(plot_dd$selection)
  plot_dd$treatment<-as.factor(plot_dd$treatment)
  plot_dd$measure2<-factor(plot_dd$measure,levels=c("F[ST]","D[XY]","Delta*pi"))
  
  g1<-ggplot(plot_dd,aes(x=GPA,y=divergence,colour=treatment))+
    geom_smooth(method="lm")+
    facet_grid(measure2~selection,labeller = label_parsed,scales = "free_y")+
    scale_colour_manual(values = rev(ghibli_palette("YesterdayMedium"))) +
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.title=element_text(size=20),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=16,angle=45,hjust=1),
          strip.text = element_text(size=20),
          legend.text = element_text(size=14),
          legend.title = element_text(size=15))+
    labs(y="Divergence",x="Generations Post Adaptation",colour="Treatment")
  
  return(g1)
})

pdf(paste0("figs/slim_output_",mut_rate,"/GPA_figures.pdf"),width=16,height=10)
print(plot_list[[1]])
print(plot_list[[2]])
print(plot_list[[3]])
dev.off()

##### What explains measures? #####

# Partition data by generation and all 32 unique treatments
treats<-unique(as.character(slim_dd$Treatment_Demo))

# Run over gens
gen_out<-data.frame(rbindlist(pblapply(1:4,function(x){
  gen_tmp<-slim_dd[slim_dd$Generation == paste0("Gen_",x),]
  
  # Run over within treatments
  treat_out<-data.frame(rbindlist(lapply(1:length(treats),function(y){
    
    treat_tmp<-gen_tmp[gen_tmp$Treatment_Demo == treats[y],]
    
    # What explains variance for each measure?
    measures<-c("FST","dXY","deltaPI")
    measures_out<-data.frame(rbindlist(lapply(1:length(measures),function(z){
      
      # Make col
      treat_tmp$measure<-treat_tmp[,measures[z]]
      
      # Fit forest over 100 trees
      fitregforest<-randomForest(measure~Gene_Size+Target_sel+Selection_Coef+Evolv_Gen+Exon_N+Pop2_PI,
                                 data=treat_tmp,
                                 ntree=100)
      
      # Get importance of variables
      out_mat<-matrix(nrow=1,ncol=length(t(fitregforest$importance)))
      props<-fitregforest$importance/sum(fitregforest$importance)
      out_mat[1,]<-(t(props))
      out_dd<-data.frame(out_mat)
      colnames(out_dd)<-rownames(fitregforest$importance)
      
      # Tidy dataframe
      out_dd$Generation<-paste0("Gen_",x)
      out_dd$Treatment<-treats[y]
      out_dd$Measure<-measures[z]
      
      return(out_dd)
      
    })))
    
    return(measures_out)
  })))
  return(treat_out)
})))

# Tidy Gen_out by incorporating raw values of treatment levels
for(i in 1:length(gen_out[,1])){
  gen_out$Migration[i]<-slim_dd[slim_dd$Treatment_Demo == gen_out$Treatment[i],"Migration"][1]
  gen_out$Contraction[i]<-slim_dd[slim_dd$Treatment_Demo == gen_out$Treatment[i],"Pop2_Size"][1]
  gen_out$Bottleneck[i]<-slim_dd[slim_dd$Treatment_Demo == gen_out$Treatment[i],"Bottleneck"][1]
}

# Save output
write.table(gen_out,
            paste0("outputs/slim_output_",mut_rate,"/randomForest_results_allgenerations_alltreatments.txt"),
            quote = F, row.names = F,sep="\t")

# What best explains FST in Gen4 heatmaps
gen4<-gen_out[gen_out$Generation == "Gen_4",]
gen4_2<-data.frame(importance=unlist(list(gen4[,1:6])),
                   variable=rep(colnames(gen4)[1:6],each=length(gen4$Gene_Size)),
                   Migration=rep(gen4$Migration,6),
                   Contraction=rep(gen4$Contraction,6),
                   measure=rep(gen4$Measure,6))

gen4_2$Migration<-as.factor(gen4_2$Migration)
gen4_2$Contraction<-as.factor(gen4_2$Contraction)
gen4_2$Metric<-factor(gen4_2$measure,levels=c("FST","dXY","deltaPI"))
levels(gen4_2$Metric)<-c("F[ST]","D[XY]","Delta*pi")

gen4_2$variable<-as.character(gen4_2$variable)

# Relabel variables
gen4_2[gen4_2$variable == "Evolv_Gen","variable"]<-"Pheno Gen"
gen4_2[gen4_2$variable == "Target_sel","variable"]<-"Target Proportion"
gen4_2[gen4_2$variable == "Gene_Size","variable"]<-"Gene Size"
gen4_2[gen4_2$variable == "Selection_Coef","variable"]<-"Selection"
gen4_2[gen4_2$variable == "Exon_N","variable"]<-"Exon N"
gen4_2[gen4_2$variable == "Pop2_PI","variable"]<-"Diverging Pi"

# Plot
mig_fig<-ggplot(gen4_2,aes(x=Migration,y=importance,fill=variable))+
  geom_boxplot(outlier.alpha = 0)+
  facet_wrap(~Metric,ncol=1,scales = "free",labeller = label_parsed)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=24),
        axis.text = element_text(size=22),
        strip.text = element_text(size=22),
        legend.title = element_text(size=20),
        legend.text = element_text(size=14))+
  ylab("Relative Importance")+
  labs(fill="Source")+
  scale_fill_manual(values = ghibli_palette("PonyoMedium")[c(4,6,7,5,3,2)])

contraction_fig<-ggplot(gen4_2,aes(x=Contraction,y=importance,fill=variable))+
  geom_boxplot(outlier.alpha = 0)+
  facet_wrap(~Metric,ncol=1,scales = "free",labeller = label_parsed)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=24),
        axis.text = element_text(size=22),
        strip.text = element_text(size=22),
        legend.title = element_text(size=20),
        legend.text = element_text(size=14))+
  ylab("Relative Importance")+
  labs(fill="Source")+
  scale_fill_manual(values = ghibli_palette("PonyoMedium")[c(4,6,7,5,3,2)])

# Export fig
pdf(paste0("figs/slim_output_",mut_rate,"/Migration_effect_on_heatmap.pdf"),
    width=8,height=12)
print(mig_fig)
print(contraction_fig)
dev.off()


}


get_blups<-function(fixedFormula,randFormula,rcovFormula,METdata,
                    removeOutliers=TRUE,...){
     # Run model first time
     mmer_output<-mmer(fixed = as.formula(fixedFormula),
                       random = as.formula(randFormula),
                       rcov= as.formula(rcovFormula),
                       data=METdata, 
                       getPEV = T)
     if(removeOutliers==TRUE){
          # Outlier Detection and Removal
          ## index observations that are defined as outliers
          outliers<-which(abs(scale(mmer_output$residuals))>3.3)
          
          # remove outliers, if any
          if(length(outliers)>0){ 
               x<-METdata[-outliers,] 
               # Refit the model
               starttime<-proc.time()[3]
               mmer_output<-mmer(fixed = as.formula(fixedFormula),
                                 random = as.formula(randFormula),
                                 rcov= as.formula(rcovFormula),
                                 data=x, 
                                 getPEV = T)
               stoptime<-proc.time()[3]; elapsed<-stoptime-starttime; elapsed/60
          }
     }
     if(length(outliers)==0 | removeOutliers==FALSE){ outliers<-NULL }
     
     # log likelihood of the model, AIC, convergence T/F
     modelfit<-summary(mmer_output)$logo
     # number of groups for factors, could be used to compute DF 
     groups<-summary(mmer_output)$groups
     # variance components
     varcomp<-summary(mmer_output)$varcomp
     
     # genotypic variance
     Vg<-mmer_output$sigma$germplasmName %>% as.vector()
     # Mean error-variance across trials 
     ## across heterogeneous error variance estimates
     meanHetVarE<-mmer_output$sigma %>% 
          .[names(.) %>% grep(":units",.,value=T)] %>% 
          unlist() %>% 
          mean()
     # Mean number of reps per accession
     meanNreps<-METdata %>%
          count(germplasmName) %$% mean(n)
     # Broad-sense heritability
     H2<-Vg/(Vg+(meanHetVarE/meanNreps))
     
     # Extract the BLUPs and PEVs, compute Reliability (REL), 
     # de-regressed BLUPs and weights for downstream analysis
     blups<-mmer_output$U$germplasmName$Value %>% 
          tibble(germplasmName=names(.),BLUP=.) %>% 
          mutate(germplasmName=gsub("germplasmName","",germplasmName),
                 PEV=diag(mmer_output$PevU$germplasmName$Value), # prediction error variance
                 REL=1-PEV/Vg, # Reliability
                 drgBLUP=BLUP/REL, # De-regressed BLUP
                 WT=(1-H2)/((0.1 + (1-REL)/REL)*H2) # weight for downstream
          )
     
     # Combine all outputs into one object the function can return() 
     out<-tibble(H2=H2,meanHetVarE=meanHetVarE,meanNreps=meanNreps,
                 modelfit=list(modelfit),
                 groups=list(groups),
                 blups=list(blups),
                 varcomp=list(varcomp),
                 outliers=list(outliers))
     return(out)
}
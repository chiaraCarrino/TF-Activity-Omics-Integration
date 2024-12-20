library (glmnet)
library(stringr)
library(caret)

###### preprocessing of TN FoldChange of phosphorylations values
df=read.table("input/Phospho_TN_FC.txt", sep=',', header=TRUE)

df_num = as.matrix(df[,2:ncol(df)])
rownames(df_num) = sapply(df$Gene_Phosphosite,function(x) strsplit(as.character(x),split = "\\\\")[[1]][1])

df_nes1=read.table("input/TF_Activity_sum_pos_sub_neg.tsv", sep='\t', header=TRUE)
                          
new_vet_paz<-c()

for (paz in df$Gene_Phosphosite){
  
 
  paz_splitted<-str_split(paz, '_')[[1]]
  
  new_paz<-(paste(paz_splitted[1], '-', paz_splitted[2], sep=''))
  new_vet_paz=append(new_vet_paz, new_paz)
}

df$case_submitter_id=new_vet_paz
                          
df_nes=df_nes1[df_nes1$case_submitter_id %in% c(df[['case_submitter_id']]), ]                          
                          
df_num_nes = as.matrix(df_nes[,2:ncol(df_nes)])
rownames(df_num_nes) = sapply(df_nes$case_submitter_id,function(x) strsplit(as.character(x),split = "\\\\")[[1]][1])

                              
#This dataset contains the transcription factors (TFs) and the corresponding phosphorylations whose substrates are involved in at least one pathway associated with the TFs.
pathway_df=read.csv('../Phosphorylations_activity_association/file_output/MSigDB_TF_gene_phosphorylations_pathways_TN_FC.txt', sep='\t')

tf_list<-c(colnames(df_num_nes))
                             
for (var_tf in tf_list){
  
  print(var_tf)

  var_tf_simple=c(str_split_i(var_tf, '_',1))
  
    
    
  phos_2_select_all=unlist(pathway_df[pathway_df['TF'] == var_tf_simple, ]['Phosphorylation'])

  phos_2_select=c()
  if (":" %in% colnames(df_num)[1]){
        
     for (e in phos_2_select_all){
          if ((grepl(':', e))){
             phos_2_select=append(phos_2_select,e)
            }

        }
  }else if (!(":" %in% colnames(df_num)[1])){
        
        for (e in phos_2_select_all){
            if (!(grepl(':', e))){
                phos_2_select=append(phos_2_select,e)
            }

        }
    }
    phos_2_select=unique(phos_2_select)
    
    print(phos_2_select)
    
  
    df=read.table("input/Phospho_TN_FC.txt", sep=',', header=TRUE)

    df_num = as.matrix(df[,phos_2_select])
    rownames(df_num) = sapply(df$case_submitter_id,function(x) strsplit(as.character(x),split = "\\\\")[[1]][1])
    tf_df=data.frame(df_num)
    tf_df$y=as.data.frame(df_num_nes)[[var_tf]]                         
                              
  if (ncol(tf_df)>1){  
 
   
  set.seed(1985)
            
  #creating training and test set
  index<-createDataPartition(tf_df$y, p=.8, list=FALSE, times=1)
  
  train_df<-tf_df[index,]

  test_df<-tf_df[-index,]
        
  #cross-validation on training set
  ctrlspecs<- trainControl(method="cv", number=5,
                           savePredictions='all')
  #lamda-vector creation
  lambda_vector<-10^seq(5,-5, length=500)
   
  #creastion of a LASSO model for the TF   
  model1<- train(#y ~ .,
                df_num, tf_df$y,
                form=tf_df$y ~ df_num,
                 data=train_df,
                 preProcess=c("center","scale"),
                 method="glmnet",
                 tuneGrid=expand.grid(alpha=1, lambda=lambda_vector),
                 trControl=ctrlspecs,
                 na.action=na.omit)
  
  #selection of the labda best value
  model1$bestTune
  
  #selection of the coefficients identified by LASSO  
  df_coef <- round(coef(model1$finalModel, model1$bestTune$lambda), 3)
  
  if (length(df_coef[df_coef[, 1] != 0, ])>1){
    
    #test of the LASSO model on the test set
    predictions1<-predict(model1, newdata=test_df)
     
    RMSE<-RMSE(predictions1, test_df$y)
    Rsquared<-R2(predictions1, test_df$y)
    mod1perf<-data.frame(RMSE,
                         Rsquared)
    
    mod1perf['TF']<-var_tf
    
      write.table(mod1perf,'output/LASSO_regression/RMSE_Rsquared_phospho_TN_FC_TF_onlyDBCK_NES_SumPosSubNeg_LassoOnPathwayPhosphorylations.txt', append=TRUE, sep='\t', row.names=FALSE, col.names=FALSE)
      
    
    df_coef_selected<-as.data.frame(df_coef[df_coef[, 1] != 0, ])
    colnames(df_coef_selected)[which(names(df_coef_selected) == "df_coef[df_coef[, 1] != 0, ]")] <- "Contribution"
    df_coef_selected$TF<-var_tf
    library(tibble)
    df_coef_selected <- tibble::rownames_to_column(df_coef_selected, "Site")
       write.table(df_coef_selected,'output/LASSO_regression/Coef_phospho_TN_FC_onlyDBCK_NES_SumPosSubNeg_LassoOnPathwayPhosphorylations.txt', append=TRUE, sep='\t', row.names=FALSE, col.names=FALSE)
       
      }
    }
}


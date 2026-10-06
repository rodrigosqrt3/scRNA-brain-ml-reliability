# ==============================================================================
# ERAMIA-RS 2026 - REPRODUCIBLE SCRNA-SEQ CLASSIFICATION BENCHMARK
# ==============================================================================
# Primary analysis: repeated stratified cross-validation with paired resamples.
# Gene filtering, variable-gene selection, scaling and PCA are fitted exclusively
# inside each training fold. Results are checkpointed after every model fit.
#
# Run with RStudio Source or: Rscript Projeto.R
# ==============================================================================

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
master_seed<-335705
n_folds<-5
n_repeats<-5
min_detected_cells<-10
n_hvg<-2000
n_pcs<-100
library_scale<-10000
detected_cores<-parallel::detectCores(logical=FALSE)
if(is.na(detected_cores))detected_cores<-1
n_threads<-max(1,detected_cores-1)
models_to_run<-c("Random Forest","SVM","XGBoost","Naive Bayes","KNN")
representations_to_run<-c("Genes","100 PCs")
resume_from_checkpoints<-TRUE
save_top_features<-50

# Prespecified settings are held constant across all folds. No model is given
# an unequal tuning budget, and validation performance is not used for tuning.
rf_ntree<-500
svm_cost<-1
xgb_nrounds<-100
xgb_max_depth<-6
xgb_eta<-0.3
xgb_subsample<-1
xgb_colsample_bytree<-1
knn_k<-5

get_script_dir<-function(){
  args<-commandArgs(trailingOnly=FALSE)
  file_arg<-grep("^--file=",args,value=TRUE)
  if(length(file_arg)>0)return(dirname(normalizePath(sub("^--file=","",file_arg[1]))))
  frame_files<-vapply(sys.frames(),function(x)if(is.null(x$ofile))"" else x$ofile,character(1))
  frame_files<-frame_files[nzchar(frame_files)]
  if(length(frame_files)>0)return(dirname(normalizePath(frame_files[length(frame_files)])))
  getwd()
}

project_dir<-get_script_dir()
source(file.path(project_dir,"pipeline_functions.R"))
analysis_id<-paste0("repeated",n_repeats,"x",n_folds,"_hvg",n_hvg,"_pc",n_pcs)
results_dir<-file.path(project_dir,"results",analysis_id)
checkpoint_dir<-file.path(results_dir,"checkpoints")
table_dir<-file.path(results_dir,"tables")
figure_dir<-file.path(results_dir,"figures")
dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)

required_packages<-c("scRNAseq","SummarizedExperiment","Matrix","matrixStats","randomForest","e1071","xgboost","caret","pROC","irlba","ggplot2","dplyr","tidyr")
missing_packages<-required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)>0)stop(paste("Missing packages:",paste(missing_packages,collapse=", ")))
`%>%`<-dplyr::`%>%`
options(stringsAsFactors=FALSE)
set.seed(master_seed)

settings<-list(rf_ntree=rf_ntree,svm_cost=svm_cost,xgb_nrounds=xgb_nrounds,xgb_max_depth=xgb_max_depth,xgb_eta=xgb_eta,xgb_subsample=xgb_subsample,xgb_colsample_bytree=xgb_colsample_bytree,knn_k=knn_k,save_top_features=save_top_features)
configuration<-data.frame(parameter=c("master_seed","n_folds","n_repeats","min_detected_cells","n_hvg","n_pcs","library_scale","rf_ntree","svm_cost","xgb_nrounds","xgb_max_depth","xgb_eta","xgb_subsample","xgb_colsample_bytree","knn_k","n_threads"),value=as.character(c(master_seed,n_folds,n_repeats,min_detected_cells,n_hvg,n_pcs,library_scale,rf_ntree,svm_cost,xgb_nrounds,xgb_max_depth,xgb_eta,xgb_subsample,xgb_colsample_bytree,knn_k,n_threads)))
configuration_path<-file.path(table_dir,"analysis_configuration.csv")
if(file.exists(configuration_path)&&length(list.files(checkpoint_dir,pattern="\\.rds$"))>0){
  previous_configuration<-read.csv(configuration_path,stringsAsFactors=FALSE)
  if(!identical(paste(previous_configuration$parameter,previous_configuration$value),paste(configuration$parameter,configuration$value)))stop("Configuration changed while checkpoints exist. Restore the previous settings or use a new results directory.")
}else{
  write.csv(configuration,configuration_path,row.names=FALSE)
}

cat("\n==============================================================================\n")
cat("ERAMIA-RS 2026 - REPEATED CROSS-VALIDATION\n")
cat("==============================================================================\n")
cat("Analysis ID:",analysis_id,"\n")
cat("Output directory:",results_dir,"\n")
cat("Resampling:",n_repeats,"repeats x",n_folds,"folds\n")
cat("Models:",paste(models_to_run,collapse=", "),"\n\n")

# ==============================================================================
# 2. DATA IMPORT AND SAMPLE-WISE NORMALIZATION
# ==============================================================================
cat("Loading Zeisel mouse brain data...\n")
sce<-scRNAseq::ZeiselBrainData()
counts_raw<-SummarizedExperiment::assay(sce,"counts")
if(!inherits(counts_raw,"sparseMatrix"))counts_raw<-Matrix::Matrix(counts_raw,sparse=TRUE)
cell_labels<-factor(sce$level1class)
valid_cells<-!is.na(cell_labels)&Matrix::colSums(counts_raw)>0
counts_raw<-counts_raw[,valid_cells,drop=FALSE]
cell_labels<-droplevels(cell_labels[valid_cells])
cell_ids<-colnames(counts_raw)
if(is.null(cell_ids))cell_ids<-sprintf("cell_%04d",seq_len(ncol(counts_raw)))

# This transformation uses only each cell's own library size. All operations
# learned across cells remain inside training folds in prepare_fold_data().
library_sizes<-Matrix::colSums(counts_raw)
logcounts<-counts_raw%*%Matrix::Diagonal(x=library_scale/library_sizes)
logcounts<-log1p(logcounts)

data_overview<-data.frame(n_cells=ncol(counts_raw),n_genes=nrow(counts_raw),n_classes=nlevels(cell_labels),median_library_size=median(library_sizes),median_detected_genes=median(Matrix::colSums(counts_raw>0)),global_sparsity=1-Matrix::nnzero(counts_raw)/(nrow(counts_raw)*ncol(counts_raw)))
class_distribution<-data.frame(class=levels(cell_labels),n=as.numeric(table(cell_labels)),proportion=as.numeric(prop.table(table(cell_labels))))
write.csv(data_overview,file.path(table_dir,"data_overview.csv"),row.names=FALSE)
write.csv(class_distribution,file.path(table_dir,"class_distribution.csv"),row.names=FALSE)
cat("Cells:",ncol(counts_raw),"| Genes:",nrow(counts_raw),"| Classes:",nlevels(cell_labels),"\n")
print(class_distribution)

# ==============================================================================
# 3. REPEATED STRATIFIED CROSS-VALIDATION
# ==============================================================================
fold_map<-create_repeated_stratified_folds(cell_labels,n_folds,n_repeats,master_seed)
fold_map$cell_id<-cell_ids[fold_map$cell_index]
fold_map$class<-as.character(cell_labels[fold_map$cell_index])
write.csv(fold_map,file.path(table_dir,"resampling_folds.csv"),row.names=FALSE)

class_levels<-levels(cell_labels)
probability_columns<-paste0("prob__",make.names(class_levels,unique=TRUE))
expected_tasks<-n_repeats*n_folds*length(representations_to_run)*length(models_to_run)
completed_tasks<-0

for(repeat_id in seq_len(n_repeats)){
  for(fold_id in seq_len(n_folds)){
    test_index<-fold_map$cell_index[fold_map$repeat_id==repeat_id&fold_map$fold==fold_id]
    train_index<-setdiff(seq_len(ncol(logcounts)),test_index)
    preprocessing_seed<-master_seed+repeat_id*1000+fold_id*100
    cat("\nPreparing repeat",repeat_id,"of",n_repeats,"- fold",fold_id,"of",n_folds,"...\n")
    fold_data<-prepare_fold_data(logcounts,train_index,test_index,min_detected_cells,n_hvg,n_pcs,preprocessing_seed)
    y_train<-cell_labels[train_index]
    y_test<-cell_labels[test_index]
    for(representation_id in seq_along(representations_to_run)){
      representation_name<-representations_to_run[representation_id]
      x_train<-fold_data[[representation_name]]$train
      x_test<-fold_data[[representation_name]]$test
      for(model_id in seq_along(models_to_run)){
        model_name<-models_to_run[model_id]
        checkpoint_name<-paste0("repeat",sprintf("%02d",repeat_id),"_fold",sprintf("%02d",fold_id),"_",gsub("[^A-Za-z0-9]","",representation_name),"_",gsub("[^A-Za-z0-9]","",model_name),".rds")
        checkpoint_path<-file.path(checkpoint_dir,checkpoint_name)
        checkpoint_is_current<-FALSE
        if(resume_from_checkpoints&&file.exists(checkpoint_path)){
          if(model_name!="SVM"){
            checkpoint_is_current<-TRUE
          }else{
            existing_checkpoint<-readRDS(checkpoint_path)
            checkpoint_is_current<-"prediction_rule"%in%colnames(existing_checkpoint$metadata)&&existing_checkpoint$metadata$prediction_rule[1]=="native_scaled_class"
            rm(existing_checkpoint)
          }
        }
        if(checkpoint_is_current){
          completed_tasks<-completed_tasks+1
          cat("  [",completed_tasks,"/",expected_tasks,"] Existing checkpoint:",model_name,"-",representation_name,"\n")
          next
        }
        model_seed<-master_seed+repeat_id*10000+fold_id*100+representation_id*10+model_id
        cat("  [",completed_tasks+1,"/",expected_tasks,"] Fitting",model_name,"on",representation_name,"...\n")
        start_time<-proc.time()[3]
        fitted<-fit_predict_model(model_name,x_train,y_train,x_test,class_levels,representation_name,model_seed,n_threads,settings)
        elapsed_seconds<-proc.time()[3]-start_time
        prediction_data<-data.frame(cell_index=test_index,cell_id=cell_ids[test_index],repeat_id=repeat_id,fold=fold_id,model=model_name,representation=representation_name,truth=as.character(y_test),predicted=as.character(fitted$predicted),stringsAsFactors=FALSE)
        prediction_data[probability_columns]<-fitted$probabilities
        importance_data<-fitted$importance
        if(nrow(importance_data)>0){
          importance_data$repeat_id<-repeat_id
          importance_data$fold<-fold_id
        }
        prediction_rule<-if(model_name=="SVM")"native_scaled_class" else "maximum_probability"
        checkpoint<-list(predictions=prediction_data,importance=importance_data,metadata=data.frame(repeat_id=repeat_id,fold=fold_id,model=model_name,representation=representation_name,train_n=length(train_index),test_n=length(test_index),eligible_genes=fold_data$metadata$eligible_genes,selected_genes=fold_data$metadata$selected_genes,pca_components=fold_data$metadata$pca_components,seed=model_seed,elapsed_seconds=elapsed_seconds,prediction_rule=prediction_rule))
        saveRDS(checkpoint,checkpoint_path)
        completed_tasks<-completed_tasks+1
        rm(fitted,checkpoint,prediction_data,importance_data)
        gc(verbose=FALSE)
      }
    }
    rm(fold_data,x_train,x_test)
    gc(verbose=FALSE)
  }
}

# ==============================================================================
# 4. AGGREGATION AND PERFORMANCE ESTIMATION
# ==============================================================================
checkpoint_files<-list.files(checkpoint_dir,pattern="\\.rds$",full.names=TRUE)
if(length(checkpoint_files)!=expected_tasks)stop(paste("Expected",expected_tasks,"checkpoints but found",length(checkpoint_files),". Re-run the script to complete the analysis."))
all_results<-lapply(checkpoint_files,readRDS)
predictions<-dplyr::bind_rows(lapply(all_results,"[[","predictions"))
importance_results<-dplyr::bind_rows(lapply(all_results,"[[","importance"))
run_metadata<-dplyr::bind_rows(lapply(all_results,"[[","metadata"))
write.csv(predictions,file.path(table_dir,"out_of_fold_predictions.csv"),row.names=FALSE)
write.csv(run_metadata,file.path(table_dir,"run_metadata.csv"),row.names=FALSE)

metric_rows<-list()
class_metric_rows<-list()
confusion_rows<-list()
group_keys<-predictions%>%dplyr::distinct(repeat_id,model,representation)%>%dplyr::arrange(repeat_id,model,representation)
for(group_id in seq_len(nrow(group_keys))){
  key<-group_keys[group_id,]
  subset_predictions<-predictions%>%dplyr::filter(repeat_id==key$repeat_id,model==key$model,representation==key$representation)%>%dplyr::arrange(cell_index)
  probability_matrix<-as.matrix(subset_predictions[,probability_columns,drop=FALSE])
  colnames(probability_matrix)<-class_levels
  metrics<-calculate_metrics(subset_predictions$truth,subset_predictions$predicted,probability_matrix,class_levels)
  metric_rows[[group_id]]<-cbind(key,metrics$overall)
  class_metric_rows[[group_id]]<-cbind(key,metrics$per_class)
  confusion_data<-as.data.frame(metrics$confusion,stringsAsFactors=FALSE)
  colnames(confusion_data)<-c("truth","predicted","count")
  confusion_rows[[group_id]]<-cbind(key,confusion_data)
}

metrics_by_repeat<-dplyr::bind_rows(metric_rows)
class_metrics_by_repeat<-dplyr::bind_rows(class_metric_rows)
confusion_by_repeat<-dplyr::bind_rows(confusion_rows)
metric_long<-metrics_by_repeat%>%tidyr::pivot_longer(cols=c(accuracy,balanced_accuracy,macro_precision,macro_recall,macro_f1,macro_auc,log_loss,brier_score),names_to="metric",values_to="value")
class_metric_long<-class_metrics_by_repeat%>%tidyr::pivot_longer(cols=c(precision,recall,f1,auc),names_to="metric",values_to="value")
metric_summary<-summarize_with_ci(metric_long,c("model","representation","metric"))
class_metric_summary<-summarize_with_ci(class_metric_long,c("model","representation","class","metric"))
model_comparisons<-paired_model_comparisons(metric_long)
representation_comparisons<-paired_representation_comparisons(metric_long)

write.csv(metrics_by_repeat,file.path(table_dir,"metrics_by_repeat.csv"),row.names=FALSE)
write.csv(class_metrics_by_repeat,file.path(table_dir,"class_metrics_by_repeat.csv"),row.names=FALSE)
write.csv(confusion_by_repeat,file.path(table_dir,"confusion_matrices_by_repeat.csv"),row.names=FALSE)
write.csv(metric_summary,file.path(table_dir,"metric_summary_95CI.csv"),row.names=FALSE)
write.csv(class_metric_summary,file.path(table_dir,"class_metric_summary_95CI.csv"),row.names=FALSE)
write.csv(model_comparisons,file.path(table_dir,"paired_model_comparisons.csv"),row.names=FALSE)
write.csv(representation_comparisons,file.path(table_dir,"paired_representation_comparisons.csv"),row.names=FALSE)

if(nrow(importance_results)>0){
  total_resamples<-n_repeats*n_folds
  importance_stability<-importance_results%>%dplyr::filter(representation=="Genes")%>%dplyr::group_by(model,feature)%>%dplyr::summarise(selection_count=dplyr::n(),selection_frequency=selection_count/total_resamples,mean_scaled_importance=mean(scaled_importance),median_rank=median(rank),.groups="drop")%>%dplyr::arrange(model,dplyr::desc(selection_frequency),dplyr::desc(mean_scaled_importance))
  write.csv(importance_results,file.path(table_dir,"feature_importance_by_fold.csv"),row.names=FALSE)
  write.csv(importance_stability,file.path(table_dir,"feature_importance_stability.csv"),row.names=FALSE)
}

# ==============================================================================
# 5. PUBLICATION-READY FIGURES
# ==============================================================================
model_order<-c("XGBoost","Random Forest","SVM","Naive Bayes","KNN")
metric_labels<-c(accuracy="Accuracy",balanced_accuracy="Balanced accuracy",macro_f1="Macro F1")
plot_metrics<-metric_long%>%dplyr::filter(metric%in%names(metric_labels))%>%dplyr::mutate(model=factor(model,levels=model_order),metric=factor(metric,levels=names(metric_labels),labels=metric_labels))
p_metrics<-ggplot2::ggplot(plot_metrics,ggplot2::aes(x=model,y=value,color=representation,group=interaction(model,representation)))+ggplot2::geom_boxplot(width=0.55,outlier.shape=NA,position=ggplot2::position_dodge(width=0.65))+ggplot2::geom_point(size=1.8,alpha=0.8,position=ggplot2::position_jitterdodge(jitter.width=0.08,dodge.width=0.65))+ggplot2::facet_wrap(~metric,ncol=1,scales="free_y")+ggplot2::scale_color_manual(values=c("Genes"="#1F77B4","100 PCs"="#D95F02"))+ggplot2::labs(x=NULL,y="Repeated-CV estimate",color="Representation")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=25,hjust=1),legend.position="top",strip.text=ggplot2::element_text(face="bold"))
save_plot(p_metrics,"repeated_cv_metrics",7.2,7.0,figure_dir)

recall_plot_data<-class_metric_summary%>%dplyr::filter(metric=="recall",representation=="Genes")%>%dplyr::mutate(model=factor(model,levels=model_order),class=factor(class,levels=rev(class_levels)))
p_recall<-ggplot2::ggplot(recall_plot_data,ggplot2::aes(x=model,y=class,fill=mean))+ggplot2::geom_tile(color="white",linewidth=0.4)+ggplot2::geom_text(ggplot2::aes(label=sprintf("%.2f",mean)),size=3)+ggplot2::scale_fill_gradient(low="#FFF5EB",high="#7F0000",limits=c(0,1),name="Recall")+ggplot2::labs(x=NULL,y=NULL)+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=25,hjust=1),legend.position="right")
save_plot(p_recall,"class_specific_recall_genes",7.5,4.4,figure_dir)

representation_plot_data<-representation_comparisons%>%dplyr::filter(metric%in%c("accuracy","macro_f1"))%>%dplyr::mutate(model=factor(model,levels=model_order),metric=factor(metric,levels=c("accuracy","macro_f1"),labels=c("Accuracy","Macro F1")))
p_representation<-ggplot2::ggplot(representation_plot_data,ggplot2::aes(x=model,y=mean_difference_pca_minus_genes))+ggplot2::geom_hline(yintercept=0,color="#777777",linetype="dashed")+ggplot2::geom_errorbar(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),width=0.15,color="#34495E")+ggplot2::geom_point(size=2.5,color="#D95F02")+ggplot2::facet_wrap(~metric,ncol=1,scales="free_y")+ggplot2::labs(x=NULL,y="Mean paired difference (100 PCs - genes)")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=25,hjust=1),strip.text=ggplot2::element_text(face="bold"))
save_plot(p_representation,"representation_paired_differences",7.0,5.5,figure_dir)

if(exists("importance_stability")&&nrow(importance_stability)>0){
  stable_features<-importance_stability%>%dplyr::group_by(model)%>%dplyr::slice_max(order_by=selection_frequency,n=15,with_ties=FALSE)%>%dplyr::ungroup()%>%dplyr::mutate(feature_model=paste(feature,model,sep="___"))
  p_importance<-ggplot2::ggplot(stable_features,ggplot2::aes(x=selection_frequency,y=reorder(feature_model,selection_frequency),fill=mean_scaled_importance))+ggplot2::geom_col()+ggplot2::facet_wrap(~model,scales="free_y",ncol=2)+ggplot2::scale_y_discrete(labels=function(x)sub("___.*$","",x))+ggplot2::scale_fill_viridis_c(option="C",name="Mean scaled\nimportance")+ggplot2::labs(x="Proportion of folds in top features",y=NULL)+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),strip.text=ggplot2::element_text(face="bold"))
  save_plot(p_importance,"feature_importance_stability",8.0,6.2,figure_dir)
}

# ==============================================================================
# 6. REPRODUCIBILITY RECORD
# ==============================================================================
method_specifications<-data.frame(model=c("Random Forest","SVM","XGBoost","Naive Bayes","KNN"),specification=c(paste0("ntree=",rf_ntree,"; mtry=floor(sqrt(p))"),paste0("radial kernel; cost=",svm_cost,"; gamma=1/p; training-fold scaling enabled"),paste0("nrounds=",xgb_nrounds,"; max_depth=",xgb_max_depth,"; eta=",xgb_eta,"; subsample=",xgb_subsample,"; colsample_bytree=",xgb_colsample_bytree),"Gaussian class-conditional densities; Laplace=0",paste0("k=",knn_k,"; Euclidean distance; uniform voting")))
write.csv(method_specifications,file.path(table_dir,"method_specifications.csv"),row.names=FALSE)
package_versions<-data.frame(package=required_packages,version=vapply(required_packages,function(x)as.character(utils::packageVersion(x)),character(1)))
write.csv(package_versions,file.path(table_dir,"package_versions.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(results_dir,"sessionInfo.txt"))

cat("\n==============================================================================\n")
cat("ANALYSIS COMPLETED\n")
cat("==============================================================================\n")
cat("Primary summary:\n")
print(metric_summary%>%dplyr::filter(metric%in%c("accuracy","balanced_accuracy","macro_f1"))%>%dplyr::arrange(metric,dplyr::desc(mean)))
cat("\nResults saved to:",results_dir,"\n")
cat("Upload the complete results folder for manuscript revision.\n")
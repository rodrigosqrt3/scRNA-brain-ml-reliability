# ==============================================================================
# ERAMIA-RS 2026 - C1-RUN-BLOCKED REPEATED CROSS-VALIDATION
# ==============================================================================
# Cell identifiers contain a run prefix before the underscore. The 3,005 cells
# form 76 such groups. Entire runs are assigned to folds, preventing cells from
# the same C1 run from appearing in training and validation simultaneously.
# All learned preprocessing operations remain restricted to training folds.
# ==============================================================================

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
master_seed<-335705
n_folds<-5
n_repeats<-5
group_optimization_restarts<-250
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
save_top_features<-30

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
analysis_id<-paste0("grouped",n_repeats,"x",n_folds,"_c1run_hvg",n_hvg,"_pc",n_pcs)
results_dir<-file.path(project_dir,"results",analysis_id)
checkpoint_dir<-file.path(results_dir,"checkpoints")
table_dir<-file.path(results_dir,"tables")
figure_dir<-file.path(results_dir,"figures")
dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)

required_packages<-c("scRNAseq","SummarizedExperiment","Matrix","matrixStats","randomForest","e1071","xgboost","caret","pROC","irlba","ggplot2","dplyr","tidyr","scales")
missing_packages<-required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)>0)stop(paste("Missing packages:",paste(missing_packages,collapse=", ")))
`%>%`<-dplyr::`%>%`
options(stringsAsFactors=FALSE)
set.seed(master_seed)

settings<-list(rf_ntree=rf_ntree,svm_cost=svm_cost,xgb_nrounds=xgb_nrounds,xgb_max_depth=xgb_max_depth,xgb_eta=xgb_eta,xgb_subsample=xgb_subsample,xgb_colsample_bytree=xgb_colsample_bytree,knn_k=knn_k,save_top_features=save_top_features)
configuration<-data.frame(parameter=c("master_seed","n_folds","n_repeats","group_optimization_restarts","min_detected_cells","n_hvg","n_pcs","library_scale","rf_ntree","svm_cost","xgb_nrounds","xgb_max_depth","xgb_eta","xgb_subsample","xgb_colsample_bytree","knn_k","n_threads"),value=as.character(c(master_seed,n_folds,n_repeats,group_optimization_restarts,min_detected_cells,n_hvg,n_pcs,library_scale,rf_ntree,svm_cost,xgb_nrounds,xgb_max_depth,xgb_eta,xgb_subsample,xgb_colsample_bytree,knn_k,n_threads)))
configuration_path<-file.path(table_dir,"analysis_configuration.csv")
if(file.exists(configuration_path)&&length(list.files(checkpoint_dir,pattern="\\.rds$"))>0){
  previous_configuration<-read.csv(configuration_path,stringsAsFactors=FALSE)
  if(!identical(paste(previous_configuration$parameter,previous_configuration$value),paste(configuration$parameter,configuration$value)))stop("Configuration changed while checkpoints exist. Use a new results directory or restore the previous settings.")
}else{
  write.csv(configuration,configuration_path,row.names=FALSE)
}

cat("\n==============================================================================\n")
cat("C1-RUN-BLOCKED REPEATED CROSS-VALIDATION\n")
cat("==============================================================================\n")
cat("Analysis ID:",analysis_id,"\n")
cat("Resampling:",n_repeats,"repeats x",n_folds,"group-blocked folds\n")
cat("Models:",paste(models_to_run,collapse=", "),"\n\n")

# ==============================================================================
# 2. DATA IMPORT AND C1 RUN IDENTIFICATION
# ==============================================================================
sce<-scRNAseq::ZeiselBrainData()
counts_raw<-SummarizedExperiment::assay(sce,"counts")
if(!inherits(counts_raw,"sparseMatrix"))counts_raw<-Matrix::Matrix(counts_raw,sparse=TRUE)
cell_labels<-factor(sce$level1class)
valid_cells<-!is.na(cell_labels)&Matrix::colSums(counts_raw)>0
counts_raw<-counts_raw[,valid_cells,drop=FALSE]
cell_labels<-droplevels(cell_labels[valid_cells])
cell_ids<-colnames(counts_raw)
if(is.null(cell_ids))stop("Cell identifiers are required to derive the C1 run groups.")
c1_run<-sub("_.*$","",cell_ids)
if(any(c1_run==cell_ids))stop("At least one cell identifier does not contain the expected run-well separator.")
if(length(unique(c1_run))<n_folds)stop("Too few C1 run groups for grouped cross-validation.")

library_sizes<-Matrix::colSums(counts_raw)
logcounts<-counts_raw%*%Matrix::Diagonal(x=library_scale/library_sizes)
logcounts<-log1p(logcounts)

run_distribution<-data.frame(c1_run=c1_run,class=as.character(cell_labels))%>%dplyr::count(c1_run,class,name="n")%>%dplyr::group_by(c1_run)%>%dplyr::mutate(run_total=sum(n),class_proportion=n/run_total)%>%dplyr::ungroup()
run_summary<-run_distribution%>%dplyr::group_by(c1_run)%>%dplyr::summarise(n_cells=sum(n),n_classes=sum(n>0),.groups="drop")
class_run_coverage<-run_distribution%>%dplyr::filter(n>0)%>%dplyr::group_by(class)%>%dplyr::summarise(n_runs=dplyr::n_distinct(c1_run),n_cells=sum(n),.groups="drop")
write.csv(run_distribution,file.path(table_dir,"c1_run_class_distribution.csv"),row.names=FALSE)
write.csv(run_summary,file.path(table_dir,"c1_run_summary.csv"),row.names=FALSE)
write.csv(class_run_coverage,file.path(table_dir,"class_run_coverage.csv"),row.names=FALSE)
cat("Cells:",ncol(counts_raw),"| C1 runs:",length(unique(c1_run)),"| Classes:",nlevels(cell_labels),"\n")
print(class_run_coverage)

# ==============================================================================
# 3. REPEATED GROUP-STRATIFIED FOLDS
# ==============================================================================
fold_objects<-create_repeated_group_stratified_folds(cell_labels,c1_run,n_folds,n_repeats,master_seed,group_optimization_restarts)
fold_map<-fold_objects$cell_folds
group_assignments<-fold_objects$group_assignments
fold_map$cell_id<-cell_ids[fold_map$cell_index]
fold_map$class<-as.character(cell_labels[fold_map$cell_index])
write.csv(fold_map,file.path(table_dir,"grouped_resampling_folds.csv"),row.names=FALSE)
write.csv(group_assignments,file.path(table_dir,"group_assignments.csv"),row.names=FALSE)

fold_balance<-fold_map%>%dplyr::count(repeat_id,fold,class,name="n")%>%dplyr::group_by(repeat_id,fold)%>%dplyr::mutate(fold_total=sum(n),class_proportion=n/fold_total)%>%dplyr::ungroup()
fold_group_counts<-fold_map%>%dplyr::distinct(repeat_id,fold,group_id)%>%dplyr::count(repeat_id,fold,name="n_groups")
write.csv(fold_balance,file.path(table_dir,"fold_class_balance.csv"),row.names=FALSE)
write.csv(fold_group_counts,file.path(table_dir,"fold_group_counts.csv"),row.names=FALSE)

leakage_check<-fold_map%>%dplyr::distinct(repeat_id,group_id,fold)%>%dplyr::count(repeat_id,group_id,name="n_folds")
if(any(leakage_check$n_folds!=1))stop("Group leakage detected in the blocked folds.")

class_levels<-levels(cell_labels)
probability_columns<-paste0("prob__",make.names(class_levels,unique=TRUE))
expected_tasks<-n_repeats*n_folds*length(representations_to_run)*length(models_to_run)
completed_tasks<-0

# ==============================================================================
# 4. MODEL FITTING WITH CHECKPOINTS
# ==============================================================================
for(repeat_id in seq_len(n_repeats)){
  for(fold_id in seq_len(n_folds)){
    test_index<-fold_map$cell_index[fold_map$repeat_id==repeat_id&fold_map$fold==fold_id]
    train_index<-setdiff(seq_len(ncol(logcounts)),test_index)
    training_groups<-unique(c1_run[train_index])
    validation_groups<-unique(c1_run[test_index])
    if(length(intersect(training_groups,validation_groups))>0)stop("C1 run leakage detected.")
    preprocessing_seed<-master_seed+repeat_id*1000+fold_id*100
    cat("\nPreparing grouped repeat",repeat_id,"of",n_repeats,"- fold",fold_id,"of",n_folds,"...\n")
    cat("Training runs:",length(training_groups),"| Validation runs:",length(validation_groups),"\n")
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
        if(resume_from_checkpoints&&file.exists(checkpoint_path)){
          completed_tasks<-completed_tasks+1
          cat("  [",completed_tasks,"/",expected_tasks,"] Existing checkpoint:",model_name,"-",representation_name,"\n")
          next
        }
        model_seed<-master_seed+repeat_id*10000+fold_id*100+representation_id*10+model_id
        cat("  [",completed_tasks+1,"/",expected_tasks,"] Fitting",model_name,"on",representation_name,"...\n")
        start_time<-proc.time()[3]
        fitted<-fit_predict_model(model_name,x_train,y_train,x_test,class_levels,representation_name,model_seed,n_threads,settings)
        elapsed_seconds<-proc.time()[3]-start_time
        prediction_data<-data.frame(cell_index=test_index,cell_id=cell_ids[test_index],c1_run=c1_run[test_index],repeat_id=repeat_id,fold=fold_id,model=model_name,representation=representation_name,truth=as.character(y_test),predicted=as.character(fitted$predicted),stringsAsFactors=FALSE)
        prediction_data[probability_columns]<-fitted$probabilities
        importance_data<-fitted$importance
        if(nrow(importance_data)>0){
          importance_data$repeat_id<-repeat_id
          importance_data$fold<-fold_id
        }
        checkpoint<-list(predictions=prediction_data,importance=importance_data,metadata=data.frame(repeat_id=repeat_id,fold=fold_id,model=model_name,representation=representation_name,train_n=length(train_index),test_n=length(test_index),training_runs=length(training_groups),validation_runs=length(validation_groups),eligible_genes=fold_data$metadata$eligible_genes,selected_genes=fold_data$metadata$selected_genes,pca_components=fold_data$metadata$pca_components,seed=model_seed,elapsed_seconds=elapsed_seconds,grouping_rule="cell_id prefix before underscore"))
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
# 5. AGGREGATION AND COMPARISON WITH CELL-STRATIFIED CV
# ==============================================================================
checkpoint_files<-list.files(checkpoint_dir,pattern="\\.rds$",full.names=TRUE)
if(length(checkpoint_files)!=expected_tasks)stop(paste("Expected",expected_tasks,"checkpoints but found",length(checkpoint_files),"."))
all_results<-lapply(checkpoint_files,readRDS)
predictions<-dplyr::bind_rows(lapply(all_results,"[[","predictions"))
importance_results<-dplyr::bind_rows(lapply(all_results,"[[","importance"))
run_metadata<-dplyr::bind_rows(lapply(all_results,"[[","metadata"))
write.csv(predictions,file.path(table_dir,"out_of_fold_predictions_grouped.csv"),row.names=FALSE)
write.csv(run_metadata,file.path(table_dir,"run_metadata.csv"),row.names=FALSE)

metric_rows<-list()
class_metric_rows<-list()
group_keys<-predictions%>%dplyr::distinct(repeat_id,model,representation)%>%dplyr::arrange(repeat_id,model,representation)
for(group_id in seq_len(nrow(group_keys))){
  key<-group_keys[group_id,]
  subset_predictions<-predictions%>%dplyr::filter(repeat_id==key$repeat_id,model==key$model,representation==key$representation)%>%dplyr::arrange(cell_index)
  probability_matrix<-as.matrix(subset_predictions[,probability_columns,drop=FALSE])
  colnames(probability_matrix)<-class_levels
  metrics<-calculate_metrics(subset_predictions$truth,subset_predictions$predicted,probability_matrix,class_levels)
  metric_rows[[group_id]]<-cbind(key,metrics$overall)
  class_metric_rows[[group_id]]<-cbind(key,metrics$per_class)
}

metrics_by_repeat<-dplyr::bind_rows(metric_rows)
class_metrics_by_repeat<-dplyr::bind_rows(class_metric_rows)
metric_long<-metrics_by_repeat%>%tidyr::pivot_longer(cols=c(accuracy,balanced_accuracy,macro_precision,macro_recall,macro_f1,macro_auc,log_loss,brier_score),names_to="metric",values_to="value")
class_metric_long<-class_metrics_by_repeat%>%tidyr::pivot_longer(cols=c(precision,recall,f1,auc),names_to="metric",values_to="value")
metric_summary<-summarize_with_ci(metric_long,c("model","representation","metric"))
class_metric_summary<-summarize_with_ci(class_metric_long,c("model","representation","class","metric"))
write.csv(metrics_by_repeat,file.path(table_dir,"metrics_by_repeat.csv"),row.names=FALSE)
write.csv(class_metrics_by_repeat,file.path(table_dir,"class_metrics_by_repeat.csv"),row.names=FALSE)
write.csv(metric_summary,file.path(table_dir,"metric_summary_95CI.csv"),row.names=FALSE)
write.csv(class_metric_summary,file.path(table_dir,"class_metric_summary_95CI.csv"),row.names=FALSE)

primary_summary_paths<-c(file.path(project_dir,"results","repeated5x5_hvg2000_pc100","tables","metric_summary_95CI.csv"),file.path(project_dir,"tables","metric_summary_95CI.csv"))
primary_summary_path<-primary_summary_paths[file.exists(primary_summary_paths)][1]
if(length(primary_summary_path)>0&&!is.na(primary_summary_path)){
  primary_summary<-read.csv(primary_summary_path,stringsAsFactors=FALSE)%>%dplyr::select(model,representation,metric,cell_stratified_mean=mean,cell_stratified_sd=sd)
  validation_comparison<-metric_summary%>%dplyr::select(model,representation,metric,grouped_mean=mean,grouped_sd=sd)%>%dplyr::left_join(primary_summary,by=c("model","representation","metric"))%>%dplyr::mutate(grouped_minus_cell_stratified=grouped_mean-cell_stratified_mean)
  write.csv(validation_comparison,file.path(table_dir,"validation_scheme_comparison.csv"),row.names=FALSE)
}

if(nrow(importance_results)>0){
  total_resamples<-n_repeats*n_folds
  importance_stability<-importance_results%>%dplyr::filter(representation=="Genes")%>%dplyr::group_by(model,feature)%>%dplyr::summarise(selection_count=dplyr::n(),selection_frequency=selection_count/total_resamples,mean_scaled_importance=mean(scaled_importance),median_rank=median(rank),.groups="drop")%>%dplyr::arrange(model,dplyr::desc(selection_frequency),dplyr::desc(mean_scaled_importance))
  write.csv(importance_results,file.path(table_dir,"feature_importance_by_fold.csv"),row.names=FALSE)
  write.csv(importance_stability,file.path(table_dir,"feature_importance_stability.csv"),row.names=FALSE)
}

# ==============================================================================
# 6. FIGURES AND REPRODUCIBILITY
# ==============================================================================
model_order<-c("XGBoost","Random Forest","SVM","Naive Bayes","KNN")
plot_data<-metric_summary%>%dplyr::filter(metric%in%c("accuracy","balanced_accuracy","macro_f1"))%>%dplyr::mutate(model=factor(model,levels=model_order),metric=factor(metric,levels=c("accuracy","balanced_accuracy","macro_f1"),labels=c("Accuracy","Balanced accuracy","Macro F1")))
p_grouped<-ggplot2::ggplot(plot_data,ggplot2::aes(x=model,y=mean,color=representation,group=representation))+ggplot2::geom_point(position=ggplot2::position_dodge(width=0.45),size=2.3)+ggplot2::geom_errorbar(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),position=ggplot2::position_dodge(width=0.45),width=0.18)+ggplot2::facet_wrap(~metric,ncol=1,scales="free_y")+ggplot2::scale_color_manual(values=c("Genes"="#1F77B4","100 PCs"="#D95F02"))+ggplot2::labs(x=NULL,y="C1-run-blocked estimate",color="Representation")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=25,hjust=1),legend.position="top",strip.text=ggplot2::element_text(face="bold"))
save_plot(p_grouped,"grouped_cv_metrics",7.2,7.0,figure_dir)

if(exists("validation_comparison")){
  comparison_plot_data<-validation_comparison%>%dplyr::filter(metric%in%c("accuracy","balanced_accuracy","macro_f1"))%>%dplyr::mutate(model=factor(model,levels=model_order),metric=factor(metric,levels=c("accuracy","balanced_accuracy","macro_f1"),labels=c("Accuracy","Balanced accuracy","Macro F1")))
  p_gap<-ggplot2::ggplot(comparison_plot_data,ggplot2::aes(x=model,y=grouped_minus_cell_stratified,fill=representation))+ggplot2::geom_hline(yintercept=0,linetype="dashed",color="#666666")+ggplot2::geom_col(position=ggplot2::position_dodge(width=0.75),width=0.65)+ggplot2::facet_wrap(~metric,ncol=1,scales="free_y")+ggplot2::scale_fill_manual(values=c("Genes"="#1F77B4","100 PCs"="#D95F02"))+ggplot2::labs(x=NULL,y="Grouped minus cell-stratified estimate",fill="Representation")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=25,hjust=1),legend.position="top",strip.text=ggplot2::element_text(face="bold"))
  save_plot(p_gap,"validation_generalization_gap",7.2,7.0,figure_dir)
}

package_versions<-data.frame(package=required_packages,version=vapply(required_packages,function(x)as.character(utils::packageVersion(x)),character(1)))
write.csv(package_versions,file.path(table_dir,"package_versions.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(results_dir,"sessionInfo.txt"))

cat("\n==============================================================================\n")
cat("GROUPED ANALYSIS COMPLETED\n")
cat("==============================================================================\n")
print(metric_summary%>%dplyr::filter(metric%in%c("accuracy","balanced_accuracy","macro_f1"))%>%dplyr::arrange(metric,dplyr::desc(mean)))
cat("\nResults saved to:",results_dir,"\n")

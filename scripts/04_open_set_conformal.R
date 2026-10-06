# ==============================================================================
# ERAMIA-RS 2026 - GROUP-BLOCKED OPEN-SET AND CONFORMAL ANALYSIS
# ==============================================================================
# Each level-one class is omitted from model training in turn and treated as an
# unknown population. C1 runs are assigned to training, calibration and test
# partitions without overlap. SVM and XGBoost are evaluated with gene features.
# Confidence rejection and class-conditional split-conformal prediction sets
# are calculated from a separate calibration fold.
# ==============================================================================

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
master_seed<-335705
n_folds<-5
n_repeats<-3
group_optimization_restarts<-250
min_detected_cells<-10
n_hvg<-2000
n_pcs<-100
library_scale<-10000
known_acceptance_target<-0.95
conformal_alpha<-c(0.05,0.10,0.20)
models_to_run<-c("SVM","XGBoost")
resume_from_checkpoints<-TRUE
detected_cores<-parallel::detectCores(logical=FALSE)
if(is.na(detected_cores))detected_cores<-1
n_threads<-max(1,detected_cores-1)

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
analysis_id<-paste0("open_set_conformal_grouped",n_repeats,"x",n_folds,"_hvg",n_hvg)
results_dir<-file.path(project_dir,"results",analysis_id)
checkpoint_dir<-file.path(results_dir,"checkpoints")
table_dir<-file.path(results_dir,"tables")
figure_dir<-file.path(results_dir,"figures")
dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)

required_packages<-c("scRNAseq","SummarizedExperiment","Matrix","matrixStats","e1071","xgboost","pROC","irlba","ggplot2","dplyr","tidyr","scales")
missing_packages<-required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)>0)stop(paste("Missing packages:",paste(missing_packages,collapse=", ")))
`%>%`<-dplyr::`%>%`
options(stringsAsFactors=FALSE)
set.seed(master_seed)

settings<-list(rf_ntree=rf_ntree,svm_cost=svm_cost,xgb_nrounds=xgb_nrounds,xgb_max_depth=xgb_max_depth,xgb_eta=xgb_eta,xgb_subsample=xgb_subsample,xgb_colsample_bytree=xgb_colsample_bytree,knn_k=knn_k,save_top_features=20)
configuration<-data.frame(parameter=c("master_seed","n_folds","n_repeats","group_optimization_restarts","min_detected_cells","n_hvg","library_scale","known_acceptance_target","conformal_alpha","models","n_threads"),value=c(master_seed,n_folds,n_repeats,group_optimization_restarts,min_detected_cells,n_hvg,library_scale,known_acceptance_target,paste(conformal_alpha,collapse=";"),paste(models_to_run,collapse=";"),n_threads))
configuration_path<-file.path(table_dir,"analysis_configuration.csv")
if(file.exists(configuration_path)&&length(list.files(checkpoint_dir,pattern="\\.rds$"))>0){
  previous_configuration<-read.csv(configuration_path,stringsAsFactors=FALSE)
  if(!identical(paste(previous_configuration$parameter,previous_configuration$value),paste(configuration$parameter,configuration$value)))stop("Configuration changed while checkpoints exist. Restore the previous settings or use a new results directory.")
}else{
  write.csv(configuration,configuration_path,row.names=FALSE)
}

average_precision<-function(labels,scores){
  labels<-as.integer(labels)
  ordering<-order(scores,decreasing=TRUE)
  labels<-labels[ordering]
  positives<-sum(labels==1)
  if(positives==0)return(NA_real_)
  precision_at_rank<-cumsum(labels==1)/seq_along(labels)
  sum(precision_at_rank[labels==1])/positives
}

calculate_open_metrics<-function(data,class_levels){
  known<-data$truth_open!="Unknown"
  unknown<-!known
  open_levels<-c(class_levels,"Unknown")
  confusion<-table(factor(data$truth_open,levels=open_levels),factor(data$open_prediction,levels=open_levels))
  class_recall<-safe_divide(diag(confusion),rowSums(confusion))
  response<-as.integer(unknown)
  unknown_auc<-tryCatch(as.numeric(pROC::auc(pROC::roc(response,data$unknown_score,levels=c(0,1),direction="<",quiet=TRUE))),error=function(e)NA_real_)
  data.frame(known_classification_accuracy=mean(data$predicted[known]==data$truth[known]),known_acceptance_rate=mean(data$accepted[known]),unknown_rejection_rate=mean(!data$accepted[unknown]),accepted_known_accuracy=if(any(data$accepted[known]))mean(data$predicted[known][data$accepted[known]]==data$truth[known][data$accepted[known]]) else NA_real_,open_set_accuracy=mean(data$open_prediction==data$truth_open),open_set_balanced_accuracy=mean(class_recall,na.rm=TRUE),unknown_detection_auc=unknown_auc,unknown_average_precision=average_precision(response,data$unknown_score),mean_known_confidence=mean(data$max_probability[known]),mean_unknown_confidence=mean(data$max_probability[unknown]))
}

conformal_pvalues<-function(calibration_probability,calibration_truth,test_probability,class_levels){
  output<-matrix(NA_real_,nrow=nrow(test_probability),ncol=length(class_levels),dimnames=list(NULL,class_levels))
  for(class_id in seq_along(class_levels)){
    class_name<-class_levels[class_id]
    calibration_indices<-which(calibration_truth==class_name)
    if(length(calibration_indices)==0)stop(paste("Calibration fold has no observations for class",class_name))
    calibration_scores<-1-calibration_probability[calibration_indices,class_id]
    test_scores<-1-test_probability[,class_id]
    output[,class_id]<-vapply(test_scores,function(score)(1+sum(calibration_scores>=score))/(length(calibration_scores)+1),numeric(1))
  }
  output
}

cat("\n==============================================================================\n")
cat("GROUP-BLOCKED OPEN-SET AND CONFORMAL ANALYSIS\n")
cat("==============================================================================\n")
cat("Models:",paste(models_to_run,collapse=", "),"\n")
cat("Unknown-class experiments: all level-one classes\n")
cat("Grouped repeats:",n_repeats,"| Outer folds:",n_folds,"\n\n")

# ==============================================================================
# 2. DATA AND GROUPED PARTITIONS
# ==============================================================================
sce<-scRNAseq::ZeiselBrainData()
counts_raw<-SummarizedExperiment::assay(sce,"counts")
if(!inherits(counts_raw,"sparseMatrix"))counts_raw<-Matrix::Matrix(counts_raw,sparse=TRUE)
cell_labels<-factor(sce$level1class)
valid_cells<-!is.na(cell_labels)&Matrix::colSums(counts_raw)>0
counts_raw<-counts_raw[,valid_cells,drop=FALSE]
cell_labels<-droplevels(cell_labels[valid_cells])
cell_ids<-colnames(counts_raw)
if(is.null(cell_ids))stop("Cell identifiers are required.")
c1_run<-sub("_.*$","",cell_ids)
if(any(c1_run==cell_ids))stop("Unexpected cell identifier format.")

library_sizes<-Matrix::colSums(counts_raw)
logcounts<-counts_raw%*%Matrix::Diagonal(x=library_scale/library_sizes)
logcounts<-log1p(logcounts)
all_class_levels<-levels(cell_labels)

fold_objects<-create_repeated_group_stratified_folds(cell_labels,c1_run,n_folds,n_repeats,master_seed,group_optimization_restarts)
fold_map<-fold_objects$cell_folds
fold_map$cell_id<-cell_ids[fold_map$cell_index]
fold_map$class<-as.character(cell_labels[fold_map$cell_index])
write.csv(fold_map,file.path(table_dir,"open_set_grouped_folds.csv"),row.names=FALSE)
write.csv(fold_objects$group_assignments,file.path(table_dir,"group_assignments.csv"),row.names=FALSE)

expected_tasks<-n_repeats*n_folds*length(all_class_levels)*length(models_to_run)
completed_tasks<-0

# ==============================================================================
# 3. LEAVE-ONE-CLASS-OUT MODEL FITTING
# ==============================================================================
for(repeat_id in seq_len(n_repeats)){
  for(test_fold in seq_len(n_folds)){
    calibration_fold<-if(test_fold==n_folds)1 else test_fold+1
    training_folds<-setdiff(seq_len(n_folds),c(test_fold,calibration_fold))
    repeat_map<-fold_map%>%dplyr::filter(repeat_id==!!repeat_id)
    outer_test_index<-repeat_map$cell_index[repeat_map$fold==test_fold]
    outer_calibration_index<-repeat_map$cell_index[repeat_map$fold==calibration_fold]
    outer_training_index<-repeat_map$cell_index[repeat_map$fold%in%training_folds]
    if(length(intersect(c1_run[outer_training_index],c1_run[c(outer_calibration_index,outer_test_index)]))>0)stop("C1 run leakage detected.")
    for(omitted_id in seq_along(all_class_levels)){
      omitted_class<-all_class_levels[omitted_id]
      known_levels<-setdiff(all_class_levels,omitted_class)
      train_index<-outer_training_index[cell_labels[outer_training_index]!=omitted_class]
      calibration_index<-outer_calibration_index[cell_labels[outer_calibration_index]!=omitted_class]
      test_index<-outer_test_index
      if(!all(known_levels%in%as.character(cell_labels[train_index])))stop("Training partition is missing a known class.")
      if(!all(known_levels%in%as.character(cell_labels[calibration_index])))stop("Calibration partition is missing a known class.")
      if(!omitted_class%in%as.character(cell_labels[test_index]))stop("Test partition is missing the omitted class.")
      evaluation_index<-c(calibration_index,test_index)
      preprocessing_seed<-master_seed+repeat_id*100000+test_fold*1000+omitted_id*10
      cat("\nRepeat",repeat_id,"/",n_repeats,"| Test fold",test_fold,"/",n_folds,"| Unknown:",omitted_class,"\n")
      fold_data<-prepare_fold_data(logcounts,train_index,evaluation_index,min_detected_cells,n_hvg,n_pcs,preprocessing_seed)
      x_train<-fold_data$Genes$train
      x_evaluation<-fold_data$Genes$test
      calibration_position<-seq_along(calibration_index)
      test_position<-(length(calibration_index)+1):length(evaluation_index)
      y_train<-factor(cell_labels[train_index],levels=known_levels)
      calibration_truth<-as.character(cell_labels[calibration_index])
      test_truth<-as.character(cell_labels[test_index])
      for(model_id in seq_along(models_to_run)){
        model_name<-models_to_run[model_id]
        checkpoint_name<-paste0("repeat",sprintf("%02d",repeat_id),"_testfold",sprintf("%02d",test_fold),"_unknown_",gsub("[^A-Za-z0-9]","",omitted_class),"_",gsub("[^A-Za-z0-9]","",model_name),".rds")
        checkpoint_path<-file.path(checkpoint_dir,checkpoint_name)
        if(resume_from_checkpoints&&file.exists(checkpoint_path)){
          completed_tasks<-completed_tasks+1
          cat("  [",completed_tasks,"/",expected_tasks,"] Existing checkpoint:",model_name,"\n")
          next
        }
        model_seed<-master_seed+repeat_id*1000000+test_fold*10000+omitted_id*100+model_id
        cat("  [",completed_tasks+1,"/",expected_tasks,"] Fitting",model_name,"...\n")
        start_time<-proc.time()[3]
        fitted<-fit_predict_model(model_name,x_train,y_train,x_evaluation,known_levels,"Genes",model_seed,n_threads,settings)
        elapsed_seconds<-proc.time()[3]-start_time
        calibration_probability<-fitted$probabilities[calibration_position,,drop=FALSE]
        test_probability<-fitted$probabilities[test_position,,drop=FALSE]
        test_predicted<-as.character(fitted$predicted[test_position])
        calibration_confidence<-apply(calibration_probability,1,max)
        confidence_threshold<-as.numeric(stats::quantile(calibration_confidence,probs=1-known_acceptance_target,type=1,names=FALSE))
        test_confidence<-apply(test_probability,1,max)
        accepted<-test_confidence>=confidence_threshold
        truth_open<-ifelse(test_truth==omitted_class,"Unknown",test_truth)
        open_prediction<-ifelse(accepted,test_predicted,"Unknown")
        unknown_score<-1-test_confidence
        pvalues<-conformal_pvalues(calibration_probability,calibration_truth,test_probability,known_levels)
        
        prediction_data<-data.frame(cell_index=test_index,cell_id=cell_ids[test_index],c1_run=c1_run[test_index],repeat_id=repeat_id,test_fold=test_fold,calibration_fold=calibration_fold,omitted_class=omitted_class,model=model_name,truth=test_truth,truth_open=truth_open,predicted=test_predicted,max_probability=test_confidence,confidence_threshold=confidence_threshold,accepted=accepted,open_prediction=open_prediction,unknown_score=unknown_score,stringsAsFactors=FALSE)
        probability_columns<-paste0("prob__",make.names(known_levels,unique=TRUE))
        prediction_data[probability_columns]<-test_probability
        pvalue_columns<-paste0("pvalue__",make.names(known_levels,unique=TRUE))
        prediction_data[pvalue_columns]<-pvalues
        
        conformal_rows<-list()
        for(alpha_id in seq_along(conformal_alpha)){
          alpha<-conformal_alpha[alpha_id]
          included<-pvalues>alpha
          set_size<-rowSums(included)
          truth_position<-match(test_truth,known_levels)
          truth_in_set<-rep(NA,length(test_index))
          known_test<-test_truth!=omitted_class
          truth_in_set[known_test]<-included[cbind(which(known_test),truth_position[known_test])]
          set_label<-apply(included,1,function(row){
            selected<-known_levels[row]
            if(length(selected)==0)"[empty]" else paste(selected,collapse=";")
          })
          conformal_rows[[alpha_id]]<-data.frame(cell_index=test_index,repeat_id=repeat_id,test_fold=test_fold,omitted_class=omitted_class,model=model_name,alpha=alpha,truth=test_truth,truth_open=truth_open,set_size=set_size,empty_set=set_size==0,singleton_set=set_size==1,truth_in_set=truth_in_set,prediction_set=set_label,stringsAsFactors=FALSE)
        }
        conformal_data<-dplyr::bind_rows(conformal_rows)
        metadata<-data.frame(repeat_id=repeat_id,test_fold=test_fold,calibration_fold=calibration_fold,omitted_class=omitted_class,model=model_name,train_n=length(train_index),calibration_n=length(calibration_index),test_n=length(test_index),unknown_test_n=sum(test_truth==omitted_class),known_test_n=sum(test_truth!=omitted_class),training_runs=length(unique(c1_run[train_index])),calibration_runs=length(unique(c1_run[calibration_index])),test_runs=length(unique(c1_run[test_index])),confidence_threshold=confidence_threshold,eligible_genes=fold_data$metadata$eligible_genes,selected_genes=fold_data$metadata$selected_genes,seed=model_seed,elapsed_seconds=elapsed_seconds)
        checkpoint<-list(predictions=prediction_data,conformal=conformal_data,metadata=metadata)
        saveRDS(checkpoint,checkpoint_path)
        completed_tasks<-completed_tasks+1
        rm(fitted,checkpoint,prediction_data,conformal_data)
        gc(verbose=FALSE)
      }
      rm(fold_data,x_train,x_evaluation)
      gc(verbose=FALSE)
    }
  }
}

# ==============================================================================
# 4. AGGREGATION ACROSS COMPLETE OUT-OF-FOLD REPEATS
# ==============================================================================
checkpoint_files<-list.files(checkpoint_dir,pattern="\\.rds$",full.names=TRUE)
if(length(checkpoint_files)!=expected_tasks)stop(paste("Expected",expected_tasks,"checkpoints but found",length(checkpoint_files),"."))
all_results<-lapply(checkpoint_files,readRDS)
predictions<-dplyr::bind_rows(lapply(all_results,"[[","predictions"))
conformal_predictions<-dplyr::bind_rows(lapply(all_results,"[[","conformal"))
run_metadata<-dplyr::bind_rows(lapply(all_results,"[[","metadata"))
write.csv(predictions,file.path(table_dir,"open_set_predictions.csv"),row.names=FALSE)
write.csv(conformal_predictions,file.path(table_dir,"conformal_prediction_sets.csv"),row.names=FALSE)
write.csv(run_metadata,file.path(table_dir,"run_metadata.csv"),row.names=FALSE)

metric_rows<-list()
conformal_metric_rows<-list()
forced_assignment_rows<-list()
group_keys<-predictions%>%dplyr::distinct(repeat_id,omitted_class,model)%>%dplyr::arrange(repeat_id,omitted_class,model)
for(group_id in seq_len(nrow(group_keys))){
  key<-group_keys[group_id,]
  subset_data<-predictions%>%dplyr::filter(repeat_id==key$repeat_id,omitted_class==key$omitted_class,model==key$model)
  known_levels<-setdiff(all_class_levels,key$omitted_class)
  metric_rows[[group_id]]<-cbind(key,calculate_open_metrics(subset_data,known_levels))
  forced_assignment_rows[[group_id]]<-subset_data%>%dplyr::filter(truth_open=="Unknown")%>%dplyr::count(predicted,name="n")%>%dplyr::mutate(proportion=n/sum(n),repeat_id=key$repeat_id,omitted_class=key$omitted_class,model=key$model)
  conformal_subset<-conformal_predictions%>%dplyr::filter(repeat_id==key$repeat_id,omitted_class==key$omitted_class,model==key$model)
  conformal_metric_rows[[group_id]]<-conformal_subset%>%dplyr::group_by(alpha)%>%dplyr::summarise(known_coverage=mean(truth_in_set[truth_open!="Unknown"],na.rm=TRUE),known_mean_set_size=mean(set_size[truth_open!="Unknown"]),known_singleton_rate=mean(singleton_set[truth_open!="Unknown"]),known_empty_rate=mean(empty_set[truth_open!="Unknown"]),unknown_empty_rate=mean(empty_set[truth_open=="Unknown"]),unknown_singleton_rate=mean(singleton_set[truth_open=="Unknown"]),unknown_mean_set_size=mean(set_size[truth_open=="Unknown"]),.groups="drop")%>%dplyr::mutate(repeat_id=key$repeat_id,omitted_class=key$omitted_class,model=key$model)
}

metrics_by_repeat<-dplyr::bind_rows(metric_rows)
conformal_metrics_by_repeat<-dplyr::bind_rows(conformal_metric_rows)
forced_assignments<-dplyr::bind_rows(forced_assignment_rows)
metric_long<-metrics_by_repeat%>%tidyr::pivot_longer(cols=c(known_classification_accuracy,known_acceptance_rate,unknown_rejection_rate,accepted_known_accuracy,open_set_accuracy,open_set_balanced_accuracy,unknown_detection_auc,unknown_average_precision,mean_known_confidence,mean_unknown_confidence),names_to="metric",values_to="value")
metric_summary<-summarize_with_ci(metric_long,c("omitted_class","model","metric"))
conformal_long<-conformal_metrics_by_repeat%>%tidyr::pivot_longer(cols=c(known_coverage,known_mean_set_size,known_singleton_rate,known_empty_rate,unknown_empty_rate,unknown_singleton_rate,unknown_mean_set_size),names_to="metric",values_to="value")
conformal_summary<-summarize_with_ci(conformal_long,c("omitted_class","model","alpha","metric"))

write.csv(metrics_by_repeat,file.path(table_dir,"open_set_metrics_by_repeat.csv"),row.names=FALSE)
write.csv(metric_summary,file.path(table_dir,"open_set_metric_summary_95CI.csv"),row.names=FALSE)
write.csv(conformal_metrics_by_repeat,file.path(table_dir,"conformal_metrics_by_repeat.csv"),row.names=FALSE)
write.csv(conformal_summary,file.path(table_dir,"conformal_metric_summary_95CI.csv"),row.names=FALSE)
write.csv(forced_assignments,file.path(table_dir,"unknown_forced_class_assignments.csv"),row.names=FALSE)

# ==============================================================================
# 5. FIGURES
# ==============================================================================
plot_detection<-metric_summary%>%dplyr::filter(metric%in%c("unknown_detection_auc","unknown_rejection_rate","known_acceptance_rate"))%>%dplyr::mutate(metric=factor(metric,levels=c("unknown_detection_auc","unknown_rejection_rate","known_acceptance_rate"),labels=c("Unknown-detection AUC","Unknown rejection","Known acceptance")))
p_detection<-ggplot2::ggplot(plot_detection,ggplot2::aes(x=omitted_class,y=mean,color=model,group=model))+ggplot2::geom_point(position=ggplot2::position_dodge(width=0.35),size=2.2)+ggplot2::geom_errorbar(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),position=ggplot2::position_dodge(width=0.35),width=0.16)+ggplot2::facet_wrap(~metric,ncol=1,scales="free_y")+ggplot2::scale_color_manual(values=c("SVM"="#1F77B4","XGBoost"="#D95F02"))+ggplot2::labs(x="Omitted cell type",y="Estimate",color="Model")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=25,hjust=1),legend.position="top",strip.text=ggplot2::element_text(face="bold"))
save_plot(p_detection,"open_set_detection",8.0,7.2,figure_dir)

plot_conformal<-conformal_summary%>%dplyr::filter(metric%in%c("known_coverage","unknown_empty_rate","known_mean_set_size"))%>%dplyr::mutate(metric=factor(metric,levels=c("known_coverage","unknown_empty_rate","known_mean_set_size"),labels=c("Known-cell coverage","Unknown empty-set rate","Known mean set size")),alpha=factor(alpha))
p_conformal<-ggplot2::ggplot(plot_conformal,ggplot2::aes(x=omitted_class,y=mean,color=model,shape=alpha,group=interaction(model,alpha)))+ggplot2::geom_point(position=ggplot2::position_dodge(width=0.5),size=2)+ggplot2::facet_wrap(~metric,ncol=1,scales="free_y")+ggplot2::scale_color_manual(values=c("SVM"="#1F77B4","XGBoost"="#D95F02"))+ggplot2::labs(x="Omitted cell type",y="Estimate",color="Model",shape=expression(alpha))+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=25,hjust=1),legend.position="top",strip.text=ggplot2::element_text(face="bold"))
save_plot(p_conformal,"conformal_open_set",8.0,7.2,figure_dir)

package_versions<-data.frame(package=required_packages,version=vapply(required_packages,function(x)as.character(utils::packageVersion(x)),character(1)))
write.csv(package_versions,file.path(table_dir,"package_versions.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(results_dir,"sessionInfo.txt"))

cat("\n==============================================================================\n")
cat("OPEN-SET ANALYSIS COMPLETED\n")
cat("==============================================================================\n")
print(metric_summary%>%dplyr::filter(metric%in%c("unknown_detection_auc","unknown_rejection_rate","known_acceptance_rate","open_set_accuracy"))%>%dplyr::arrange(metric,dplyr::desc(mean)))
cat("\nResults saved to:",results_dir,"\n")
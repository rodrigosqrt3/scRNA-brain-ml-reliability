# ==============================================================================
# ERAMIA-RS 2026 - FINAL GROUPED UNCERTAINTY AND OSCR POSTPROCESSING
# ==============================================================================
# This script does not refit any classifier. It performs two final analyses:
# 1. calibration and selective classification using group-blocked OOF results;
# 2. threshold-free open-set evaluation using OSCR curves and AUOSCR.
# If the large open-set CSV is absent, predictions are rebuilt from checkpoints.
# ==============================================================================

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
master_seed<-335705
n_calibration_bins<-10
coverage_grid<-c(1,0.95,0.90,0.85,0.80,0.70,0.60,0.50)
oscr_fpr_grid<-seq(0,1,by=0.01)
oscr_report_fpr<-c(0.05,0.10,0.20)
grouped_analysis_id<-"grouped5x5_c1run_hvg2000_pc100"
open_set_analysis_id<-"open_set_conformal_grouped3x5_hvg2000"
output_analysis_id<-"final_grouped_uncertainty_oscr"

get_script_dir<-function(){
  args<-commandArgs(trailingOnly=FALSE)
  file_arg<-grep("^--file=",args,value=TRUE)
  if(length(file_arg)>0)return(dirname(normalizePath(sub("^--file=","",file_arg[1]))))
  frame_files<-vapply(sys.frames(),function(x)if(is.null(x$ofile))"" else x$ofile,character(1))
  frame_files<-frame_files[nzchar(frame_files)]
  if(length(frame_files)>0)return(dirname(normalizePath(frame_files[length(frame_files)])))
  getwd()
}

resolve_analysis_dir<-function(project_dir,analysis_id){
  candidates<-c(file.path(project_dir,"results",analysis_id),file.path(project_dir,analysis_id))
  existing<-candidates[dir.exists(candidates)]
  if(length(existing)==0)stop(paste("Analysis directory not found for",analysis_id))
  existing[1]
}

project_dir<-get_script_dir()
source(file.path(project_dir,"pipeline_functions.R"))
grouped_dir<-resolve_analysis_dir(project_dir,grouped_analysis_id)
open_set_dir<-resolve_analysis_dir(project_dir,open_set_analysis_id)
results_dir<-file.path(project_dir,"results",output_analysis_id)
table_dir<-file.path(results_dir,"tables")
figure_dir<-file.path(results_dir,"figures")
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)

required_packages<-c("ggplot2","dplyr","tidyr","scales")
missing_packages<-required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)>0)stop(paste("Missing packages:",paste(missing_packages,collapse=", ")))
`%>%`<-dplyr::`%>%`
options(stringsAsFactors=FALSE)
set.seed(master_seed)

configuration<-data.frame(parameter=c("master_seed","n_calibration_bins","coverage_grid","oscr_fpr_grid","oscr_report_fpr","grouped_analysis_id","open_set_analysis_id"),value=c(master_seed,n_calibration_bins,paste(coverage_grid,collapse=";"),paste(oscr_fpr_grid,collapse=";"),paste(oscr_report_fpr,collapse=";"),grouped_analysis_id,open_set_analysis_id))
write.csv(configuration,file.path(table_dir,"analysis_configuration.csv"),row.names=FALSE)

summarize_values<-function(data,group_columns,value_column="value"){
  data%>%dplyr::group_by(dplyr::across(dplyr::all_of(group_columns)))%>%dplyr::summarise(n=dplyr::n(),mean=mean(.data[[value_column]],na.rm=TRUE),sd=stats::sd(.data[[value_column]],na.rm=TRUE),se=sd/sqrt(n),ci_lower=mean-stats::qt(0.975,df=pmax(n-1,1))*se,ci_upper=mean+stats::qt(0.975,df=pmax(n-1,1))*se,.groups="drop")
}

multiclass_brier<-function(probability_matrix,truth,class_levels){
  truth_index<-match(truth,class_levels)
  one_hot<-matrix(0,nrow=length(truth),ncol=length(class_levels))
  one_hot[cbind(seq_along(truth_index),truth_index)]<-1
  mean(rowSums((probability_matrix-one_hot)^2))
}

macro_f1_local<-function(truth,predicted,class_levels){
  confusion<-table(factor(truth,levels=class_levels),factor(predicted,levels=class_levels))
  true_positive<-diag(confusion)
  precision<-safe_divide(true_positive,colSums(confusion))
  recall<-safe_divide(true_positive,rowSums(confusion))
  f1<-safe_divide(2*precision*recall,precision+recall)
  mean(f1,na.rm=TRUE)
}

save_plot_local<-function(plot_object,file_stem,width,height){
  save_plot(plot_object,file_stem,width,height,figure_dir)
}

# ==============================================================================
# 2. GROUP-BLOCKED CALIBRATION AND SELECTIVE CLASSIFICATION
# ==============================================================================
grouped_prediction_path<-file.path(grouped_dir,"tables","out_of_fold_predictions_grouped.csv")
if(!file.exists(grouped_prediction_path))stop(paste("Grouped OOF predictions not found at",grouped_prediction_path))

cat("\n==============================================================================\n")
cat("GROUP-BLOCKED CALIBRATION AND SELECTIVE CLASSIFICATION\n")
cat("==============================================================================\n")
cat("Reading:",grouped_prediction_path,"\n")

grouped_predictions<-read.csv(grouped_prediction_path,check.names=FALSE,stringsAsFactors=FALSE)
grouped_probability_columns<-grep("^prob__",names(grouped_predictions),value=TRUE)
grouped_class_levels<-sort(unique(grouped_predictions$truth))
grouped_expected_columns<-paste0("prob__",make.names(grouped_class_levels,unique=TRUE))
if(!all(grouped_expected_columns%in%grouped_probability_columns))stop("Grouped probability columns do not match the observed classes.")
grouped_probability_columns<-grouped_expected_columns

expected_grouped_rows<-length(unique(grouped_predictions$repeat_id))*length(unique(grouped_predictions$model))*length(unique(grouped_predictions$representation))*length(unique(grouped_predictions$cell_index))
if(nrow(grouped_predictions)!=expected_grouped_rows)stop("Grouped predictions are incomplete.")
if(anyDuplicated(grouped_predictions[,c("repeat_id","model","representation","cell_index")]))stop("Duplicated grouped OOF predictions were detected.")

grouped_calibration_rows<-list()
grouped_metric_rows<-list()
grouped_risk_rows<-list()
grouped_keys<-grouped_predictions%>%dplyr::distinct(repeat_id,model,representation)%>%dplyr::arrange(repeat_id,model,representation)

for(group_id in seq_len(nrow(grouped_keys))){
  key<-grouped_keys[group_id,]
  subset_data<-grouped_predictions%>%dplyr::filter(repeat_id==key$repeat_id,model==key$model,representation==key$representation)%>%dplyr::arrange(cell_index)
  probability_matrix<-as.matrix(subset_data[,grouped_probability_columns,drop=FALSE])
  storage.mode(probability_matrix)<-"double"
  probability_matrix[!is.finite(probability_matrix)]<-0
  probability_matrix<-pmax(probability_matrix,0)
  probability_matrix<-probability_matrix/pmax(rowSums(probability_matrix),1e-15)
  colnames(probability_matrix)<-grouped_class_levels
  predicted_index<-match(subset_data$predicted,grouped_class_levels)
  truth_index<-match(subset_data$truth,grouped_class_levels)
  if(anyNA(predicted_index)|anyNA(truth_index))stop("Unknown grouped class label detected.")
  confidence<-probability_matrix[cbind(seq_len(nrow(probability_matrix)),predicted_index)]
  correct<-as.integer(subset_data$predicted==subset_data$truth)
  probability_argmax<-grouped_class_levels[max.col(probability_matrix,ties.method="first")]
  ordered_probability<-t(apply(probability_matrix,1,sort,decreasing=TRUE))
  margin<-ordered_probability[,1]-ordered_probability[,2]
  entropy<--rowSums(probability_matrix*log(pmax(probability_matrix,1e-15)))/log(length(grouped_class_levels))
  calibration_bin<-cut(confidence,breaks=seq(0,1,length.out=n_calibration_bins+1),include.lowest=TRUE,right=TRUE,labels=FALSE)
  calibration_data<-data.frame(calibration_bin=seq_len(n_calibration_bins))%>%dplyr::left_join(data.frame(calibration_bin=calibration_bin,confidence=confidence,correct=correct)%>%dplyr::group_by(calibration_bin)%>%dplyr::summarise(n=dplyr::n(),mean_confidence=mean(confidence),observed_accuracy=mean(correct),absolute_gap=abs(mean_confidence-observed_accuracy),.groups="drop"),by="calibration_bin")%>%dplyr::mutate(repeat_id=key$repeat_id,model=key$model,representation=key$representation)
  grouped_calibration_rows[[group_id]]<-calibration_data
  valid_calibration<-calibration_data%>%dplyr::filter(!is.na(n))
  ece<-sum(valid_calibration$n*valid_calibration$absolute_gap)/sum(valid_calibration$n)
  observed_probability<-probability_matrix[cbind(seq_along(truth_index),truth_index)]
  grouped_metric_rows[[group_id]]<-data.frame(repeat_id=key$repeat_id,model=key$model,representation=key$representation,accuracy=mean(correct),mean_confidence=mean(confidence),confidence_bias=mean(confidence)-mean(correct),ece=ece,maximum_calibration_error=max(valid_calibration$absolute_gap),log_loss=-mean(log(pmax(observed_probability,1e-15))),brier_score=multiclass_brier(probability_matrix,subset_data$truth,grouped_class_levels),mean_entropy=mean(entropy),mean_margin=mean(margin),native_argmax_agreement=mean(subset_data$predicted==probability_argmax))
  grouped_risk_rows[[group_id]]<-dplyr::bind_rows(lapply(coverage_grid,function(target_coverage){
    retained_n<-max(1,ceiling(target_coverage*nrow(subset_data)))
    retained_index<-order(confidence,decreasing=TRUE)[seq_len(retained_n)]
    data.frame(repeat_id=key$repeat_id,model=key$model,representation=key$representation,target_coverage=target_coverage,observed_coverage=retained_n/nrow(subset_data),retained_n=retained_n,selective_accuracy=mean(correct[retained_index]),selective_error=1-mean(correct[retained_index]),selective_macro_f1=macro_f1_local(subset_data$truth[retained_index],subset_data$predicted[retained_index],grouped_class_levels),mean_retained_confidence=mean(confidence[retained_index]))
  }))
}

grouped_calibration<-dplyr::bind_rows(grouped_calibration_rows)
grouped_metrics<-dplyr::bind_rows(grouped_metric_rows)
grouped_risk<-dplyr::bind_rows(grouped_risk_rows)
grouped_metric_long<-grouped_metrics%>%tidyr::pivot_longer(cols=c(accuracy,mean_confidence,confidence_bias,ece,maximum_calibration_error,log_loss,brier_score,mean_entropy,mean_margin,native_argmax_agreement),names_to="metric",values_to="value")
grouped_metric_summary<-summarize_values(grouped_metric_long,c("model","representation","metric"))
grouped_risk_long<-grouped_risk%>%tidyr::pivot_longer(cols=c(selective_accuracy,selective_error,selective_macro_f1,mean_retained_confidence),names_to="metric",values_to="value")
grouped_risk_summary<-summarize_values(grouped_risk_long,c("model","representation","target_coverage","metric"))
grouped_reliability<-grouped_calibration%>%dplyr::filter(!is.na(n))%>%dplyr::group_by(model,representation,calibration_bin)%>%dplyr::summarise(total_n=sum(n),mean_confidence=stats::weighted.mean(mean_confidence,n),observed_accuracy=stats::weighted.mean(observed_accuracy,n),absolute_gap=abs(mean_confidence-observed_accuracy),.groups="drop")

write.csv(grouped_calibration,file.path(table_dir,"grouped_calibration_by_repeat_and_bin.csv"),row.names=FALSE)
write.csv(grouped_metrics,file.path(table_dir,"grouped_calibration_metrics_by_repeat.csv"),row.names=FALSE)
write.csv(grouped_metric_summary,file.path(table_dir,"grouped_calibration_summary_95CI.csv"),row.names=FALSE)
write.csv(grouped_risk,file.path(table_dir,"grouped_risk_coverage_by_repeat.csv"),row.names=FALSE)
write.csv(grouped_risk_summary,file.path(table_dir,"grouped_risk_coverage_summary_95CI.csv"),row.names=FALSE)
write.csv(grouped_reliability,file.path(table_dir,"grouped_reliability_summary.csv"),row.names=FALSE)

model_order<-c("XGBoost","Random Forest","SVM","Naive Bayes","KNN")
plot_grouped_risk<-grouped_risk_summary%>%dplyr::filter(metric=="selective_accuracy")%>%dplyr::mutate(model=factor(model,levels=model_order))
p_grouped_risk<-ggplot2::ggplot(plot_grouped_risk,ggplot2::aes(x=target_coverage,y=mean,color=model,group=model))+ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower,ymax=ci_upper,fill=model),alpha=0.08,color=NA)+ggplot2::geom_line(linewidth=0.8)+ggplot2::geom_point(size=1.8)+ggplot2::facet_wrap(~representation,ncol=2)+ggplot2::scale_x_reverse(breaks=sort(unique(coverage_grid),decreasing=TRUE),labels=scales::percent_format())+ggplot2::scale_y_continuous(labels=scales::percent_format(accuracy=0.1))+ggplot2::labs(x="Retained coverage",y="Group-blocked accuracy after abstention",color="Model",fill="Model")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom",strip.text=ggplot2::element_text(face="bold"))
save_plot_local(p_grouped_risk,"grouped_selective_accuracy_coverage",8.2,4.5)

plot_grouped_reliability<-grouped_reliability%>%dplyr::mutate(model=factor(model,levels=model_order))
p_grouped_reliability<-ggplot2::ggplot(plot_grouped_reliability,ggplot2::aes(x=mean_confidence,y=observed_accuracy,color=model,group=model,size=total_n))+ggplot2::geom_abline(slope=1,intercept=0,linetype="dashed",color="#777777")+ggplot2::geom_line(linewidth=0.6)+ggplot2::geom_point(alpha=0.9)+ggplot2::facet_wrap(~representation,ncol=2)+ggplot2::coord_equal(xlim=c(0,1),ylim=c(0,1))+ggplot2::scale_size_continuous(range=c(1.5,5),guide="none")+ggplot2::labs(x="Mean predicted confidence",y="Observed group-blocked accuracy",color="Model")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom",strip.text=ggplot2::element_text(face="bold"))
save_plot_local(p_grouped_reliability,"grouped_reliability_diagram",8.2,4.5)

# ==============================================================================
# 3. OPEN-SET CLASSIFICATION RATE AND AUOSCR
# ==============================================================================
open_prediction_path<-file.path(open_set_dir,"tables","open_set_predictions.csv")
if(file.exists(open_prediction_path)){
  cat("\nReading open-set predictions:",open_prediction_path,"\n")
  open_predictions<-read.csv(open_prediction_path,check.names=FALSE,stringsAsFactors=FALSE)
}else{
  checkpoint_dir<-file.path(open_set_dir,"checkpoints")
  checkpoint_files<-list.files(checkpoint_dir,pattern="\\.rds$",full.names=TRUE)
  if(length(checkpoint_files)==0)stop("Neither open_set_predictions.csv nor open-set checkpoints were found.")
  cat("\nReconstructing open-set predictions from",length(checkpoint_files),"checkpoints.\n")
  open_results<-lapply(checkpoint_files,readRDS)
  open_predictions<-dplyr::bind_rows(lapply(open_results,"[[","predictions"))
  rm(open_results)
  gc(verbose=FALSE)
}

required_open_columns<-c("cell_index","repeat_id","omitted_class","model","truth","truth_open","predicted","max_probability")
if(!all(required_open_columns%in%names(open_predictions)))stop("Open-set prediction columns are incomplete.")
if(anyDuplicated(open_predictions[,c("repeat_id","omitted_class","model","cell_index")]))stop("Duplicated open-set predictions were detected.")
expected_open_rows<-length(unique(open_predictions$repeat_id))*length(unique(open_predictions$omitted_class))*length(unique(open_predictions$model))*length(unique(open_predictions$cell_index))
if(nrow(open_predictions)!=expected_open_rows)stop(paste("Open-set predictions are incomplete: expected",expected_open_rows,"rows but found",nrow(open_predictions)))

build_oscr_curve<-function(data){
  known<-data$truth_open!="Unknown"
  unknown<-!known
  n_known<-sum(known)
  n_unknown<-sum(unknown)
  if(n_known==0|n_unknown==0)stop("OSCR requires known and unknown observations.")
  event_data<-data.frame(confidence=data$max_probability,correct_known=as.integer(known&data$predicted==data$truth),unknown_accept=as.integer(unknown))%>%dplyr::group_by(confidence)%>%dplyr::summarise(correct_known=sum(correct_known),unknown_accept=sum(unknown_accept),.groups="drop")%>%dplyr::arrange(dplyr::desc(confidence))%>%dplyr::mutate(ccr=cumsum(correct_known)/n_known,fpr=cumsum(unknown_accept)/n_unknown)
  dplyr::bind_rows(data.frame(confidence=Inf,correct_known=0,unknown_accept=0,ccr=0,fpr=0),event_data)
}

calculate_auoscr<-function(curve){
  ordering<-order(curve$fpr,curve$ccr)
  x<-curve$fpr[ordering]
  y<-curve$ccr[ordering]
  sum(diff(x)*(head(y,-1)+tail(y,-1))/2)
}

oscr_curve_rows<-list()
oscr_metric_rows<-list()
oscr_grid_rows<-list()
open_keys<-open_predictions%>%dplyr::distinct(repeat_id,omitted_class,model)%>%dplyr::arrange(repeat_id,omitted_class,model)

for(group_id in seq_len(nrow(open_keys))){
  key<-open_keys[group_id,]
  subset_data<-open_predictions%>%dplyr::filter(repeat_id==key$repeat_id,omitted_class==key$omitted_class,model==key$model)
  curve<-build_oscr_curve(subset_data)
  curve$repeat_id<-key$repeat_id
  curve$omitted_class<-key$omitted_class
  curve$model<-key$model
  oscr_curve_rows[[group_id]]<-curve
  known<-subset_data$truth_open!="Unknown"
  report_values<-vapply(oscr_report_fpr,function(fpr_limit)max(curve$ccr[curve$fpr<=fpr_limit],na.rm=TRUE),numeric(1))
  metric_row<-data.frame(repeat_id=key$repeat_id,omitted_class=key$omitted_class,model=key$model,auoscr=calculate_auoscr(curve),known_closed_set_accuracy=mean(subset_data$predicted[known]==subset_data$truth[known]))
  for(fpr_id in seq_along(oscr_report_fpr))metric_row[[paste0("ccr_at_fpr_",sprintf("%02d",round(100*oscr_report_fpr[fpr_id])))]]<-report_values[fpr_id]
  oscr_metric_rows[[group_id]]<-metric_row
  oscr_grid_rows[[group_id]]<-data.frame(repeat_id=key$repeat_id,omitted_class=key$omitted_class,model=key$model,fpr=oscr_fpr_grid,ccr=vapply(oscr_fpr_grid,function(fpr_limit)max(curve$ccr[curve$fpr<=fpr_limit],na.rm=TRUE),numeric(1)))
}

oscr_curves<-dplyr::bind_rows(oscr_curve_rows)
oscr_metrics<-dplyr::bind_rows(oscr_metric_rows)
oscr_grid<-dplyr::bind_rows(oscr_grid_rows)
oscr_metric_columns<-setdiff(names(oscr_metrics),c("repeat_id","omitted_class","model"))
oscr_metric_long<-oscr_metrics%>%tidyr::pivot_longer(cols=dplyr::all_of(oscr_metric_columns),names_to="metric",values_to="value")
oscr_metric_summary<-summarize_values(oscr_metric_long,c("omitted_class","model","metric"))
oscr_overall_by_repeat<-oscr_metric_long%>%dplyr::group_by(repeat_id,model,metric)%>%dplyr::summarise(value=mean(value),.groups="drop")
oscr_overall_summary<-summarize_values(oscr_overall_by_repeat,c("model","metric"))
oscr_grid_summary<-summarize_values(oscr_grid,c("omitted_class","model","fpr"),"ccr")

write.csv(oscr_metrics,file.path(table_dir,"oscr_metrics_by_repeat.csv"),row.names=FALSE)
write.csv(oscr_metric_summary,file.path(table_dir,"oscr_metric_summary_95CI.csv"),row.names=FALSE)
write.csv(oscr_overall_by_repeat,file.path(table_dir,"oscr_overall_by_repeat.csv"),row.names=FALSE)
write.csv(oscr_overall_summary,file.path(table_dir,"oscr_overall_summary_95CI.csv"),row.names=FALSE)
write.csv(oscr_grid_summary,file.path(table_dir,"oscr_curve_summary_95CI.csv"),row.names=FALSE)

plot_oscr<-oscr_grid_summary%>%dplyr::mutate(model=factor(model,levels=c("SVM","XGBoost")))
p_oscr<-ggplot2::ggplot(plot_oscr,ggplot2::aes(x=fpr,y=mean,color=model,fill=model,group=model))+ggplot2::geom_ribbon(ggplot2::aes(ymin=pmax(ci_lower,0),ymax=pmin(ci_upper,1)),alpha=0.10,color=NA)+ggplot2::geom_line(linewidth=0.8)+ggplot2::facet_wrap(~omitted_class,ncol=3)+ggplot2::scale_x_continuous(labels=scales::percent_format(),breaks=seq(0,1,by=0.25))+ggplot2::scale_y_continuous(labels=scales::percent_format(),limits=c(0,1))+ggplot2::scale_color_manual(values=c("SVM"="#1F77B4","XGBoost"="#D95F02"))+ggplot2::scale_fill_manual(values=c("SVM"="#1F77B4","XGBoost"="#D95F02"))+ggplot2::labs(x="False-positive rate for unknown cells",y="Correct classification rate for known cells",color="Model",fill="Model")+ggplot2::theme_minimal(base_size=9)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom",strip.text=ggplot2::element_text(face="bold"))
save_plot_local(p_oscr,"open_set_oscr_curves",8.2,7.0)

# ==============================================================================
# 4. REPRODUCIBILITY AND CONSOLE SUMMARY
# ==============================================================================
package_versions<-data.frame(package=required_packages,version=vapply(required_packages,function(x)as.character(utils::packageVersion(x)),character(1)))
write.csv(package_versions,file.path(table_dir,"package_versions.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(results_dir,"sessionInfo.txt"))

cat("\n==============================================================================\n")
cat("FINAL POSTPROCESSING COMPLETED\n")
cat("==============================================================================\n")
cat("Grouped prediction rows:",nrow(grouped_predictions),"\n")
cat("Grouped model-representation-repeat groups:",nrow(grouped_keys),"\n")
cat("Open-set prediction rows:",nrow(open_predictions),"\n")
cat("Open-set class-model-repeat groups:",nrow(open_keys),"\n\n")
cat("Group-blocked calibration summary:\n")
print(grouped_metric_summary%>%dplyr::filter(metric%in%c("ece","log_loss","brier_score"))%>%dplyr::arrange(metric,mean))
cat("\nOverall OSCR summary:\n")
print(oscr_overall_summary%>%dplyr::arrange(metric,dplyr::desc(mean)))
cat("\nUpload the complete folder:\n")
cat("  ",results_dir,"\n")
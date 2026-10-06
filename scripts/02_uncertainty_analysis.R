# ==============================================================================
# ERAMIA-RS 2026 - CALIBRATION AND SELECTIVE CLASSIFICATION
# ==============================================================================
# This script reuses the completed out-of-fold predictions. It does not refit
# any classifier. Uncertainty estimates are calculated independently for each
# complete repeated-CV prediction set and are summarized across repeats.
# ==============================================================================

master_seed<-335705
n_calibration_bins<-10
coverage_grid<-c(1,0.95,0.90,0.85,0.80,0.70,0.60,0.50)
primary_analysis_id<-"repeated5x5_hvg2000_pc100"

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
primary_table_dir<-file.path(project_dir,"results",primary_analysis_id,"tables")
if(!dir.exists(primary_table_dir)&&dir.exists(file.path(project_dir,"tables")))primary_table_dir<-file.path(project_dir,"tables")
prediction_path<-file.path(primary_table_dir,"out_of_fold_predictions.csv")
if(!file.exists(prediction_path))stop(paste("Out-of-fold predictions not found at",prediction_path))

results_dir<-file.path(project_dir,"results","uncertainty_analysis")
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

safe_divide_local<-function(numerator,denominator){
  result<-numerator/denominator
  result[!is.finite(result)]<-NA_real_
  result
}

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
  precision<-safe_divide_local(true_positive,colSums(confusion))
  recall<-safe_divide_local(true_positive,rowSums(confusion))
  f1<-safe_divide_local(2*precision*recall,precision+recall)
  mean(f1,na.rm=TRUE)
}

cat("\n==============================================================================\n")
cat("CALIBRATION AND SELECTIVE CLASSIFICATION\n")
cat("==============================================================================\n")
cat("Reading:",prediction_path,"\n")

predictions<-read.csv(prediction_path,check.names=FALSE,stringsAsFactors=FALSE)
probability_columns<-grep("^prob__",names(predictions),value=TRUE)
if(length(probability_columns)<2)stop("Probability columns were not found.")
class_levels<-sort(unique(predictions$truth))
expected_probability_columns<-paste0("prob__",make.names(class_levels,unique=TRUE))
if(!all(expected_probability_columns%in%probability_columns))stop("Probability columns do not match the observed classes.")
probability_columns<-expected_probability_columns

enriched_rows<-list()
calibration_rows<-list()
calibration_metric_rows<-list()
risk_coverage_rows<-list()
class_confidence_rows<-list()
group_keys<-predictions%>%dplyr::distinct(repeat_id,model,representation)%>%dplyr::arrange(repeat_id,model,representation)

for(group_id in seq_len(nrow(group_keys))){
  key<-group_keys[group_id,]
  subset_data<-predictions%>%dplyr::filter(repeat_id==key$repeat_id,model==key$model,representation==key$representation)%>%dplyr::arrange(cell_index)
  probability_matrix<-as.matrix(subset_data[,probability_columns,drop=FALSE])
  storage.mode(probability_matrix)<-"double"
  probability_matrix[!is.finite(probability_matrix)]<-0
  probability_matrix<-probability_matrix/pmax(rowSums(probability_matrix),1e-15)
  probability_matrix<-pmin(pmax(probability_matrix,0),1)
  colnames(probability_matrix)<-class_levels
  predicted_index<-match(subset_data$predicted,class_levels)
  truth_index<-match(subset_data$truth,class_levels)
  prediction_confidence<-probability_matrix[cbind(seq_len(nrow(probability_matrix)),predicted_index)]
  probability_argmax<-class_levels[max.col(probability_matrix,ties.method="first")]
  ordered_probabilities<-t(apply(probability_matrix,1,sort,decreasing=TRUE))
  probability_margin<-ordered_probabilities[,1]-ordered_probabilities[,2]
  normalized_entropy<--rowSums(probability_matrix*log(pmax(probability_matrix,1e-15)))/log(length(class_levels))
  correct<-as.integer(subset_data$truth==subset_data$predicted)
  calibration_bin<-cut(prediction_confidence,breaks=seq(0,1,length.out=n_calibration_bins+1),include.lowest=TRUE,right=TRUE,labels=FALSE)
  enriched<-subset_data
  enriched$correct<-correct
  enriched$probability_argmax<-probability_argmax
  enriched$native_argmax_agreement<-as.integer(subset_data$predicted==probability_argmax)
  enriched$prediction_confidence<-prediction_confidence
  enriched$probability_margin<-probability_margin
  enriched$normalized_entropy<-normalized_entropy
  enriched$calibration_bin<-calibration_bin
  enriched_rows[[group_id]]<-enriched
  
  calibration_data<-data.frame(calibration_bin=seq_len(n_calibration_bins))%>%dplyr::left_join(data.frame(calibration_bin=calibration_bin,confidence=prediction_confidence,correct=correct)%>%dplyr::group_by(calibration_bin)%>%dplyr::summarise(n=dplyr::n(),mean_confidence=mean(confidence),observed_accuracy=mean(correct),absolute_gap=abs(mean_confidence-observed_accuracy),.groups="drop"),by="calibration_bin")%>%dplyr::mutate(repeat_id=key$repeat_id,model=key$model,representation=key$representation)
  calibration_rows[[group_id]]<-calibration_data
  valid_calibration<-calibration_data%>%dplyr::filter(!is.na(n))
  ece<-sum(valid_calibration$n*valid_calibration$absolute_gap)/sum(valid_calibration$n)
  observed_probability<-probability_matrix[cbind(seq_along(truth_index),truth_index)]
  calibration_metric_rows[[group_id]]<-data.frame(repeat_id=key$repeat_id,model=key$model,representation=key$representation,accuracy=mean(correct),mean_confidence=mean(prediction_confidence),confidence_bias=mean(prediction_confidence)-mean(correct),ece=ece,maximum_calibration_error=max(valid_calibration$absolute_gap),log_loss=-mean(log(pmax(observed_probability,1e-15))),brier_score=multiclass_brier(probability_matrix,subset_data$truth,class_levels),mean_entropy=mean(normalized_entropy),mean_margin=mean(probability_margin),native_argmax_agreement=mean(subset_data$predicted==probability_argmax))
  
  risk_rows_group<-lapply(coverage_grid,function(target_coverage){
    retained_n<-max(1,ceiling(target_coverage*nrow(subset_data)))
    retained_index<-order(prediction_confidence,decreasing=TRUE)[seq_len(retained_n)]
    data.frame(repeat_id=key$repeat_id,model=key$model,representation=key$representation,target_coverage=target_coverage,observed_coverage=retained_n/nrow(subset_data),retained_n=retained_n,selective_accuracy=mean(correct[retained_index]),selective_error=1-mean(correct[retained_index]),selective_macro_f1=macro_f1_local(subset_data$truth[retained_index],subset_data$predicted[retained_index],class_levels),mean_retained_confidence=mean(prediction_confidence[retained_index]))
  })
  risk_coverage_rows[[group_id]]<-dplyr::bind_rows(risk_rows_group)
  
  class_confidence_rows[[group_id]]<-data.frame(truth=subset_data$truth,correct=correct,confidence=prediction_confidence,entropy=normalized_entropy)%>%dplyr::group_by(truth)%>%dplyr::summarise(n=dplyr::n(),accuracy=mean(correct),mean_confidence=mean(confidence),median_confidence=median(confidence),mean_entropy=mean(entropy),.groups="drop")%>%dplyr::mutate(repeat_id=key$repeat_id,model=key$model,representation=key$representation)
}

enriched_predictions<-dplyr::bind_rows(enriched_rows)
calibration_by_bin<-dplyr::bind_rows(calibration_rows)
calibration_metrics<-dplyr::bind_rows(calibration_metric_rows)
risk_coverage<-dplyr::bind_rows(risk_coverage_rows)
class_confidence<-dplyr::bind_rows(class_confidence_rows)

calibration_long<-calibration_metrics%>%tidyr::pivot_longer(cols=c(accuracy,mean_confidence,confidence_bias,ece,maximum_calibration_error,log_loss,brier_score,mean_entropy,mean_margin,native_argmax_agreement),names_to="metric",values_to="value")
calibration_summary<-summarize_values(calibration_long,c("model","representation","metric"))
risk_coverage_long<-risk_coverage%>%tidyr::pivot_longer(cols=c(selective_accuracy,selective_error,selective_macro_f1,mean_retained_confidence),names_to="metric",values_to="value")
risk_coverage_summary<-summarize_values(risk_coverage_long,c("model","representation","target_coverage","metric"))
class_confidence_long<-class_confidence%>%tidyr::pivot_longer(cols=c(accuracy,mean_confidence,median_confidence,mean_entropy),names_to="metric",values_to="value")
class_confidence_summary<-summarize_values(class_confidence_long,c("model","representation","truth","metric"))
reliability_summary<-calibration_by_bin%>%dplyr::filter(!is.na(n))%>%dplyr::group_by(model,representation,calibration_bin)%>%dplyr::summarise(total_n=sum(n),mean_confidence=stats::weighted.mean(mean_confidence,n),observed_accuracy=stats::weighted.mean(observed_accuracy,n),absolute_gap=abs(mean_confidence-observed_accuracy),.groups="drop")

write.csv(enriched_predictions,file.path(table_dir,"predictions_with_uncertainty.csv"),row.names=FALSE)
write.csv(calibration_by_bin,file.path(table_dir,"calibration_by_repeat_and_bin.csv"),row.names=FALSE)
write.csv(calibration_metrics,file.path(table_dir,"calibration_metrics_by_repeat.csv"),row.names=FALSE)
write.csv(calibration_summary,file.path(table_dir,"calibration_summary_95CI.csv"),row.names=FALSE)
write.csv(risk_coverage,file.path(table_dir,"risk_coverage_by_repeat.csv"),row.names=FALSE)
write.csv(risk_coverage_summary,file.path(table_dir,"risk_coverage_summary_95CI.csv"),row.names=FALSE)
write.csv(class_confidence,file.path(table_dir,"class_confidence_by_repeat.csv"),row.names=FALSE)
write.csv(class_confidence_summary,file.path(table_dir,"class_confidence_summary_95CI.csv"),row.names=FALSE)
write.csv(reliability_summary,file.path(table_dir,"reliability_summary.csv"),row.names=FALSE)

model_order<-c("XGBoost","Random Forest","SVM","Naive Bayes","KNN")
plot_risk<-risk_coverage_summary%>%dplyr::filter(metric=="selective_accuracy")%>%dplyr::mutate(model=factor(model,levels=model_order))
p_risk<-ggplot2::ggplot(plot_risk,ggplot2::aes(x=target_coverage,y=mean,color=model,group=model))+ggplot2::geom_line(linewidth=0.8)+ggplot2::geom_point(size=1.8)+ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower,ymax=ci_upper,fill=model),alpha=0.08,color=NA)+ggplot2::facet_wrap(~representation,ncol=2)+ggplot2::scale_x_reverse(breaks=sort(unique(coverage_grid),decreasing=TRUE),labels=scales::percent_format())+ggplot2::scale_y_continuous(labels=scales::percent_format(accuracy=0.1))+ggplot2::labs(x="Retained coverage",y="Accuracy after abstention",color="Model",fill="Model")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom",strip.text=ggplot2::element_text(face="bold"))
ggplot2::ggsave(file.path(figure_dir,"selective_accuracy_coverage.pdf"),p_risk,width=8.2,height=4.5,device=if(capabilities("cairo"))grDevices::cairo_pdf else "pdf")
ggplot2::ggsave(file.path(figure_dir,"selective_accuracy_coverage.png"),p_risk,width=8.2,height=4.5,dpi=300,bg="white")

plot_reliability<-reliability_summary%>%dplyr::mutate(model=factor(model,levels=model_order))
p_reliability<-ggplot2::ggplot(plot_reliability,ggplot2::aes(x=mean_confidence,y=observed_accuracy,color=model,group=model,size=total_n))+ggplot2::geom_abline(slope=1,intercept=0,linetype="dashed",color="#777777")+ggplot2::geom_line(linewidth=0.6)+ggplot2::geom_point(alpha=0.9)+ggplot2::facet_wrap(~representation,ncol=2)+ggplot2::coord_equal(xlim=c(0,1),ylim=c(0,1))+ggplot2::scale_size_continuous(range=c(1.5,5),guide="none")+ggplot2::labs(x="Mean predicted confidence",y="Observed accuracy",color="Model")+ggplot2::theme_minimal(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom",strip.text=ggplot2::element_text(face="bold"))
ggplot2::ggsave(file.path(figure_dir,"reliability_diagram.pdf"),p_reliability,width=8.2,height=4.5,device=if(capabilities("cairo"))grDevices::cairo_pdf else "pdf")
ggplot2::ggsave(file.path(figure_dir,"reliability_diagram.png"),p_reliability,width=8.2,height=4.5,dpi=300,bg="white")

capture.output(sessionInfo(),file=file.path(results_dir,"sessionInfo.txt"))

cat("Prediction rows:",nrow(predictions),"\n")
cat("Completed model-representation-repeat groups:",nrow(group_keys),"\n\n")
cat("Calibration summary:\n")
print(calibration_summary%>%dplyr::filter(metric%in%c("ece","log_loss","brier_score"))%>%dplyr::arrange(metric,mean))
cat("\nUncertainty analysis completed. Upload the complete folder:\n")
cat("  ",results_dir,"\n")
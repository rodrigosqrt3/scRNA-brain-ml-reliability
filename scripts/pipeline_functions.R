# ==============================================================================
# FUNCTIONS FOR THE ERAMIA-RS REPEATED CROSS-VALIDATION PIPELINE
# ==============================================================================

safe_divide<-function(numerator,denominator){
  result<-numerator/denominator
  result[!is.finite(result)]<-NA_real_
  result
}

normalize_probabilities<-function(prob,class_levels){
  prob<-as.matrix(prob)
  aligned<-matrix(0,nrow=nrow(prob),ncol=length(class_levels),dimnames=list(NULL,class_levels))
  common<-intersect(colnames(prob),class_levels)
  if(length(common)>0)aligned[,common]<-prob[,common,drop=FALSE]
  aligned[!is.finite(aligned)]<-0
  aligned[aligned<0]<-0
  row_total<-rowSums(aligned)
  zero_rows<-row_total<=0
  if(any(zero_rows))aligned[zero_rows,]<-1/length(class_levels)
  if(any(!zero_rows))aligned[!zero_rows,]<-aligned[!zero_rows,,drop=FALSE]/row_total[!zero_rows]
  aligned
}

create_repeated_stratified_folds<-function(y,k,repeats,seed){
  y<-factor(y)
  output<-vector("list",repeats)
  for(repeat_id in seq_len(repeats)){
    set.seed(seed+repeat_id)
    fold_id<-integer(length(y))
    for(class_name in levels(y)){
      indices<-sample(which(y==class_name))
      fold_id[indices]<-rep(seq_len(k),length.out=length(indices))
    }
    output[[repeat_id]]<-data.frame(cell_index=seq_along(y),repeat_id=repeat_id,fold=fold_id)
  }
  dplyr::bind_rows(output)
}

create_repeated_group_stratified_folds<-function(y,groups,k,repeats,seed,n_restarts=250){
  y<-factor(y)
  groups<-as.character(groups)
  if(length(y)!=length(groups))stop("Labels and groups must have the same length.")
  if(any(is.na(groups)|!nzchar(groups)))stop("Grouping variable contains missing or empty values.")
  group_class<-table(groups,y)
  group_names<-rownames(group_class)
  if(length(group_names)<k)stop("The number of groups is smaller than the requested number of folds.")
  groups_per_class<-colSums(group_class>0)
  if(any(groups_per_class<k))stop(paste("At least one class occurs in fewer than",k,"groups."))
  group_sizes<-rowSums(group_class)
  target_class<-colSums(group_class)/k
  target_size<-sum(group_class)/k
  evaluate_assignment<-function(assignment){
    fold_class<-matrix(0,nrow=k,ncol=ncol(group_class))
    for(group_id in seq_along(assignment)){
      fold_id<-assignment[group_id]
      if(!is.na(fold_id))fold_class[fold_id,]<-fold_class[fold_id,]+group_class[group_id,]
    }
    fold_size<-rowSums(fold_class)
    class_loss<-sum(((fold_class-matrix(target_class,nrow=k,ncol=ncol(group_class),byrow=TRUE))/pmax(matrix(target_class,nrow=k,ncol=ncol(group_class),byrow=TRUE),1))^2)
    size_loss<-sum(((fold_size-target_size)/pmax(target_size,1))^2)
    empty_penalty<-sum(fold_class==0)*25
    class_loss+0.25*size_loss+empty_penalty
  }
  output<-vector("list",repeats)
  assignment_output<-vector("list",repeats)
  for(repeat_id in seq_len(repeats)){
    best_assignment<-NULL
    best_score<-Inf
    for(restart_id in seq_len(n_restarts)){
      set.seed(seed+repeat_id*100000+restart_id)
      ordering<-order(-group_sizes+stats::runif(length(group_sizes),-1e-6,1e-6))
      assignment<-rep(NA_integer_,length(group_names))
      initial_folds<-sample(seq_len(k))
      assignment[ordering[seq_len(k)]]<-initial_folds
      if(length(ordering)>k){
        for(group_id in ordering[(k+1):length(ordering)]){
          candidate_scores<-vapply(seq_len(k),function(fold_id){
            candidate_assignment<-assignment
            candidate_assignment[group_id]<-fold_id
            evaluate_assignment(candidate_assignment)
          },numeric(1))
          best_folds<-which(candidate_scores==min(candidate_scores))
          selected_position<-if(length(best_folds)==1)1 else sample.int(length(best_folds),1)
          assignment[group_id]<-best_folds[selected_position]
        }
      }
      score<-evaluate_assignment(assignment)
      if(score<best_score){
        best_score<-score
        best_assignment<-assignment
      }
    }
    fold_class<-matrix(0,nrow=k,ncol=ncol(group_class),dimnames=list(seq_len(k),colnames(group_class)))
    for(group_id in seq_along(best_assignment))fold_class[best_assignment[group_id],]<-fold_class[best_assignment[group_id],]+group_class[group_id,]
    if(any(fold_class==0))stop(paste("Unable to create class-complete grouped folds for repeat",repeat_id,"."))
    group_assignment<-data.frame(group_id=group_names,repeat_id=repeat_id,fold=best_assignment,group_size=as.numeric(group_sizes),optimization_score=best_score)
    cell_fold<-best_assignment[match(groups,group_names)]
    output[[repeat_id]]<-data.frame(cell_index=seq_along(y),repeat_id=repeat_id,fold=cell_fold,group_id=groups)
    assignment_output[[repeat_id]]<-group_assignment
  }
  list(cell_folds=dplyr::bind_rows(output),group_assignments=dplyr::bind_rows(assignment_output),group_class_matrix=group_class)
}

sparse_row_variances<-function(x){
  n<-ncol(x)
  means<-Matrix::rowMeans(x)
  squared_means<-Matrix::rowMeans(x^2)
  variances<-(squared_means-means^2)*n/(n-1)
  variances[!is.finite(variances)]<-0
  pmax(variances,0)
}

prepare_fold_data<-function(logcounts,train_index,test_index,min_cells,n_features,n_components,seed){
  detected_in_training<-Matrix::rowSums(logcounts[,train_index,drop=FALSE]>0)
  eligible<-which(detected_in_training>=min_cells)
  if(length(eligible)<n_features)stop("Fewer eligible genes than requested features in this fold.")
  training_sparse<-logcounts[eligible,train_index,drop=FALSE]
  gene_variance<-sparse_row_variances(training_sparse)
  selected_local<-order(gene_variance,decreasing=TRUE)[seq_len(n_features)]
  selected_genes<-rownames(training_sparse)[selected_local]
  x_train<-as.matrix(Matrix::t(training_sparse[selected_local,,drop=FALSE]))
  x_test<-as.matrix(Matrix::t(logcounts[selected_genes,test_index,drop=FALSE]))
  centers<-matrixStats::colMeans2(x_train)
  scales<-matrixStats::colSds(x_train,center=centers)
  valid<-is.finite(scales)&scales>0
  x_train<-x_train[,valid,drop=FALSE]
  x_test<-x_test[,valid,drop=FALSE]
  selected_genes<-selected_genes[valid]
  centers<-centers[valid]
  scales<-scales[valid]
  x_train<-sweep(sweep(x_train,2,centers,"-"),2,scales,"/")
  x_test<-sweep(sweep(x_test,2,centers,"-"),2,scales,"/")
  colnames(x_train)<-selected_genes
  colnames(x_test)<-selected_genes
  rank_pca<-min(n_components,nrow(x_train)-1,ncol(x_train)-1)
  set.seed(seed)
  pca_fit<-irlba::prcomp_irlba(x_train,n=rank_pca,center=FALSE,scale.=FALSE)
  pca_train<-pca_fit$x[,seq_len(rank_pca),drop=FALSE]
  pca_test<-x_test%*%pca_fit$rotation[,seq_len(rank_pca),drop=FALSE]
  pc_names<-sprintf("PC%03d",seq_len(rank_pca))
  colnames(pca_train)<-pc_names
  colnames(pca_test)<-pc_names
  list(Genes=list(train=x_train,test=x_test,features=selected_genes),`100 PCs`=list(train=pca_train,test=pca_test,features=pc_names),metadata=data.frame(eligible_genes=length(eligible),selected_genes=length(selected_genes),pca_components=rank_pca))
}

extract_top_features<-function(values,feature_names,model,representation,top_n){
  values<-as.numeric(values)
  valid<-is.finite(values)&values>=0
  if(!any(valid))return(data.frame())
  values<-values[valid]
  feature_names<-feature_names[valid]
  ordering<-order(values,decreasing=TRUE)[seq_len(min(top_n,length(values)))]
  selected_values<-values[ordering]
  maximum<-max(selected_values)
  data.frame(model=model,representation=representation,feature=feature_names[ordering],importance=selected_values,scaled_importance=if(maximum>0)selected_values/maximum else 0,rank=seq_along(ordering))
}

fit_predict_model<-function(model,x_train,y_train,x_test,class_levels,representation,seed,n_threads,settings){
  set.seed(seed)
  y_train<-factor(y_train,levels=class_levels)
  safe_feature_names<-make.names(colnames(x_train),unique=TRUE)
  original_feature_names<-colnames(x_train)
  colnames(x_train)<-safe_feature_names
  colnames(x_test)<-safe_feature_names
  importance_result<-data.frame()
  native_prediction<-NULL
  if(model=="Random Forest"){
    fit<-randomForest::randomForest(x=x_train,y=y_train,ntree=settings$rf_ntree,mtry=max(1,floor(sqrt(ncol(x_train)))),importance=TRUE)
    prob<-predict(fit,x_test,type="prob")
    importance_values<-randomForest::importance(fit,type=2)
    importance_result<-extract_top_features(importance_values[,1],original_feature_names,model,representation,settings$save_top_features)
  }else if(model=="SVM"){
    fit<-e1071::svm(x=x_train,y=y_train,kernel="radial",cost=settings$svm_cost,gamma=1/ncol(x_train),scale=TRUE,probability=TRUE)
    svm_prediction<-predict(fit,x_test,probability=TRUE)
    native_prediction<-factor(as.character(svm_prediction),levels=class_levels)
    prob<-attr(svm_prediction,"probabilities")
  }else if(model=="XGBoost"){
    labels<-as.integer(y_train)-1L
    dtrain<-xgboost::xgb.DMatrix(data=x_train,label=labels)
    dtest<-xgboost::xgb.DMatrix(data=x_test)
    parameters<-list(objective="multi:softprob",eval_metric="mlogloss",num_class=length(class_levels),max_depth=settings$xgb_max_depth,eta=settings$xgb_eta,subsample=settings$xgb_subsample,colsample_bytree=settings$xgb_colsample_bytree,nthread=n_threads,seed=seed)
    fit<-xgboost::xgb.train(params=parameters,data=dtrain,nrounds=settings$xgb_nrounds,verbose=0)
    prob<-matrix(predict(fit,dtest),ncol=length(class_levels),byrow=TRUE,dimnames=list(NULL,class_levels))
    importance_table<-xgboost::xgb.importance(feature_names=safe_feature_names,model=fit)
    if(nrow(importance_table)>0){
      original_names<-original_feature_names[match(importance_table$Feature,safe_feature_names)]
      importance_result<-extract_top_features(importance_table$Gain,original_names,model,representation,settings$save_top_features)
    }
  }else if(model=="Naive Bayes"){
    fit<-e1071::naiveBayes(x=x_train,y=y_train,laplace=0)
    prob<-predict(fit,x_test,type="raw")
  }else if(model=="KNN"){
    fit<-caret::knn3(x=x_train,y=y_train,k=settings$knn_k)
    prob<-predict(fit,x_test,type="prob")
  }else{
    stop(paste("Unknown model:",model))
  }
  prob<-normalize_probabilities(prob,class_levels)
  if(is.null(native_prediction)){
    predicted<-factor(class_levels[max.col(prob,ties.method="first")],levels=class_levels)
  }else{
    predicted<-native_prediction
  }
  list(predicted=predicted,probabilities=prob,importance=importance_result)
}

calculate_metrics<-function(truth,predicted,prob,class_levels){
  truth<-factor(truth,levels=class_levels)
  predicted<-factor(predicted,levels=class_levels)
  confusion<-table(truth,predicted)
  true_positive<-diag(confusion)
  support<-rowSums(confusion)
  predicted_count<-colSums(confusion)
  precision<-safe_divide(true_positive,predicted_count)
  recall<-safe_divide(true_positive,support)
  f1<-safe_divide(2*precision*recall,precision+recall)
  accuracy<-mean(truth==predicted)
  truth_column<-match(as.character(truth),class_levels)
  observed_probability<-prob[cbind(seq_along(truth_column),truth_column)]
  log_loss<--mean(log(pmax(observed_probability,1e-15)))
  one_hot<-matrix(0,nrow=length(truth),ncol=length(class_levels))
  one_hot[cbind(seq_along(truth_column),truth_column)]<-1
  brier_score<-mean(rowSums((prob-one_hot)^2))
  class_auc<-vapply(seq_along(class_levels),function(i){
    response<-as.integer(truth==class_levels[i])
    tryCatch(as.numeric(pROC::auc(pROC::roc(response,prob[,i],levels=c(0,1),direction="<",quiet=TRUE))),error=function(e)NA_real_)
  },numeric(1))
  overall<-data.frame(accuracy=accuracy,balanced_accuracy=mean(recall,na.rm=TRUE),macro_precision=mean(precision,na.rm=TRUE),macro_recall=mean(recall,na.rm=TRUE),macro_f1=mean(f1,na.rm=TRUE),macro_auc=mean(class_auc,na.rm=TRUE),log_loss=log_loss,brier_score=brier_score)
  per_class<-data.frame(class=class_levels,support=as.numeric(support),precision=as.numeric(precision),recall=as.numeric(recall),f1=as.numeric(f1),auc=class_auc)
  list(overall=overall,per_class=per_class,confusion=confusion)
}

summarize_with_ci<-function(data,group_columns,value_column="value"){
  data%>%dplyr::group_by(dplyr::across(dplyr::all_of(group_columns)))%>%dplyr::summarise(n=dplyr::n(),mean=mean(.data[[value_column]],na.rm=TRUE),sd=sd(.data[[value_column]],na.rm=TRUE),se=sd/sqrt(n),ci_lower=mean-qt(0.975,df=pmax(n-1,1))*se,ci_upper=mean+qt(0.975,df=pmax(n-1,1))*se,.groups="drop")
}

paired_model_comparisons<-function(metric_long){
  rows<-list()
  row_id<-1
  for(representation_name in unique(metric_long$representation)){
    for(metric_name in c("accuracy","balanced_accuracy","macro_f1")){
      subset_data<-metric_long%>%dplyr::filter(representation==representation_name,metric==metric_name)
      pairs<-combn(sort(unique(subset_data$model)),2,simplify=FALSE)
      for(pair in pairs){
        wide<-subset_data%>%dplyr::filter(model%in%pair)%>%dplyr::select(repeat_id,model,value)%>%tidyr::pivot_wider(names_from=model,values_from=value)
        difference<-wide[[pair[1]]]-wide[[pair[2]]]
        difference<-difference[is.finite(difference)]
        if(length(difference)>1){
          test<-wilcox.test(difference,mu=0,exact=FALSE)
          se<-sd(difference)/sqrt(length(difference))
          margin<-qt(0.975,df=length(difference)-1)*se
          rows[[row_id]]<-data.frame(representation=representation_name,metric=metric_name,model_1=pair[1],model_2=pair[2],n=length(difference),mean_difference=mean(difference),ci_lower=mean(difference)-margin,ci_upper=mean(difference)+margin,p_value=test$p.value)
          row_id<-row_id+1
        }
      }
    }
  }
  result<-dplyr::bind_rows(rows)
  if(nrow(result)>0)result<-result%>%dplyr::group_by(representation,metric)%>%dplyr::mutate(p_holm=p.adjust(p_value,method="holm"))%>%dplyr::ungroup()
  result
}

paired_representation_comparisons<-function(metric_long){
  rows<-list()
  row_id<-1
  for(model_name in unique(metric_long$model)){
    for(metric_name in c("accuracy","balanced_accuracy","macro_f1")){
      wide<-metric_long%>%dplyr::filter(model==model_name,metric==metric_name)%>%dplyr::select(repeat_id,representation,value)%>%tidyr::pivot_wider(names_from=representation,values_from=value)
      if(all(c("Genes","100 PCs")%in%colnames(wide))){
        difference<-wide[["100 PCs"]]-wide[["Genes"]]
        difference<-difference[is.finite(difference)]
        if(length(difference)>1){
          test<-wilcox.test(difference,mu=0,exact=FALSE)
          se<-sd(difference)/sqrt(length(difference))
          margin<-qt(0.975,df=length(difference)-1)*se
          rows[[row_id]]<-data.frame(model=model_name,metric=metric_name,n=length(difference),mean_difference_pca_minus_genes=mean(difference),ci_lower=mean(difference)-margin,ci_upper=mean(difference)+margin,p_value=test$p.value)
          row_id<-row_id+1
        }
      }
    }
  }
  result<-dplyr::bind_rows(rows)
  if(nrow(result)>0)result<-result%>%dplyr::group_by(metric)%>%dplyr::mutate(p_holm=p.adjust(p_value,method="holm"))%>%dplyr::ungroup()
  result
}

save_plot<-function(plot_object,file_stem,width,height,figure_dir){
  pdf_path<-file.path(figure_dir,paste0(file_stem,".pdf"))
  png_path<-file.path(figure_dir,paste0(file_stem,".png"))
  pdf_device<-if(capabilities("cairo"))grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(pdf_path,plot_object,width=width,height=height,device=pdf_device)
  ggplot2::ggsave(png_path,plot_object,width=width,height=height,dpi=300,bg="white")
}
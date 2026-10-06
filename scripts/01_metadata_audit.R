# ==============================================================================
# ERAMIA-RS 2026 - ZEISEL METADATA AND GROUPED-VALIDATION AUDIT
# ==============================================================================
# This script does not fit classification models. It documents the available
# cell metadata and evaluates which variables may support blocked resampling.
# Run this script before defining a group-aware cross-validation experiment.
# ==============================================================================

master_seed<-335705

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
results_dir<-file.path(project_dir,"results","metadata_audit")
table_dir<-file.path(results_dir,"tables")
figure_dir<-file.path(results_dir,"figures")
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)

required_packages<-c("scRNAseq","SummarizedExperiment","Matrix","ggplot2","dplyr","tidyr","scales")
missing_packages<-required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)>0)stop(paste("Missing packages:",paste(missing_packages,collapse=", ")))
`%>%`<-dplyr::`%>%`
options(stringsAsFactors=FALSE)
set.seed(master_seed)

collapse_examples<-function(x,n=5){
  values<-unique(as.character(x[!is.na(x)]))
  if(length(values)==0)return(NA_character_)
  paste(utils::head(values,n),collapse=" | ")
}

cramers_v<-function(x,y){
  complete<-stats::complete.cases(x,y)
  contingency<-table(x[complete],y[complete])
  if(nrow(contingency)<2||ncol(contingency)<2)return(NA_real_)
  test<-suppressWarnings(stats::chisq.test(contingency,correct=FALSE))
  denominator<-sum(contingency)*min(nrow(contingency)-1,ncol(contingency)-1)
  if(!is.finite(denominator)||denominator<=0)return(NA_real_)
  sqrt(as.numeric(test$statistic)/denominator)
}

safe_filename<-function(x){
  output<-gsub("[^A-Za-z0-9]+","_",x)
  output<-gsub("^_+|_+$","",output)
  if(!nzchar(output))output<-"metadata"
  output
}

cat("\n==============================================================================\n")
cat("ZEISEL METADATA AUDIT\n")
cat("==============================================================================\n")
cat("Loading data...\n")

sce<-scRNAseq::ZeiselBrainData()
counts_raw<-SummarizedExperiment::assay(sce,"counts")
metadata_df<-as.data.frame(SummarizedExperiment::colData(sce),stringsAsFactors=FALSE)
cell_labels<-factor(sce$level1class)
valid_cells<-!is.na(cell_labels)&Matrix::colSums(counts_raw)>0
metadata_df<-metadata_df[valid_cells,,drop=FALSE]
cell_labels<-droplevels(cell_labels[valid_cells])
cell_ids<-colnames(sce)[valid_cells]
if(is.null(cell_ids))cell_ids<-sprintf("cell_%04d",seq_len(nrow(metadata_df)))

if(!"level1class"%in%names(metadata_df))metadata_df$level1class<-as.character(cell_labels)
metadata_df$c1_run<-sub("_.*$","",cell_ids)
metadata_export<-data.frame(cell_id=cell_ids,metadata_df,check.names=FALSE)
write.csv(metadata_export,file.path(table_dir,"cell_metadata.csv"),row.names=FALSE,na="")

column_summary<-dplyr::bind_rows(lapply(names(metadata_df),function(column_name){
  values<-metadata_df[[column_name]]
  data.frame(column=column_name,class=paste(class(values),collapse="/"),n_nonmissing=sum(!is.na(values)),n_missing=sum(is.na(values)),n_unique=dplyr::n_distinct(values,na.rm=TRUE),examples=collapse_examples(values))
}))
write.csv(column_summary,file.path(table_dir,"metadata_columns.csv"),row.names=FALSE,na="")

excluded_pattern<-"level[0-9]*class|cell.?type|cluster|annotation|label|color|colour"
priority_pattern<-"group|batch|plate|mouse|animal|donor|sample|run|well|chip|individual|subject|replicate"
candidate_names<-column_summary%>%dplyr::filter(n_unique>=2,n_unique<=min(250,floor(nrow(metadata_df)/3)),!grepl(excluded_pattern,column,ignore.case=TRUE))%>%dplyr::pull(column)

candidate_rows<-list()
group_class_rows<-list()
candidate_id<-1
for(column_name in candidate_names){
  group_values<-metadata_df[[column_name]]
  group_values<-as.character(group_values)
  group_values[is.na(group_values)|!nzchar(group_values)]<-"[missing]"
  group_sizes<-table(group_values)
  contingency<-table(group_values,cell_labels)
  classes_per_group<-rowSums(contingency>0)
  groups_per_class<-colSums(contingency>0)
  missing_training_classes<-vapply(seq_len(nrow(contingency)),function(i){
    training_counts<-colSums(contingency[-i,,drop=FALSE])
    sum(contingency[i,]>0&training_counts==0)
  },numeric(1))
  candidate_rows[[candidate_id]]<-data.frame(column=column_name,priority_name_match=grepl(priority_pattern,column_name,ignore.case=TRUE),n_groups=length(group_sizes),min_group_size=min(group_sizes),median_group_size=median(as.numeric(group_sizes)),max_group_size=max(group_sizes),min_classes_per_group=min(classes_per_group),median_classes_per_group=median(classes_per_group),max_classes_per_group=max(classes_per_group),min_groups_per_class=min(groups_per_class),median_groups_per_class=median(groups_per_class),max_groups_per_class=max(groups_per_class),cramers_v_with_level1=cramers_v(group_values,cell_labels),groups_with_unseen_test_classes=sum(missing_training_classes>0),max_unseen_test_classes=max(missing_training_classes),stringsAsFactors=FALSE)
  group_class_data<-as.data.frame(contingency,stringsAsFactors=FALSE)
  colnames(group_class_data)<-c("group","class","n")
  group_class_data$column<-column_name
  group_class_rows[[candidate_id]]<-group_class_data
  candidate_id<-candidate_id+1
}

candidate_summary<-dplyr::bind_rows(candidate_rows)
group_class_counts<-dplyr::bind_rows(group_class_rows)
if(nrow(candidate_summary)>0){
  candidate_summary<-candidate_summary%>%dplyr::arrange(dplyr::desc(priority_name_match),groups_with_unseen_test_classes,cramers_v_with_level1,dplyr::desc(median_group_size))
  group_class_counts<-group_class_counts%>%dplyr::mutate(proportion_within_group=n/ave(n,column,group,FUN=sum))
}else{
  candidate_summary<-data.frame(column=character(),priority_name_match=logical(),n_groups=integer(),min_group_size=integer(),median_group_size=numeric(),max_group_size=integer(),min_classes_per_group=integer(),median_classes_per_group=numeric(),max_classes_per_group=integer(),min_groups_per_class=integer(),median_groups_per_class=numeric(),max_groups_per_class=integer(),cramers_v_with_level1=numeric(),groups_with_unseen_test_classes=integer(),max_unseen_test_classes=integer())
  group_class_counts<-data.frame(group=character(),class=character(),n=integer(),column=character(),proportion_within_group=numeric())
}
write.csv(candidate_summary,file.path(table_dir,"group_candidate_summary.csv"),row.names=FALSE,na="")
write.csv(group_class_counts,file.path(table_dir,"group_class_counts.csv"),row.names=FALSE,na="")

priority_candidates<-candidate_summary%>%dplyr::filter(priority_name_match)%>%dplyr::slice_head(n=8)%>%dplyr::pull(column)
if(length(priority_candidates)==0&&nrow(candidate_summary)>0)priority_candidates<-candidate_summary%>%dplyr::slice_head(n=4)%>%dplyr::pull(column)

if(length(priority_candidates)>0){
  plot_data<-group_class_counts%>%dplyr::filter(column%in%priority_candidates)%>%dplyr::group_by(column,group)%>%dplyr::mutate(group_total=sum(n))%>%dplyr::ungroup()%>%dplyr::filter(group_total>0)%>%dplyr::mutate(group=factor(group,levels=unique(group[order(column,group_total)])))
  p_composition<-ggplot2::ggplot(plot_data,ggplot2::aes(x=group,y=n,fill=class))+ggplot2::geom_col(position="fill",width=0.85)+ggplot2::facet_wrap(~column,scales="free_x",ncol=1)+ggplot2::scale_y_continuous(labels=scales::percent_format())+ggplot2::labs(x="Metadata group",y="Cell-type composition",fill="Level-one class")+ggplot2::theme_minimal(base_size=9)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=90,vjust=0.5,hjust=1,size=6),legend.position="bottom",strip.text=ggplot2::element_text(face="bold"))
  ggplot2::ggsave(file.path(figure_dir,"candidate_group_composition.pdf"),p_composition,width=9,height=max(4,2.2*length(priority_candidates)),device=if(capabilities("cairo"))grDevices::cairo_pdf else "pdf",limitsize=FALSE)
  ggplot2::ggsave(file.path(figure_dir,"candidate_group_composition.png"),p_composition,width=9,height=max(4,2.2*length(priority_candidates)),dpi=250,bg="white",limitsize=FALSE)
}

capture.output(sessionInfo(),file=file.path(results_dir,"sessionInfo.txt"))

cat("Cells audited:",nrow(metadata_df),"\n")
cat("Metadata columns:",ncol(metadata_df),"\n")
cat("Candidate grouping variables:",nrow(candidate_summary),"\n\n")
if(nrow(candidate_summary)>0)print(candidate_summary%>%dplyr::slice_head(n=15))
cat("\nAudit completed. Upload these files before blocked resampling is defined:\n")
cat("  ",file.path(table_dir,"metadata_columns.csv"),"\n")
cat("  ",file.path(table_dir,"group_candidate_summary.csv"),"\n")
cat("  ",file.path(table_dir,"group_class_counts.csv"),"\n")
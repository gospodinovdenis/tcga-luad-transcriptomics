
options(timeout=600)
suppressPackageStartupMessages({library(GEOquery); library(limma)})
dir.create("results",showWarnings=FALSE)
gset <- getGEO("GSE32863",GSEMatrix=TRUE,AnnotGPL=FALSE)
es <- gset[[1]]
expr <- exprs(es)
pd <- pData(es)
write.csv(pd,"results/geo_sample_metadata_raw.csv",row.names=TRUE)
write.csv(fData(es),"results/geo_feature_metadata_raw.csv",row.names=TRUE)
cat("Expression:",dim(expr),"\n")
cat("Metadata columns:",paste(names(pd),collapse=" | "),"\n")
# Use GEO sample titles: 'Adjacent non-tumor lung tissue 05L12_N (expression)'
titles <- as.character(pd$title)
type <- ifelse(grepl("adjacent non-tumor",titles,ignore.case=TRUE),"Normal",
               ifelse(grepl("lung adenocarcinoma",titles,ignore.case=TRUE),"Tumor",NA))
patient <- sub(".*?([[:alnum:]]+_[NT]).*","\\1",titles)
patient <- sub("_[NT]$","",patient)
if(anyNA(type) || any(!grepl("^[[:alnum:]]+$",patient))) stop("Unexpected sample titles; inspect results/geo_sample_metadata_raw.csv")
meta <- data.frame(gsm=colnames(expr),title=titles,patient=patient,type=type)
# Retain matched pairs only; paired design controls inter-patient variation.
tab <- table(meta$patient,meta$type)
keep_patients <- rownames(tab)[tab[,"Tumor"]==1 & tab[,"Normal"]==1]
meta <- meta[meta$patient %in% keep_patients,,drop=FALSE]
expr <- expr[,meta$gsm,drop=FALSE]
if(length(keep_patients)<20) stop("Too few matched pairs; inspect GEO metadata")
if(!identical(colnames(expr),meta$gsm)) stop("Sample order mismatch")
# GEO reports robust-spline-normalized log2 intensities; no DESeq2 on microarrays.
design <- model.matrix(~ factor(patient) + factor(type,levels=c("Normal","Tumor")),meta)
fit <- eBayes(lmFit(expr,design),trend=TRUE)
coef_name <- grep("^factor\\(type.*Tumor$",colnames(design),value=TRUE)
if(length(coef_name)!=1) stop("Could not locate tumor coefficient")
tt <- topTable(fit,coef=coef_name,number=Inf,sort.by="none")
tt$probe_id <- rownames(tt)
# Platform mapping: use GPL6884 annotation; fail loudly if gene symbol unavailable.
gpl <- getGEO("GPL6884",AnnotGPL=FALSE)
ann <- Table(gpl)
write.csv(head(ann,50),"results/platform_annotation_preview.csv",row.names=FALSE)
symbol_candidates <- grep("^(Symbol|Gene Symbol|Gene_Symbol|SYMBOL|gene_symbol)$",
                          names(ann),ignore.case=TRUE,value=TRUE)
if(length(symbol_candidates)==0) symbol_candidates <- grep("symbol",names(ann),ignore.case=TRUE,value=TRUE)
if(length(symbol_candidates)==0) stop("No gene symbol field found in GPL6884 annotation; inspect platform_annotation_preview.csv")
id_candidates <- grep("^(ID|ID_REF|Probe_Id|ProbeID)$",names(ann),ignore.case=TRUE,value=TRUE)
if(length(id_candidates)==0) stop("No probe ID column in GPL6884")
mapping <- data.frame(probe_id=as.character(ann[[id_candidates[1]]]),
                      gene_symbol=as.character(ann[[symbol_candidates[1]]]))
mapping$gene_symbol <- trimws(mapping$gene_symbol)
mapping <- mapping[!is.na(mapping$gene_symbol) & grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$",mapping$gene_symbol),]
mapping <- mapping[!duplicated(mapping$probe_id),]
tt <- merge(tt,mapping,by="probe_id")
# Multiple probes per gene: select probe with highest across-sample variance,
# independent of the tumor-normal contrast.
probe_var <- apply(expr,1,var,na.rm=TRUE)
tt$probe_variance <- probe_var[tt$probe_id]
tt <- tt[order(tt$gene_symbol,-tt$probe_variance,tt$probe_id),]
tt <- tt[!duplicated(tt$gene_symbol),]
tt <- tt[order(tt$adj.P.Val),]
write.csv(tt,"results/external_gene_differential_expression.csv",row.names=FALSE)
write.csv(meta,"results/external_sample_metadata.csv",row.names=FALSE)
selected <- tt$probe_id
external <- expr[selected,,drop=FALSE]
rownames(external) <- tt$gene_symbol[match(selected,tt$probe_id)]
write.csv(external,"results/external_log2_expression.csv",row.names=TRUE)
cat("Paired patients:",length(keep_patients),"Mapped genes:",nrow(tt),"\n")
writeLines(capture.output(sessionInfo()),"results/R_session_info.txt")

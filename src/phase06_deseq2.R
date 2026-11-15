suppressPackageStartupMessages({library(DESeq2);library(ggplot2);library(pheatmap)})
a<-commandArgs(trailingOnly=TRUE); if(length(a)<2)stop('Usage: Rscript phase6_deseq2.R <phase6_output> <deseq2_output>')
indir<-a[1]; out<-a[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
ct<-read.csv(gzfile(file.path(indir,'processed/raw_counts.csv.gz')),row.names=1,check.names=FALSE)
meta<-read.csv(file.path(indir,'processed/aligned_metadata.csv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(colnames(ct)),!anyDuplicated(meta$sample_barcode),setequal(colnames(ct),meta$sample_barcode))
meta<-meta[match(colnames(ct),meta$sample_barcode),,drop=FALSE];rownames(meta)<-meta$sample_barcode
stopifnot(identical(colnames(ct),rownames(meta)),all(is.finite(as.matrix(ct))),all(as.matrix(ct)>=0),all(as.matrix(ct)==round(as.matrix(ct))))
meta$condition<-factor(meta$sample_type,levels=c('Solid Tissue Normal','Primary Tumor'))
meta$tss<-factor(meta$tss)
keep<-rowSums(as.matrix(ct)>=10)>=max(3,min(table(meta$condition)))
ct<-ct[keep,,drop=FALSE]
dds<-DESeqDataSetFromMatrix(round(as.matrix(ct)),meta,design=~1)
dds<-estimateSizeFactors(dds)
write.csv(data.frame(gene_id=rownames(dds),counts(dds,normalized=TRUE),check.names=FALSE),gzfile(file.path(out,'normalized_counts.csv.gz')),row.names=FALSE)
write.csv(data.frame(sample_barcode=colnames(dds),size_factor=as.numeric(sizeFactors(dds))),file.path(out,'size_factors.csv'),row.names=FALSE)
vsd<-vst(dds,blind=TRUE)
write.csv(data.frame(gene_id=rownames(vsd),assay(vsd),check.names=FALSE),gzfile(file.path(out,'vst_matrix.csv.gz')),row.names=FALSE)
pdf(file.path(out,'pca_vst.pdf'),width=9,height=7);print(plotPCA(vsd,intgroup=c('condition','tss')));dev.off()
pdf(file.path(out,'sample_distances.pdf'),width=11,height=10);pheatmap(as.matrix(dist(t(assay(vsd)))),main='VST sample distances');dev.off()
ctab<-table(meta$tss,meta$condition)
strong<-sum(rowSums(ctab>=2)==2)
# Conservative safety gate. Additional technical covariates still need manual review.
if(min(table(meta$condition))>=3 && strong>=2 && !anyDuplicated(meta$case_barcode)){
  design(dds)<-~tss+condition
  if(qr(model.matrix(design(dds),as.data.frame(colData(dds))))$rank==ncol(model.matrix(design(dds),as.data.frame(colData(dds))))){
    dds<-DESeq(dds)
    res<-results(dds,contrast=c('condition','Primary Tumor','Solid Tissue Normal'),alpha=.05)
    write.csv(data.frame(gene_id=rownames(res),as.data.frame(res)),file.path(out,'exploratory_dge_tss_adjusted.csv'),row.names=FALSE)
    writeLines('EXPLORATORY differential expression only; TSS adjusted. Investigate other covariates, tumor purity and site effects before publication.',file.path(out,'DGE_STATUS.txt'))
  }else writeLines('DGE BLOCKED: model matrix not full rank',file.path(out,'DGE_STATUS.txt'))
}else writeLines('DGE BLOCKED: insufficient overlap or replication across tissue-source sites. Normalization/VST completed.',file.path(out,'DGE_STATUS.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))

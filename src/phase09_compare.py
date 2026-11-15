
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from scipy.stats import spearmanr, binomtest, wilcoxon
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import json
out=Path("results");out.mkdir(exist_ok=True)
tcga=pd.read_csv("input/annotated_dge.csv")
ext=pd.read_csv(out/"external_gene_differential_expression.csv")
tcga=tcga.rename(columns={"gene_name":"gene_symbol","log2FoldChange":"tcga_log2FC","padj":"tcga_FDR"})
ext=ext.rename(columns={"logFC":"external_log2FC","adj.P.Val":"external_FDR"})
tcga=tcga.dropna(subset=["gene_symbol"]).sort_values("tcga_FDR").drop_duplicates("gene_symbol")
both=tcga.merge(ext[["gene_symbol","external_log2FC","external_FDR"]],on="gene_symbol")
both=both.replace([np.inf,-np.inf],np.nan).dropna(subset=["tcga_log2FC","external_log2FC"])
both["direction_agreement"]=np.sign(both.tcga_log2FC)==np.sign(both.external_log2FC)
both["external_replication"]=(both.external_FDR<0.05)&both.direction_agreement
both.to_csv(out/"cross_cohort_gene_comparison.csv",index=False)
candidates=["COL1A1","COL3A1","MMP1","TOP2A","MKI67","AGER","SFTPC","SFTPA1"]
both[both.gene_symbol.isin(candidates)].to_csv(out/"candidate_gene_validation.csv",index=False)
# Selected before external analysis using original TCGA criteria
sig=both[(both.tcga_FDR<0.05)&(both.tcga_log2FC.abs()>1)]
n=int(len(sig));k=int(sig.direction_agreement.sum())
rho,p=spearmanr(both.tcga_log2FC,both.external_log2FC)
summary={"external_dataset":"GSE32863","matched_genes":len(both),
         "tcga_significant_mapped":n,"direction_concordance":k/n if n else None,
         "direction_concordance_n":k,"external_replication_fraction":float(sig.external_replication.mean()) if n else None,
         "spearman_rho_all_matched":float(rho),"spearman_p":float(p),
         "direction_binomial_p":float(binomtest(k,n,0.5,alternative="greater").pvalue) if n else None}
# Leading-edge sets from original TCGA GSEA, held fixed for external test.
# Per-sample score = average within-gene standardized expression, not GSVA/GSEA NES.
expr=pd.read_csv(out/"external_log2_expression.csv",index_col=0)
meta=pd.read_csv(out/"external_sample_metadata.csv")
assert list(expr.columns)==list(meta.gsm)
z=expr.sub(expr.mean(axis=1),axis=0).div(expr.std(axis=1).replace(0,np.nan),axis=0)
pathways=[]
for lib in ["Reactome","GO_BP"]:
    gsea=pd.read_csv(f"input/gsea_{lib}.csv").sort_values("FDR q-val")
    gsea=gsea[gsea["FDR q-val"]<0.05]
    # Select top 5 positive and 5 negative per library, preselected from TCGA
    chosen=pd.concat([gsea[gsea.NES>0].nlargest(5,"NES"),gsea[gsea.NES<0].nsmallest(5,"NES")])
    for _,row in chosen.iterrows():
        genes=set(str(row.Lead_genes).split(";"))
        matched=sorted(genes.intersection(z.index))
        if len(matched)<5: continue
        score=z.loc[matched].mean(axis=0)
        vals=pd.DataFrame({"gsm":score.index,"score":score.values}).merge(meta,on="gsm")
        paired=vals.pivot(index="patient",columns="type",values="score").dropna()
        diff=paired["Tumor"]-paired["Normal"]
        if len(diff)<10: continue
        pval=wilcoxon(diff).pvalue if (diff!=0).any() else 1.0
        pathways.append({"library":lib,"pathway":row.Term,"tcga_NES":row.NES,
                         "tcga_FDR":row["FDR q-val"],"external_genes_scored":len(matched),
                         "paired_patients":len(diff),"external_mean_paired_score_difference":float(diff.mean()),
                         "external_p":float(pval),"direction_agreement":bool(np.sign(diff.mean())==np.sign(row.NES))})
path=pd.DataFrame(pathways)
if len(path):
    path["external_FDR"]=multipletests(path.external_p,method="fdr_bh")[1]
    path["external_replication"]=(path.external_FDR<0.05)&path.direction_agreement
path.to_csv(out/"pathway_signature_validation.csv",index=False)
summary["pathways_evaluated"]=len(path)
summary["pathways_replicated"]=int(path.external_replication.sum()) if len(path) else 0
(out/"validation_summary.json").write_text(json.dumps(summary,indent=2))
fig,ax=plt.subplots(figsize=(7,6))
ax.scatter(both.tcga_log2FC,both.external_log2FC,s=5,alpha=.2)
ax.axhline(0,color="grey",linewidth=.6);ax.axvline(0,color="grey",linewidth=.6)
ax.set(xlabel="TCGA log2 fold change (RNA-seq)",ylabel="GSE32863 log2 fold change (microarray)",
       title=f"Cross-cohort gene effect concordance (Spearman rho={rho:.2f})")
fig.tight_layout();fig.savefig(out/"cross_cohort_concordance.png",dpi=220);plt.close(fig)
if len(path):
    fig,ax=plt.subplots(figsize=(9,max(4,len(path)*.33)))
    ordered=path.sort_values("external_mean_paired_score_difference")
    ax.barh(range(len(ordered)),ordered.external_mean_paired_score_difference)
    ax.set_yticks(range(len(ordered)),[str(x)[:52] for x in ordered.pathway])
    ax.set_xlabel("Mean paired tumor − normal standardized signature score")
    fig.tight_layout();fig.savefig(out/"pathway_signature_validation.png",dpi=220);plt.close(fig)
print(json.dumps(summary,indent=2))

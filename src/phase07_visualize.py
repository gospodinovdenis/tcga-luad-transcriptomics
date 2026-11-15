import argparse, zipfile, json, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

def readzip(path, suffix):
    with zipfile.ZipFile(path) as z:
        matches=[n for n in z.namelist() if n.endswith(suffix)]
        if len(matches)!=1: raise ValueError(f'Expected one {suffix}, found {matches}')
        with z.open(matches[0]) as f: return pd.read_csv(f)

def main(args):
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    d=readzip(args.deseq,'exploratory_dge_tss_adjusted.csv')
    a=readzip(args.results,'gene_annotations.csv')
    m=readzip(args.results,'aligned_metadata.csv')
    with zipfile.ZipFile(args.deseq) as z:
        vname=next(n for n in z.namelist() if n.endswith('/vst_matrix.csv.gz'))
        with gzip.GzipFile(fileobj=z.open(vname)) as f: vst=pd.read_csv(f,index_col=0)
    assert d.gene_id.is_unique and a.gene_id.is_unique
    assert set(m.sample_barcode)==set(vst.columns)
    assert not vst.isna().any().any()
    d=d.merge(a,on='gene_id',how='left',validate='one_to_one')
    d['significance']=np.select([(d.padj<.05)&(d.log2FoldChange>1),(d.padj<.05)&(d.log2FoldChange< -1)],['Higher in tumor','Lower in tumor'],default='Other')
    d.sort_values(['padj','pvalue'],na_position='last').to_csv(out/'annotated_dge.csv',index=False)
    top=d[(d.padj<.05)&(d.log2FoldChange.abs()>1)].sort_values('padj').head(50)
    top.to_csv(out/'top50_dge.csv',index=False)
    fig,ax=plt.subplots(figsize=(9,7))
    y=-np.log10(d.padj.clip(lower=1e-300))
    for cat,color in [('Other','#9ca3af'),('Higher in tumor','#c05621'),('Lower in tumor','#2563eb')]:
        ix=d.significance.eq(cat)&y.notna();ax.scatter(d.loc[ix,'log2FoldChange'],y[ix],s=5,alpha=.5,label=f'{cat} ({ix.sum()})',c=color,rasterized=True)
    ax.axvline(-1,ls='--',c='gray');ax.axvline(1,ls='--',c='gray');ax.axhline(-np.log10(.05),ls='--',c='gray')
    ax.set(xlabel='log2 fold change (tumor vs normal)',ylabel='-log10 adjusted p-value',title='TCGA-LUAD exploratory site-adjusted DESeq2');ax.legend(markerscale=2,fontsize=8);fig.tight_layout();fig.savefig(out/'volcano.png',dpi=180);plt.close(fig)
    genes=[g for g in top.gene_id if g in vst.index]
    X=vst.loc[genes,m.sample_barcode].copy();X=X.sub(X.mean(axis=1),axis=0).div(X.std(axis=1).replace(0,1),axis=0)
    order=leaves_list(linkage(pdist(X.T),method='average'))
    X=X.iloc[:,order]
    fig,ax=plt.subplots(figsize=(11,12));im=ax.imshow(X,aspect='auto',cmap='coolwarm',vmin=-2.5,vmax=2.5)
    labels=d.set_index('gene_id').reindex(X.index).gene_name.fillna(pd.Series(X.index,index=X.index))
    ax.set_yticks(range(len(X)),labels,fontsize=6);ax.set_xticks(range(X.shape[1]),[x[5:12] for x in X.columns],rotation=90,fontsize=7)
    fig.colorbar(im,ax=ax,label='Gene-wise VST z-score',shrink=.5);fig.tight_layout();fig.savefig(out/'top50_heatmap.png',dpi=180);plt.close(fig)
    summary={'samples':len(m),'genes_tested':len(d),'fdr_0_05':int((d.padj<.05).sum()),'higher_tumor_fdr_and_lfc':int((d.significance=='Higher in tumor').sum()),'lower_tumor_fdr_and_lfc':int((d.significance=='Lower in tumor').sum()),'missing_padj':int(d.padj.isna().sum()),'tss_by_type':pd.crosstab(m.tss,m.sample_type).to_dict(),'limitations':['Exploratory association; site-adjusted model, not causal','No independent cohort validation','Pathway enrichment and leave-one-site-out DESeq2 reruns not performed by this script']}
    (out/'phase7_summary.json').write_text(json.dumps(summary,indent=2))
    (out/'pathway_input_ensembl.txt').write_text('\n'.join(g.split('.')[0] for g in d.loc[d.padj.lt(.05)&d.log2FoldChange.abs().gt(1),'gene_id'])+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',required=True);p.add_argument('--deseq',required=True);p.add_argument('--out',default='phase7_output');main(p.parse_args())

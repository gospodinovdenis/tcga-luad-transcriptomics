"""TCGA-LUAD Phase 7B: reproducible GO/Reactome ORA and preranked GSEA."""
from pathlib import Path
import argparse, json, platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def prepare(dge, annotations, out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(dge); a=pd.read_csv(annotations)
    assert {'gene_id','stat','padj','log2FoldChange'}.issubset(d.columns)
    assert {'gene_id','gene_name'}.issubset(a.columns)
    assert d.gene_id.is_unique and a.gene_id.is_unique
    m=d.merge(a[['gene_id','gene_name','gene_type']],on='gene_id',how='left',validate='one_to_one')
    m.to_csv(out/'annotated_dge.csv',index=False)
    valid=m.loc[m.stat.notna() & np.isfinite(m.stat) & m.gene_name.notna() & m.gene_name.ne('')].copy()
    valid=valid[valid.gene_name.str.match(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')]
    # Collapse duplicate HGNC symbols deterministically by largest absolute Wald statistic.
    valid['abs_stat']=valid.stat.abs()
    valid=valid.sort_values(['abs_stat','gene_id'],ascending=[False,True]).drop_duplicates('gene_name')
    rank=valid[['gene_name','stat']].sort_values('stat',ascending=False)
    rank.to_csv(out/'gsea_ranked_genes.rnk',sep='\t',index=False,header=False)
    universe=sorted(valid.gene_name.unique())
    hits=valid[(valid.padj<.05)&(valid.log2FoldChange.abs()>1)]
    up=sorted(hits.loc[hits.log2FoldChange>0,'gene_name'].unique())
    down=sorted(hits.loc[hits.log2FoldChange<0,'gene_name'].unique())
    for name,genes in [('tested_gene_universe.txt',universe),('significant_up.txt',up),('significant_down.txt',down)]:
        (out/name).write_text('\n'.join(genes)+'\n')
    info={'n_dge_rows':len(d),'n_ranked_symbols':len(rank),'n_universe_symbols':len(universe),'n_up':len(up),'n_down':len(down),'python':platform.python_version(),'rank_metric':'DESeq2 Wald statistic','significance':'padj < 0.05 and abs(log2FoldChange) > 1','duplicate_symbols':'max absolute Wald statistic, then gene_id'}
    (out/'preflight_summary.json').write_text(json.dumps(info,indent=2))
    return info

def analyze(out,permutations=1000,seed=42):
    import gseapy as gp
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    rank=pd.read_csv(out/'gsea_ranked_genes.rnk',sep='\t',header=None,names=['gene_name','stat'])
    universe=set((out/'tested_gene_universe.txt').read_text().splitlines())
    libs={'GO_BP':'GO_Biological_Process_2023','Reactome':'Reactome_2022'}
    log={}
    for label,lib in libs.items():
        genesets=gp.get_library(name=lib,organism='Human')
        # Explicitly restrict gene sets to the DESeq2 tested universe for ORA.
        restricted={term:sorted(set(genes)&universe) for term,genes in genesets.items()}
        restricted={term:genes for term,genes in restricted.items() if 10<=len(genes)<=500}
        log[label]={'library':lib,'library_terms':len(genesets),'tested_universe_terms':len(restricted)}
        for direction in ('up','down'):
            genes=(out/f'significant_{direction}.txt').read_text().splitlines()
            if not genes: continue
            enr=gp.enrich(gene_list=genes,gene_sets=restricted,background=len(universe),outdir=None,verbose=False)
            result=enr.results.sort_values('Adjusted P-value')
            result.to_csv(out/f'ora_{label}_{direction}.csv',index=False)
            top=result.head(15).iloc[::-1]
            if len(top):
                fig,ax=plt.subplots(figsize=(10,max(5,len(top)*.43)))
                ax.barh(top.Term, -np.log10(top['Adjusted P-value'].clip(lower=1e-300)))
                ax.set_xlabel('-log10 adjusted p-value');ax.set_title(f'{label} ORA: {direction} in tumors')
                fig.tight_layout();fig.savefig(out/f'ora_{label}_{direction}.png',dpi=160,bbox_inches='tight');plt.close(fig)
        # Preranked GSEA uses all ranked tested genes, not only significant hits.
        pre=gp.prerank(rnk=rank,gene_sets=restricted,permutation_num=permutations,min_size=10,max_size=500,seed=seed,threads=2,outdir=None,verbose=False)
        gres=pre.res2d.sort_values('FDR q-val')
        gres.to_csv(out/f'gsea_{label}.csv',index=False)
        if len(gres):
            top=gres.head(20).iloc[::-1].copy()
            top['NES']=pd.to_numeric(top['NES'])
            fig,ax=plt.subplots(figsize=(10,max(5,len(top)*.38)))
            ax.barh(top.Term,top.NES)
            ax.set_xlabel('Normalized enrichment score (NES)');ax.set_title(f'{label} preranked GSEA (top FDR terms)')
            fig.tight_layout();fig.savefig(out/f'gsea_{label}_nes.png',dpi=160,bbox_inches='tight');plt.close(fig)
    (out/'enrichment_run.json').write_text(json.dumps({'libraries':log,'permutations':permutations,'seed':seed,'gseapy_version':gp.__version__},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dge',required=True);p.add_argument('--annotations',required=True);p.add_argument('--out',required=True);p.add_argument('--enrich',action='store_true');p.add_argument('--permutations',type=int,default=1000)
    a=p.parse_args();print(json.dumps(prepare(a.dge,a.annotations,a.out),indent=2))
    if a.enrich: analyze(a.out,a.permutations)

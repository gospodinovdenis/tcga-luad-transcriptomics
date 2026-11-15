import argparse, hashlib, io, json, os, platform, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import requests
BASE='https://api.gdc.cancer.gov'
TYPES=['Primary Tumor','Solid Tissue Normal']
def request_json(endpoint,payload=None):
    for i in range(5):
        try:
            r=requests.post(BASE+endpoint,json=payload,timeout=120) if payload else requests.get(BASE+endpoint,timeout=120)
            r.raise_for_status(); return r.json()
        except (requests.RequestException,ValueError):
            if i==4: raise
            time.sleep(2**i)
def discover():
    filt={'op':'and','content':[{'op':'in','content':{'field':'cases.project.project_id','value':['TCGA-LUAD']}},{'op':'in','content':{'field':'data_type','value':['Gene Expression Quantification']}},{'op':'in','content':{'field':'analysis.workflow_type','value':['STAR - Counts']}},{'op':'in','content':{'field':'access','value':['open']}}]}
    fields='file_id,file_name,md5sum,file_size,analysis.workflow_type,cases.case_id,cases.submitter_id,cases.samples.sample_id,cases.samples.submitter_id,cases.samples.sample_type'
    result=[]; start=0
    while True:
        d=request_json('/files',{'filters':filt,'fields':fields,'format':'JSON','size':500,'from':start})['data']
        result.extend(d['hits']); start+=len(d['hits'])
        if start>=d['pagination']['total'] or not d['hits']: break
    records=[]
    for f in result:
        # Files with multiple case/sample associations are not silently assigned.
        cases=f.get('cases',[])
        if len(cases)!=1: continue
        case=cases[0]; samples=case.get('samples',[])
        if len(samples)!=1: continue
        sample=samples[0]
        if sample.get('sample_type') not in TYPES: continue
        barcode=sample.get('submitter_id',''); patient=case.get('submitter_id','')
        if not barcode or not patient: continue
        records.append({'file_id':f['file_id'],'file_name':f.get('file_name'),'md5sum':f.get('md5sum'),'file_size':f.get('file_size'),'case_id':case.get('case_id'),'case_barcode':patient,'sample_id':sample.get('sample_id'),'sample_barcode':barcode,'sample_type':sample['sample_type'],'tss':patient.split('-')[1] if patient.startswith('TCGA-') else ''})
    return pd.DataFrame(records)
def select(df,max_per_group=60,seed=42):
    if df.empty: raise ValueError('No STAR count files discovered; inspect GDC query and API response')
    df=df.sort_values(['sample_barcode','file_id']).drop_duplicates('sample_barcode')
    # One sample per patient, deterministic: prefer normal if both exist to preserve scarce normals.
    df['priority']=df.sample_type.map({'Solid Tissue Normal':0,'Primary Tumor':1})
    df=df.sort_values(['priority','sample_barcode','file_id']).drop_duplicates('case_barcode').drop(columns='priority')
    ct=pd.crosstab(df.tss,df.sample_type).reindex(columns=TYPES,fill_value=0)
    overlapping=ct.index[(ct>0).all(axis=1)].tolist()
    eligible=df[df.tss.isin(overlapping)].copy()
    selected=[]
    for typ in TYPES:
        group=eligible[eligible.sample_type==typ]
        # Prefer site diversity, then fill to requested maximum; stable randomized ordering.
        group=group.sample(frac=1,random_state=seed).sort_values('tss',kind='stable')
        pools=[g for _,g in group.groupby('tss',sort=True)]
        round_robin=[]
        for idx in range(max((len(g) for g in pools),default=0)):
            round_robin.extend([g.iloc[idx] for g in pools if idx<len(g)])
        selected.extend(round_robin[:max_per_group])
    out=pd.DataFrame(selected)
    return out.reset_index(drop=True),ct,overlapping

def fetch_file(row,dest):
    path=dest/(row.file_id+'.tsv')
    if path.exists() and row.md5sum and hashlib.md5(path.read_bytes()).hexdigest()==row.md5sum: return path
    for i in range(4):
        try:
            with requests.get(BASE+'/data/'+row.file_id,stream=True,timeout=180) as r:
                r.raise_for_status()
                with open(str(path)+'.part','wb') as h:
                    for chunk in r.iter_content(1024*1024):
                        if chunk:h.write(chunk)
            os.replace(str(path)+'.part',path)
            if row.md5sum and hashlib.md5(path.read_bytes()).hexdigest()!=row.md5sum: raise ValueError('MD5 mismatch '+row.file_id)
            return path
        except (requests.RequestException,ValueError):
            if i==3: raise
            time.sleep(2**i)

def parse(path):
    df=pd.read_csv(path,sep='\t',comment='#',low_memory=False)
    required=['gene_id','gene_name','gene_type','unstranded','tpm_unstranded','fpkm_unstranded']
    missing=set(required)-set(df.columns)
    if missing: raise ValueError(f'Unexpected GDC STAR TSV schema: missing {missing}; got {df.columns.tolist()}')
    df=df[df.gene_id.astype(str).str.startswith('ENSG')].copy()
    if df.gene_id.duplicated().any(): raise ValueError('Duplicate gene IDs')
    for field in ['unstranded','tpm_unstranded','fpkm_unstranded']:
        df[field]=pd.to_numeric(df[field],errors='raise')
        if not np.isfinite(df[field]).all() or (df[field]<0).any(): raise ValueError('Invalid '+field)
    if not np.equal(df.unstranded, np.floor(df.unstranded)).all():raise ValueError('Non-integer raw counts')
    return df.set_index('gene_id')
def run(args):
    out=Path(args.output); (out/'raw').mkdir(parents=True,exist_ok=True); (out/'processed').mkdir(exist_ok=True); (out/'provenance').mkdir(exist_ok=True); (out/'qc').mkdir(exist_ok=True)
    discovered=discover(); discovered.to_csv(out/'provenance'/'discovered_files.csv',index=False)
    selected,site_table,overlap=select(discovered,args.max_per_group,args.seed)
    site_table.to_csv(out/'provenance'/'available_tss_by_tissue.csv')
    selected.to_csv(out/'provenance'/'selected_manifest.csv',index=False)
    if selected.empty: raise RuntimeError('No overlapping TSS groups: cohort expansion cannot support site-aware comparison; no DE test will be run')
    print('Selected',len(selected),'samples from overlapping sites',overlap,flush=True)
    tables=[]; annotations=None; kept=[]; failures=[]
    for idx,row in enumerate(selected.itertuples(index=False),1):
        try:
            t=parse(fetch_file(row,out/'raw'))
            if annotations is None: annotations=t[['gene_name','gene_type']].copy()
            if not t.index.equals(annotations.index): raise ValueError('Gene identifiers/order mismatch across files')
            tables.append((row.sample_barcode,t[['unstranded','tpm_unstranded','fpkm_unstranded']]))
            kept.append(row.sample_barcode)
            print(idx,'/',len(selected),row.sample_barcode,flush=True)
        except Exception as e:
            failures.append({'file_id':row.file_id,'sample_barcode':row.sample_barcode,'error':str(e)})
            print('FAILED',row.sample_barcode,str(e),flush=True)
    if not tables:raise RuntimeError('No files successfully parsed')
    meta=selected.set_index('sample_barcode').loc[kept].reset_index()
    matrices={}
    for field,name in [('unstranded','raw_counts'),('tpm_unstranded','tpm'),('fpkm_unstranded','fpkm')]:
        mat=pd.concat([t[field].rename(s) for s,t in tables],axis=1)
        if mat.isna().any().any():raise ValueError('Unexpected missing values in '+field)
        mat.to_csv(out/'processed'/(name+'.csv.gz'),compression='gzip',index_label='gene_id')
        matrices[name]=mat
    annotations.to_csv(out/'processed'/'gene_annotations.csv',index_label='gene_id')
    meta.to_csv(out/'processed'/'aligned_metadata.csv',index=False)
    raw=matrices['raw_counts']; sums=raw.sum(axis=0); detected=(raw>0).sum(axis=0)
    qc=meta.set_index('sample_barcode')[['case_barcode','sample_type','tss']].copy()
    qc['assigned_gene_counts']=sums; qc['detected_genes']=detected; qc['zero_fraction']=(raw==0).mean(axis=0)
    qc['within_group_count_ratio']=qc.assigned_gene_counts/qc.groupby('sample_type').assigned_gene_counts.transform('median')
    qc['review_low_counts']=qc.within_group_count_ratio<0.25
    qc.to_csv(out/'qc'/'sample_qc.csv')
    ct=pd.crosstab(meta.tss,meta.sample_type).reindex(columns=TYPES,fill_value=0)
    ct.to_csv(out/'qc'/'selected_tss_by_tissue.csv')
    # Conservative gate: >=3 samples/group, >=2 overlapping sites with >=2 per tissue, no duplicate patients.
    strong_sites=int(((ct>=2).all(axis=1)).sum())
    n_by=meta.sample_type.value_counts()
    ready=all(n_by.get(t,0)>=3 for t in TYPES) and strong_sites>=2 and not meta.case_barcode.duplicated().any() and not failures
    report={'status':'eligible_for_exploratory_dge' if ready else 'exploratory_normalization_only','n_samples':len(meta),'n_genes':len(raw),'group_counts':n_by.to_dict(),'overlapping_sites_selected':int(((ct>0).all(axis=1)).sum()),'sites_with_at_least_two_per_group':strong_sites,'flagged_low_count_samples':int(qc.review_low_counts.sum()),'failed_downloads':failures,'warning':'TSS balance does not establish absence of batch confounding; review other covariates and GDC QC before biological claims','seed':args.seed,'max_per_group':args.max_per_group,'python':sys.version,'pandas':pd.__version__}
    (out/'qc'/'readiness.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='phase6_output');p.add_argument('--max-per-group',type=int,default=60);p.add_argument('--seed',type=int,default=42);run(p.parse_args())

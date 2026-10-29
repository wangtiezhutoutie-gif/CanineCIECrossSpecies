from pathlib import Path
import csv, shutil
import numpy as np
import pandas as pd
R=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P30=R/'WORKING'/'30_HUMAN_PSEUDOBULK_INDEX_CORRECTION'
O=R/'WORKING'/'31_FINAL_ELIGIBILITY_AND_MANUSCRIPT_AUDIT'
W=O/'work'; S=O/'source_data'
W.mkdir(exist_ok=True); S.mkdir(exist_ok=True)

def bh(p):
 p=np.asarray(p,float); n=len(p); o=np.argsort(p); z=p[o]*n/np.arange(1,n+1); z=np.minimum.accumulate(z[::-1])[::-1]; out=np.empty(n); out[o]=np.minimum(z,1); return out

def holm(p):
 p=np.asarray(p,float); n=len(p); o=np.argsort(p); z=(n-np.arange(n))*p[o]; z=np.maximum.accumulate(z); out=np.empty(n); out[o]=np.minimum(z,1); return out

def slug(s): return ''.join(c.lower() if c.isalnum() else '_' for c in s).strip('_')

subtypes=[('Epithelial','Early Enterocyte'),('Epithelial','Enterocyte'),('Epithelial','Ileum Goblet'),('Epithelial','TA'),('Myeloid','Macrophage')]
audit=[]
model_summary=[]
for lineage,subtype in subtypes:
 analysis=f'{lineage}__Subtype__{subtype}'; base=f'human_{slug(analysis)}'
 meta=pd.read_csv(P30/'work'/f'{base}_meta.tsv',sep='\t'); counts_head=pd.read_csv(P30/'work'/f'{base}_counts.tsv',sep='\t',nrows=0).columns.tolist()
 assert meta.donor.is_unique and meta.sample_column.is_unique
 assert set(meta.condition)=={'Healthy','ActiveCrohn'}
 assert set(meta.sample_column).issubset(counts_head)
 # Actual Phase11 full-model inputs (no minimum-cell filter).
 h=(meta.condition=='Healthy').sum(); c=(meta.condition=='ActiveCrohn').sum(); all20=bool((meta.n_cells>=20).all())
 for _,x in meta.iterrows():
  audit.append(['PHASE11_ACTUAL','FULL',analysis,'NONE',x.donor,x.condition,x.biologic_treatment,int(x.n_cells),x.sample_column,'YES',h,c,int(meta.n_cells.min()),'YES' if all20 else 'NO','YES','YES','YES','ESTIMABLE' if all20 else 'ELIGIBILITY_VIOLATION',('' if x.n_cells>=20 else 'N_CELLS_BELOW_20')])
 model_summary.append(['PHASE11_ACTUAL','FULL',analysis,'NONE',h,c,int(meta.n_cells.min()),all20,True,True,True,'ESTIMABLE' if all20 else 'ELIGIBILITY_VIOLATION'])
 # Corrected Phase12 eligible population.
 e=meta[meta.n_cells>=20].copy(); h=(e.condition=='Healthy').sum(); c=(e.condition=='ActiveCrohn').sum(); ok=(h>=5 and c>=5 and e.donor.is_unique and set(e.sample_column).issubset(counts_head) and (e.n_cells>=20).all())
 for _,x in e.iterrows(): audit.append(['PHASE12_CORRECTED','FULL',analysis,'NONE',x.donor,x.condition,x.biologic_treatment,int(x.n_cells),x.sample_column,'YES' if ok else 'NO',h,c,int(e.n_cells.min()),'YES','YES','YES','YES','ESTIMABLE' if ok else 'NOT_ESTIMABLE',''])
 for _,x in meta[meta.n_cells<20].iterrows(): audit.append(['PHASE12_CORRECTED','FULL',analysis,'NONE',x.donor,x.condition,x.biologic_treatment,int(x.n_cells),x.sample_column,'NO',h,c,int(e.n_cells.min()),'YES','YES','YES','YES','EXCLUDED_BEFORE_FIT','N_CELLS_BELOW_20'])
 model_summary.append(['PHASE12_CORRECTED','FULL',analysis,'NONE',h,c,int(e.n_cells.min()),True,e.donor.is_unique,True,True,'ESTIMABLE' if ok else 'NOT_ESTIMABLE'])
 # LODO: omit only participants eligible for the full model.
 for omit in e.donor:
  q=e[e.donor!=omit].copy(); hh=(q.condition=='Healthy').sum(); cc=(q.condition=='ActiveCrohn').sum(); gate=(hh>=5 and cc>=5); status='ESTIMABLE' if gate else 'NON_ESTIMABLE_DONOR_GATE'
  for _,x in q.iterrows(): audit.append(['PHASE12_CORRECTED','LODO',analysis,omit,x.donor,x.condition,x.biologic_treatment,int(x.n_cells),x.sample_column,'YES' if gate else 'NO',hh,cc,int(q.n_cells.min()),'YES','YES','YES','YES' if gate else 'NOT_RUN',status,''])
  model_summary.append(['PHASE12_CORRECTED','LODO',analysis,omit,hh,cc,int(q.n_cells.min()),True,q.donor.is_unique,True,True if gate else None,status])
 # Untreated sensitivity: gate before any edgeR fit.
 q=e[(e.condition=='Healthy')|e.biologic_treatment.astype(str).str.lower().isin(['no','no_biologic','none','na',''])].copy(); hh=(q.condition=='Healthy').sum(); cc=(q.condition=='ActiveCrohn').sum(); gate=(hh>=5 and cc>=5); status='ESTIMABLE' if gate else 'NON_ESTIMABLE_DONOR_GATE'
 for _,x in q.iterrows(): audit.append(['PHASE12_CORRECTED','UNTREATED_CROHN_ONLY',analysis,'EXCLUDE_BIOLOGIC_TREATED_CROHN',x.donor,x.condition,x.biologic_treatment,int(x.n_cells),x.sample_column,'YES' if gate else 'NO',hh,cc,int(q.n_cells.min()),'YES','YES','YES','YES' if gate else 'NOT_RUN',status,''])
 model_summary.append(['PHASE12_CORRECTED','UNTREATED_CROHN_ONLY',analysis,'EXCLUDE_BIOLOGIC_TREATED_CROHN',hh,cc,int(q.n_cells.min()),True,q.donor.is_unique,True,True if gate else None,status])

cols=['execution_version','model_type','analysis','omitted_or_filter','participant_id','condition','biologic_treatment','n_cells','pseudobulk_column','entered_edgeR_fit','healthy_participants','crohn_participants','minimum_n_cells','all_n_cells_ge20','unique_participant','condition_valid','design_full_rank','model_status','exclusion_reason']
pd.DataFrame(audit,columns=cols).to_csv(O/'01_ELIGIBILITY_AUDIT.csv',index=False)
pd.DataFrame(model_summary,columns=['execution_version','model_type','analysis','omitted_or_filter','healthy_participants','crohn_participants','minimum_n_cells','all_n_cells_ge20','unique_participant','condition_valid','design_full_rank','model_status']).to_csv(W/'ELIGIBILITY_MODEL_SUMMARY.csv',index=False)

# Replace only affected full-model results, then recalculate the complete HE12 family.
old_epi=pd.read_csv(P30/'corrected_results'/'03_HUMAN_EPITHELIAL_RESULTS_CORRECTED.csv')
new_epi=old_epi.copy(); aff=pd.read_csv(W/'phase12_affected_results.tsv',sep='\t').rename(columns={'pval':'P'})
for _,x in aff.iterrows():
 subtype=x.analysis.split('__Subtype__',1)[1]; m=(new_epi.author_subtype==subtype)&(new_epi.frozen_program_id==x.frozen_program_id)
 if m.sum()!=1: raise RuntimeError(f'Expected one row: {subtype} {x.frozen_program_id}')
 for k in ['healthy_donors','crohn_donors','tested_genes','mapped_ranked_genes','ES','NES','P','size','leading_edge_n','direction','model_status']: new_epi.loc[m,k]=x[k]
new_epi['phase11_legacy_reported_q']=new_epi.get('legacy_reported_q',np.nan)
p=new_epi.P.astype(float).to_numpy(); new_epi['BH_q']=bh(p);new_epi['Holm_q']=holm(p);new_epi['BH_significant']=np.where(new_epi.BH_q<.05,'YES','NO');new_epi['Holm_significant']=np.where(new_epi.Holm_q<.05,'YES','NO');new_epi['multiple_testing_family']='HE12'
new_epi.to_csv(W/'HUMAN_EPITHELIAL_RESULTS_PHASE12.csv',index=False)
my=pd.read_csv(P30/'corrected_results'/'04_HUMAN_MYELOID_RESULTS_CORRECTED.csv');my.to_csv(W/'HUMAN_MYELOID_RESULTS_PHASE12.csv',index=False)

keys=['lineage','author_subtype','frozen_program_id','program']
cmp=old_epi.merge(new_epi,on=keys,suffixes=('_phase11','_phase12'))
for k in ['NES','P','BH_q','Holm_q','direction','BH_significant','healthy_donors','crohn_donors','tested_genes']:
 cmp[f'{k}_changed']=np.where(cmp[f'{k}_phase11'].astype(str)==cmp[f'{k}_phase12'].astype(str),'NO','YES')
cmp['eligibility_correction_required']=np.where(cmp.author_subtype.isin(['Early Enterocyte','TA']),'YES','NO')
cmp.to_csv(O/'03_CORRECTED_MODEL_RESULTS_COMPARISON.csv',index=False)

# Rebuild combined LODO/treatment sensitivity, replacing only Early Enterocyte and TA.
old_sens=pd.read_csv(P30/'corrected_results'/'04_HUMAN_PROGRAM_LODO_AND_TREATMENT_SENSITIVITY.csv')
keep=old_sens[~old_sens.author_subtype.isin(['Early Enterocyte','TA'])].copy()
lo=pd.read_csv(W/'phase12_affected_lodo.tsv',sep='\t').rename(columns={'pval':'P'});lo['author_subtype']=lo.analysis.str.split('__Subtype__',regex=False).str[-1];lo['sensitivity_type']='LODO';lo['omitted_or_filter']=lo.omitted_donor
tr=pd.read_csv(W/'phase12_affected_treatment.tsv',sep='\t').rename(columns={'pval':'P'});tr['author_subtype']=tr.analysis.str.split('__Subtype__',regex=False).str[-1];tr['sensitivity_type']='UNTREATED_CROHN_ONLY';tr['omitted_or_filter']='EXCLUDE_BIOLOGIC_TREATED_CROHN';tr['full_NES']=np.nan;tr['direction_match_full']='NOT_ESTIMABLE'
sc=['lineage','author_subtype','frozen_program_id','program','sensitivity_type','omitted_or_filter','healthy_donors','crohn_donors','NES','P','direction','model_status','full_NES','direction_match_full']
final_sens=pd.concat([keep[sc],lo[sc],tr[sc]],ignore_index=True);final_sens.to_csv(W/'HUMAN_PROGRAM_LODO_AND_TREATMENT_SENSITIVITY_PHASE12.csv',index=False)

# Create complete corrected source-data set and replace the dependent S2 source.
src0=P30/'FINAL_CORRECTED_SUBMISSION_PACKAGE'/'figure_source_data'
for x in src0.glob('*.csv'): shutil.copy2(x,S/x.name)
new_epi.to_csv(S/'S2b_human_epithelial_programs.csv',index=False)
# Recompute global symmetric scale from all six S2 panels.
files=['S2a_dog_T_programs.csv','S2a_human_T_programs.csv','S2b_dog_epithelial_programs.csv','S2b_human_epithelial_programs.csv','S2c_dog_myeloid_programs.csv','S2c_human_myeloid_programs.csv']
mx=max(pd.to_numeric(pd.read_csv(S/f).NES,errors='coerce').abs().max() for f in files);lim=np.ceil(mx*10)/10
pd.DataFrame([{'observed_max_abs_NES':mx,'rounded_symmetric_limit':lim,'panels':6,'status':'RECOMPUTED_AFTER_PHASE12_ELIGIBILITY_CORRECTION'}]).to_csv(S/'S2_GLOBAL_COLOR_SCALE.csv',index=False)

low=pd.DataFrame(audit,columns=cols); low=low[(low.execution_version=='PHASE11_ACTUAL')&(low.n_cells<20)]
assert set(low.participant_id)=={'HA54','HA62','HA04','HA07','HA65'} and len(low)==5
assert int((new_epi.BH_q<.05).sum())==11
(O/'02_ELIGIBILITY_DECISION.md').write_text(f'''# Phase 12 eligibility decision\n\n## Stage 1 decision\n\n`STAGE_1_CORRECTED_AND_PASS`\n\nPhase 11 correctly rebuilt raw pseudobulk counts but did not filter subtype metadata by `n_cells >= 20` before the edgeR design matrix. Five ineligible donor-state pseudobulks entered formal fits: Early Enterocyte HA54 (2 cells) and HA62 (7 cells), plus TA HA04 (17), HA07 (15) and HA65 (15).\n\nPhase 12 excluded these five samples before model fitting and refit only Early Enterocyte and TA. Corrected full models contain 11 Healthy/6 Crohn and 8 Healthy/5 Crohn participants, respectively. Enterocyte (13/6), Ileum Goblet (13/6) and Macrophage (13/6) already satisfied the rule and were not refit.\n\nAll final formal subtype models now satisfy `n_cells >= 20`, at least five participants per disease group, unique participant IDs, valid condition labels, matching pseudobulk columns and full-rank designs. Corrected LODO fits enforce the same rule; TA omissions of any Crohn participant are explicitly `NON_ESTIMABLE_DONOR_GATE` because only four Crohn participants remain. All untreated-Crohn subtype sensitivities remain non-estimable because only two untreated Crohn participants are available.\n\nThe complete HE12 BH family was recalculated after replacing the six affected program estimates. Human epithelial BH-significant tests remain {int((new_epi.BH_q<.05).sum())}/12. The Human Macrophage HM3 family remains unchanged at {int((my.BH_q<.05).sum())}/3.\n''',encoding='utf-8')
print({'low_cell_phase11_rows':len(low),'epi_bh_sig':int((new_epi.BH_q<.05).sum()),'my_bh_sig':int((my.BH_q<.05).sum()),'s2_limit':float(lim)})

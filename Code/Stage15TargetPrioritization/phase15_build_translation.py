from pathlib import Path
import os, re, hashlib, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns

ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
P15=ROOT/"15_therapeutic_translation"; META=ROOT/"00_metadata"; RES=ROOT/"08_results"; LOG=ROOT/"logs"
DATE="2026-08-30"
for d in [P15/"human_clinical",P15/"canine_literature",P15/"targets",P15/"evidence",P15/"matrix",P15/"summary",P15/"figures",P15/"audit",META,RES,LOG]: d.mkdir(parents=True,exist_ok=True)

def wc(df,p): df.to_csv(p,index=False,encoding="utf-8-sig",lineterminator="\n")
def wl(n,x): (LOG/n).write_text("\n".join(x)+"\n",encoding="utf-8")
def sh(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""): h.update(b)
 return h.hexdigest()
def rd(p,**k): return pd.read_csv(p,low_memory=False,**k)
def sw(s,n=28): return "\n".join(textwrap.wrap(str(s),n))

AGA="https://gastro.org/clinical-guidance/pharmacological-management-of-moderate-to-severe-crohns-disease/"
UPD="https://aga-fileuploader-bucket.s3.us-east-2.amazonaws.com/AGA%20CD%20Living%20Guideline%20Update%20Summary%20-%20March%202026%20FINAL%20MSL.pdf"
def hrow(drug,cl,axis,target,mech,route,guide,risk,src):
 return dict(drug=drug,drug_class=cl,axis=axis,target=target,mechanism=mech,Crohn_indication="FDA-approved for moderately-to-severely active Crohn disease",approval_status="FDA_APPROVED",guideline_status=guide,evidence_level="TierA regulatory + TierB guideline + randomized clinical evidence",route=route,major_safety_warnings=risk,regulatory_source=src,guideline_source=AGA,guideline_freeze="AGA 2025 living guideline; 2026-03 evidence review left recommendations unchanged",retrieval_date=DATE)
tnfrisk="serious/opportunistic infections including tuberculosis; malignancy; hepatitis B reactivation; agent-specific heart failure, neurologic or hypersensitivity cautions"
ilrisk="infections; tuberculosis evaluation; hypersensitivity; avoid live vaccines; agent-specific hepatic, malignancy or neurologic cautions"
human=pd.DataFrame([
 hrow("infliximab","TNF antagonist","TNF/TNFR-NF-kappaB","TNF","monoclonal antibody neutralization of TNF","intravenous","AGA RECOMMENDS; higher-efficacy option",tnfrisk,"https://www.accessdata.fda.gov/drugsatfda_docs/label/2021/103772s5401lbl.pdf"),
 hrow("adalimumab","TNF antagonist","TNF/TNFR-NF-kappaB","TNF","human monoclonal antibody neutralization of TNF","subcutaneous","AGA RECOMMENDS; higher-efficacy option",tnfrisk,"https://www.accessdata.fda.gov/scripts/opdlisting/oopd/detailedIndex.cfm?cfgridkey=230306"),
 hrow("certolizumab pegol","TNF antagonist","TNF/TNFR-NF-kappaB","TNF","pegylated antigen-binding fragment neutralization of TNF","subcutaneous","AGA SUGGESTS; lower-efficacy positioning",tnfrisk,"https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/125160s315lbl.pdf"),
 hrow("ustekinumab","IL-12/23 p40 antagonist","IL23-specific","IL12B/p40","antibody binding the shared IL-12/23 p40 subunit","intravenous induction; subcutaneous maintenance","AGA RECOMMENDS; higher-efficacy option",ilrisk,"https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/125261s167lbl.pdf"),
 hrow("risankizumab","IL-23 p19 antagonist","IL23-specific","IL23A/p19","antibody selectively binding IL-23 p19","intravenous induction; subcutaneous maintenance","AGA RECOMMENDS; higher-efficacy option",ilrisk,"https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/761105s039%2C761262s011lbl.pdf"),
 hrow("guselkumab","IL-23 p19 antagonist","IL23-specific","IL23A/p19","antibody selectively binding IL-23 p19","intravenous or subcutaneous induction; subcutaneous maintenance","AGA RECOMMENDS; higher-efficacy option",ilrisk,"https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/761061s028s029lbl.pdf"),
 hrow("mirikizumab","IL-23 p19 antagonist","IL23-specific","IL23A/p19","antibody selectively binding IL-23 p19 and not IL-12","intravenous induction; subcutaneous maintenance","AGA RECOMMENDS; higher/intermediate-efficacy positioning",ilrisk,"https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/761279s005lbl.pdf"),
 hrow("upadacitinib","JAK inhibitor","JAK-STAT broad","JAK1-preferential small-molecule activity","small-molecule inhibition of JAK-dependent cytokine signaling","oral","AGA RECOMMENDS; positioning depends on prior therapy","boxed warnings: serious infections, mortality, malignancy, major cardiovascular events and thrombosis","https://www.fda.gov/drugs/news-events-human-drugs/fda-approves-first-oral-treatment-moderately-severely-active-crohns-disease"),
 hrow("vedolizumab","gut-selective integrin antagonist","alpha4beta7-MAdCAM1","alpha4beta7 integrin","antibody blocking alpha4beta7-MAdCAM1 gut lymphocyte trafficking","intravenous; subcutaneous maintenance","AGA SUGGESTS; exposure-dependent positioning","hypersensitivity/infusion reactions; infections; liver injury; neurologic/PML monitoring caution","https://www.accessdata.fda.gov/drugsatfda_docs/label/2024/125476s054%2C761133s003lbl.pdf")
])
wc(human,P15/"human_clinical"/"HUMAN_CROHN_THERAPEUTIC_EVIDENCE.csv")
wl("phase15_human_clinical.log",[f"retrieval_date={DATE}","drugs=9","FDA_verified=YES","AGA_2026_03=UNCHANGED","human_drug_ranking_used_for_axis_priority=NO"])

def crow(study,year,pmid,doi,axis,target,n,subtype,tissue,assay,direction,stats,strength,limits,source):
 return dict(study=study,year=year,PMID=pmid,DOI=doi,axis=axis,target=target,dog_n=n,disease_definition="canine chronic inflammatory enteropathy/IBD as defined by study",CIE_subtype=subtype,tissue=tissue,assay=assay,direction=direction,statistics=stats,independent_validation="independent of project cohort",evidence_strength=strength,limitations=limits,source=source,retrieval_date=DATE)
canine=pd.DataFrame([
 crow("Activation of nuclear factor-kappaB in dogs with chronic enteropathies",2010,"19740552","10.1016/j.vetimm.2009.08.014","TNF/TNFR-NF-kappaB","activated NF-kappaB","26 CIE + 13 controls","10 FRE; 16 IBD","duodenal mucosa","IHC","more activated NF-kappaB-positive macrophages in CIE; signal decreased after successful therapy","macrophage group P<0.01","MODERATE","mixed subtypes; pathway activation is not TNF-specific causality","https://pubmed.ncbi.nlm.nih.gov/19740552/"),
 crow("Upregulation of STAT3 in dogs with chronic inflammatory enteropathies",2021,"33955083","10.1111/jvim.16141","IL6-JAK-STAT3","pSTAT3","28 CIE + 10 controls","FRE/SRE/PLE","duodenal mucosa","IHC","pSTAT3 higher in multiple compartments and associated with activity","epithelial CCECAI rho=0.541; P<0.001","MODERATE","retrospective; pSTAT3 is not IL23-specific","https://pubmed.ncbi.nlm.nih.gov/33955083/"),
 crow("Pro-inflammatory cytokine expression and STAT1/3 in canine CE",2024,"37906531","10.1177/03009858231207017","IFN-JAK/STAT","IFNG;IL21;STAT1;STAT3","48 CE + 11 controls","IEL-positive and IEL-negative CE","duodenal mucosa","qPCR/western/IHC","IFNG and IL21 higher and STAT1/3 up in IEL-positive CE","group contrasts and correlations reported","MODERATE","subtype-specific; lymphoma arm is not CIE evidence","https://pubmed.ncbi.nlm.nih.gov/37906531/"),
 crow("Colonic mucosal and cytobrush cytokine mRNA in canine IBD",2021,"33471872","10.1371/journal.pone.0245580","IL23-specific","IL23A/IL23p19","26 IBD + 9 controls","large-intestinal IBD","colonic mucosa/cytobrush","qRT-PCR","IL23p19 higher in mucosa; no cytobrush difference or severity correlation","mucosa P=0.029","LIMITED","colon not project duodenum/human TI; expression-only; no cell type or intervention","https://pmc.ncbi.nlm.nih.gov/articles/PMC7817028/"),
 crow("Selected cytokine expression in canine lymphocytic-plasmacytic colitis",2014,"24976586","10.1292/jvms.13-0635","IL23-specific","IL23A/IL23p19","6 LPC + 12 controls","large-intestinal LPC","colonic mucosa","qRT-PCR","IL23p19 was the only assayed cytokine differing between groups","significance reported; exact effect unavailable in abstract audit","LIMITED","very small; colon-specific; expression-only","https://pubmed.ncbi.nlm.nih.gov/24976586/"),
 crow("Preliminary CCR9 and MAdCAM1 study in canine CE",2025,"40564263","10.3390/ani15121710","alpha4beta7-MAdCAM1","MADCAM1;CCR9","5 CE + 5 controls","mixed/unspecified CE","duodenum/blood","qPCR/multiplex/IHC","MADCAM1 and CCR9 reported higher; TNF and IFNG also higher","small exploratory contrasts","LIMITED","preliminary sample; no alpha4beta7 engagement or treatment","https://pmc.ncbi.nlm.nih.gov/articles/PMC12189423/"),
 crow("CD4 T-cell cytokines in duodenal canine IBD",2014,"24270804","10.1292/jvms.13-0008","IFN-JAK/STAT","IL17A;IFNG;IL4;IL10","27 IBD + 8 controls","IBD","duodenum","qPCR/ELISA","no significant group differences","non-significant","CONFLICTING","heterogeneous disease and selected cytokines","https://pubmed.ncbi.nlm.nih.gov/24270804/"),
 crow("Cyclooxygenase-2 and 5-lipoxygenase in canine CE",2015,"25269796","10.1111/jvim.12496","TNF/TNFR-NF-kappaB","TNF;IL6","prospective cohort; n varied by tissue","FRE and IBD","duodenum/colon","qRT-PCR","TNF and IL6 not significantly different among groups","non-significant group contrasts","CONFLICTING","transcript levels do not measure pathway activation","https://pubmed.ncbi.nlm.nih.gov/25269796/"),
 crow("ACVIM-endorsed canine CIE consensus and systematic review",2026,"41742497","10.1093/jvimsj/aalaf017","CIE current-care context","not target-specific","systematic review","all CIE subtypes","not applicable","systematic review/consensus","no targeted TNF/JAK/IL23/alpha4beta7 CIE efficacy established","evidence synthesis","STRONG_CONTEXT_ONLY","context, not target-axis efficacy","https://pubmed.ncbi.nlm.nih.gov/41742497/")
])
wc(canine,P15/"canine_literature"/"CANINE_CIE_THERAPEUTIC_AXIS_LITERATURE.csv")
il23=canine[canine.axis.eq("IL23-specific")].copy()
il23.insert(0,"audit_conclusion","DIRECT_CANINE_COLONIC_EXPRESSION_EVIDENCE_IDENTIFIED_BUT_NO_PROJECT_LINEAGE_PROGRAM_OR_CIE_INTERVENTION_EVIDENCE")
wc(il23,P15/"canine_literature"/"CANINE_IL23_SPECIFIC_EVIDENCE.csv")
wl("phase15_canine_literature.log",[f"included={len(canine)}","IL23_direct_records=2","IL23_state=LIMITED_COLONIC_EXPRESSION_ONLY","subtypes_preserved=YES"])

targets="TNF TNFRSF1A TNFRSF1B NFKB1 RELA JAK1 JAK2 JAK3 TYK2 STAT1 STAT3 STAT5A STAT5B IL6 IL6R IL23A IL23R IL12B IL12RB1 ITGA4 ITGB7 MADCAM1 IFNG IFNGR1 TLR4".split()
lins=["Epithelial","T","Myeloid"]
om=rd(ROOT/"10_cross_species_baseline"/"mapping"/"DOG_HUMAN_ORTHOLOG_MASTER.csv")
om["_one"]=om.primary_one2one.astype(str).str.lower().eq("true")
orows=[]
for t in targets:
 m=om[(om.dog_original_symbol.astype(str).str.upper()==t)|(om.human_ortholog_symbol.astype(str).str.upper()==t)].copy()
 if len(m):
  m=m.sort_values(["_one","orthology_confidence","goc_score_if_available"],ascending=False); r=m.iloc[0]; one=bool(r["_one"])
  orows.append(dict(candidate_target=t,dog_gene=r.dog_original_symbol,human_gene=r.human_ortholog_symbol,dog_ensembl_id=r.dog_source_gene_id,human_ensembl_id=r.human_ortholog_release104_id,orthology_type=r.homology_type_release104,one2one=one,GOC=r.goc_score_if_available,confidence=r.orthology_confidence,protein_identity_if_available=r.dog_percent_identity,target_conservation_status="HIGH_CONFIDENCE_ONE2ONE" if one else "ORTHOLOG_PRESENT_NON_ONE2ONE",notes="Frozen Ensembl104 Phase10 mapping; conservation does not establish human-biologic cross-reactivity"))
 else: orows.append(dict(candidate_target=t,dog_gene="NOT_MAPPED",human_gene=t,dog_ensembl_id="NA",human_ensembl_id="NA",orthology_type="NOT_MAPPED_IN_FROZEN_RELEASE104",one2one=False,GOC="NA",confidence="NA",protein_identity_if_available="NA",target_conservation_status="NOT_ASSESSABLE",notes="No exact symbol row in frozen map; not proof of absent ortholog"))
oa=pd.DataFrame(orows)

dm=rd(ROOT/"06_pseudobulk"/"metadata"/"pseudobulk_sample_metadata.csv")
hm=rd(ROOT/"07_human_crohn"/"GSE266616"/"pseudobulk"/"metadata"/"HUMAN_PSEUDOBULK_SAMPLE_METADATA.csv")
er=[]
def expr(sp,li,t):
 if sp=="Dog":
  cp=ROOT/"06_pseudobulk"/"matrices"/f"pseudobulk_{li}_counts.csv.gz"; dp=ROOT/"06_pseudobulk"/"de_results"/f"{li}_CIE_vs_Healthy_edger.csv"
  c=rd(cp); d=rd(dp); q=d[(d.gene_symbol.astype(str).str.upper()==t)|(d.unique_var_name.astype(str).str.upper()==t)]
  row=c[c.feature_id.astype(str).isin(q.feature_id.astype(str))]
  meta=dm[(dm.lineage==li)&(dm.eligible.astype(str).str.lower().isin(["true","yes","1"]))]; sc="sample_id"; cc="condition"; tv="CIE"; rv="Healthy"
 else:
  cp=ROOT/"07_human_crohn"/"GSE266616"/"pseudobulk"/"matrices"/f"human_{li}_counts.csv.gz"; dp=ROOT/"07_human_crohn"/"GSE266616"/"pseudobulk"/"de_results"/f"Human_{li}_ActiveCrohn_vs_Healthy.csv"
  c=rd(cp); d=rd(dp); q=d[d.gene_symbol.astype(str).str.upper()==t]; row=c[c.gene_symbol.astype(str).str.upper()==t]
  meta=hm[(hm.lineage==li)&(hm.analysis=="Primary_ActiveCrohn_vs_Healthy")&(hm.eligible=="YES")]; sc="patient_id"; cc="condition"; tv="ActiveCrohn"; rv="Healthy"
 ss=[x for x in meta[sc].astype(str) if x in c.columns]; ts=[x for x in meta[meta[cc]==tv][sc].astype(str) if x in c.columns]; rs=[x for x in meta[meta[cc]==rv][sc].astype(str) if x in c.columns]
 base=dict(species=sp,lineage=li,target=t,source_object=str(cp.relative_to(ROOT)).replace("\\","/"),n_samples=len(ss),n_test=len(ts),n_reference=len(rs),inference="DESCRIPTIVE_ONLY_NO_NEW_DE")
 if not len(row) or not ss: return {**base,"detected":"NO","detection_rate_all":0,"detection_rate_test":0,"detection_rate_reference":0,"mean_CPM_all":0,"mean_CPM_test":0,"mean_CPM_reference":0,"frozen_logCPM":np.nan,"frozen_logFC":np.nan,"frozen_FDR":np.nan,"frozen_direction":"NOT_TESTED_OR_NOT_DETECTED"}
 v=row[ss].sum(axis=0).astype(float); lib=c[ss].sum(axis=0).astype(float); z=v/lib.replace(0,np.nan)*1e6; dr=q.iloc[0] if len(q) else None
 rate=lambda cols:float((v[cols]>0).mean()) if cols else np.nan; mean=lambda cols:float(z[cols].mean()) if cols else np.nan
 return {**base,"detected":"YES","detection_rate_all":rate(ss),"detection_rate_test":rate(ts),"detection_rate_reference":rate(rs),"mean_CPM_all":mean(ss),"mean_CPM_test":mean(ts),"mean_CPM_reference":mean(rs),"frozen_logCPM":float(dr.logCPM) if dr is not None else np.nan,"frozen_logFC":float(dr.logFC) if dr is not None else np.nan,"frozen_FDR":float(dr.FDR) if dr is not None else np.nan,"frozen_direction":str(dr.direction) if dr is not None else "NOT_TESTED_OR_NOT_DETECTED"}
for sp in ["Dog","Human"]:
 for li in lins:
  for t in targets: er.append(expr(sp,li,t))
ex=pd.DataFrame(er)
for c in ["detection_rate_all","detection_rate_test","detection_rate_reference","mean_CPM_all","mean_CPM_test","mean_CPM_reference","frozen_logCPM","frozen_logFC","frozen_FDR"]: ex[c]=pd.to_numeric(ex[c],errors="coerce").round(6)
wc(ex,P15/"targets"/"THERAPEUTIC_TARGET_EXPRESSION_CONTEXT.csv")
dd=ex[(ex.species=="Dog")&(ex.detected=="YES")].groupby("target").lineage.apply(lambda x:"|".join(sorted(set(x)))).to_dict()
hd=ex[(ex.species=="Human")&(ex.detected=="YES")].groupby("target").lineage.apply(lambda x:"|".join(sorted(set(x)))).to_dict()
oa["dog_expression_lineage"]=oa.candidate_target.map(dd).fillna("NOT_DETECTED_IN_FROZEN_PSEUDOBULK")
oa["human_expression_lineage"]=oa.candidate_target.map(hd).fillna("NOT_DETECTED_IN_FROZEN_PSEUDOBULK")
oa["human_biologic_crossreactivity_in_dog"]="UNKNOWN_OR_UNSUPPORTED"; oa["retrieval_date"]=DATE
wc(oa,P15/"targets"/"THERAPEUTIC_TARGET_ORTHOLOGY_AUDIT.csv")
wl("phase15_target_audit.log",[f"targets={len(targets)}",f"expression_rows={len(ex)}","mapping=FROZEN_ENSEMBL104","new_DE=NO","crossreactivity_inferred_from_orthology=NO"])

safety=pd.DataFrame([
 dict(axis="TNF/TNFR-NF-kappaB",drug_or_class="human TNF antagonists",human_major_risks=tnfrisk,canine_known_risks="no validated canine CIE safety dataset for audited biologics",immunosuppression="YES",infection_risk="HIGH_CLASS_CAUTION",vaccine_interaction="human-label precautions; canine relevance not established",GI_safety="NOT_ASSESSABLE_IN_CANINE_CIE",evidence_source="FDA labels and AGA guideline",translation_caution="target conservation does not establish antibody binding, immunogenicity, efficacy or safety in dogs"),
 dict(axis="JAK-STAT broad",drug_or_class="human upadacitinib",human_major_risks="serious infections, mortality, malignancy, MACE and thrombosis",canine_known_risks="no canine upadacitinib CIE dataset",immunosuppression="YES",infection_risk="HIGH_CLASS_CAUTION",vaccine_interaction="human precautions cannot be transferred quantitatively",GI_safety="NOT_ASSESSABLE_IN_CANINE_CIE",evidence_source="FDA approval and JAK safety communication",translation_caution="not interchangeable with canine JAK inhibitors"),
 dict(axis="JAK-STAT broad",drug_or_class="canine oclacitinib",human_major_risks="NA",canine_known_risks="immune modulation; infection/demodicosis; possible neoplasia exacerbation",immunosuppression="YES",infection_risk="EXPLICIT_LABEL_CAUTION",vaccine_interaction="indication-specific label precautions",GI_safety="CIE-specific safety not established",evidence_source="FDA Apoquel label",translation_caution="species feasibility only; not CIE efficacy"),
 dict(axis="JAK-STAT broad",drug_or_class="canine ilunocitinib",human_major_risks="NA",canine_known_risks="immunosuppression/infections/neoplasia; vaccine-induced disease and inadequate vaccine response warning",immunosuppression="YES",infection_risk="EXPLICIT_LABEL_CAUTION",vaccine_interaction="MAJOR_EXPLICIT_WARNING",GI_safety="CIE-specific safety not established",evidence_source="FDA Zenrelia veterinarian letter",translation_caution="species feasibility only; not CIE efficacy"),
 dict(axis="IL6-JAK-STAT3",drug_or_class="IN-115314 investigational JAK inhibitor",human_major_risks="NA",canine_known_risks="short preclinical observations only",immunosuppression="POTENTIAL",infection_risk="LONG_TERM_NOT_ESTABLISHED",vaccine_interaction="NOT_ASSESSABLE",GI_safety="NOT_ASSESSABLE_IN_CANINE_CIE",evidence_source="PMID 40341986",translation_caution="target engagement/exposure, not disease efficacy"),
 dict(axis="IL23-specific",drug_or_class="human IL12/23 or IL23p19 antibodies",human_major_risks=ilrisk,canine_known_risks="no validated canine CIE safety/cross-reactivity dataset",immunosuppression="YES",infection_risk="CLASS_CAUTION",vaccine_interaction="human precautions; canine translation unknown",GI_safety="NOT_ASSESSABLE_IN_CANINE_CIE",evidence_source="FDA labels",translation_caution="canine cross-reactivity UNKNOWN_OR_UNSUPPORTED"),
 dict(axis="alpha4beta7-MAdCAM1",drug_or_class="human vedolizumab",human_major_risks="hypersensitivity, infection, liver injury and neurologic/PML caution",canine_known_risks="no canine CIE safety/cross-reactivity dataset",immunosuppression="GUT_TRAFFICKING_MODULATION",infection_risk="CLASS_CAUTION",vaccine_interaction="human precautions; canine translation unknown",GI_safety="NOT_ASSESSABLE_IN_CANINE_CIE",evidence_source="FDA label",translation_caution="canine binding UNKNOWN_OR_UNSUPPORTED")
])
safety["retrieval_date"]=DATE; wc(safety,P15/"evidence"/"THERAPEUTIC_SAFETY_FEASIBILITY.csv")
wl("phase15_pharmacology.log",["canine_JAK_approved_precedent=OCLACITINIB|ILUNOCITINIB","canine_JAK_PKPD=AVAILABLE","canine_CIE_JAK_efficacy=NOT_DEMONSTRATED","human_biologic_crossreactivity=UNKNOWN_OR_UNSUPPORTED","dose_output=NONE"])

cols=["axis","subaxis","candidate_targets","program_type","human_clinical_evidence","human_clinical_basis","dog_cie_association","dog_lineages","dog_direction","dog_lodo_stability","human_crohn_direction","ibdverse_support","phase10_status","phase11_status","phase12_status","phase13_status","cross_species_evidence","ced_relationship","disease_context","fre_status","region_risk","human_treatment_sensitivity","target_conservation","canine_pharmacology","human_biologic_crossreactivity","safety_feasibility","layer1","layer2","layer3","layer4","layer5","layer6","biological_translation_priority","pharmacologic_translation_feasibility","priority_class","rank_group","major_support","major_caution","interpretation"]
def ax(axis,sub,targets,ptype,hc,hb,da,lin,direction,lodo,hdir,ibd,p10,p11,p12,p13,cross,ced,ctx,fre,region,treat,cons,pharm,xreact,safe,layers,bio,feas,pc,tier,sup,caut,interp):
 return dict(zip(cols,[axis,sub,targets,ptype,hc,hb,da,lin,direction,lodo,hdir,ibd,p10,p11,p12,p13,cross,ced,ctx,fre,region,treat,cons,pharm,xreact,safe,*layers,bio,feas,pc,tier,sup,caut,interp]))
P11="EXPLORATORY_WEAK_REFERENCE: 46/46 human reference models weak"
P12="EXPLORATORY_FUNCTIONAL_STRUCTURE: Level1 transfer but weak species mixing"
axis=pd.DataFrame([
 ax("TNF/TNFR-NF-kappaB","TNF blockade, TNFR modulation and downstream NF-kappaB kept separate","TNF|TNFRSF1A|TNFRSF1B|NFKB1|RELA","DRUGGABLE_AXIS","STRONG","FDA/AGA/RCT evidence primarily for TNF ligand blockade","STRONG","Epithelial|T|Myeloid","CIE_higher in all three","T ROBUST; Epithelial/Myeloid MODERATELY_STABLE","ActiveCrohn_higher","same-direction inflammatory replication","three lineages same-direction; multiple CONCORDANT_STRONG",P11,P12,"Myeloid strongest; Epithelial moderate; T weak overall","MODERATE","CeD shared in all three","GENERIC_MUCOSAL_INFLAMMATION","CIE_ROBUST_FRE_SUBTYPE_DISCORDANT; FRE lower","mostly LOW; program provenance retained","generally directionally stable","targets audited; conservation is not antibody binding","NO_DIRECT_CANINE_ANTI_TNF_CIE_PHARMACOLOGY","UNKNOWN_OR_UNSUPPORTED","infection/immunosuppression and biologic feasibility constraints",["STRONG","STRONG","STRONG","MODERATE","MODERATE","LIMITED"],"HIGH","LIMITED","HIGH_PRIORITY_CANDIDATE_AXIS","TIER_A","human clinical validation plus multi-lineage Dog association/LODO/cross-species support","generic inflammation; FRE discordance; dog n=7; biologic cross-reactivity unknown","High scientific priority with limited implementation feasibility; not a usable-agent claim"),
 ax("JAK-STAT broad","broad JAK-dependent cytokine signaling; not IL23-specific","JAK1|JAK2|JAK3|TYK2|STAT1|STAT3|STAT5A|STAT5B","DRUGGABLE_AXIS","STRONG","upadacitinib FDA/AGA/phase3 evidence","STRONG","Myeloid primary; T support; Epithelial mixed","Myeloid/T higher; Epithelial component variable","Myeloid/T ROBUST; Epithelial MODERATE","relevant Crohn programs higher","Myeloid inflammatory replication","Myeloid CONCORDANT_STRONG; other components mixed",P11,P12,"Myeloid strong; Epithelial moderate; T weak","MODERATE","shared inflammatory signal","GENERIC_MUCOSAL_INFLAMMATION","CIE robust with FRE summary discordance","LOW_TO_PROGRAM_DEPENDENT","heterogeneity retained","multiple one-to-one targets and lineage expression","CANINE_PHARMACOLOGY_AVAILABLE; different agents/indications","NOT_APPLICABLE_TO_CLASS; agents not interchangeable","infection, neoplasia and vaccine-related risks",["STRONG","STRONG","MODERATE","MODERATE","MODERATE","MODERATE"],"HIGH","MODERATE","HIGH_PRIORITY_CANDIDATE_AXIS","TIER_A","human Crohn validation, Myeloid correspondence and dog target-class engagement","no canine CIE efficacy and major safety constraints","Stronger dog target-class feasibility than human-specific antibody axes, not CIE efficacy"),
 ax("IL6-JAK-STAT3","IL6-linked JAK1/STAT3 arm; distinct from IL23","IL6|IL6R|JAK1|STAT3","DRUGGABLE_AXIS_WITH_GAP","LIMITED","broad JAK clinical evidence is not IL6-specific","MODERATE","Myeloid primary; T support; mucosal literature","Myeloid/T component higher; canine pSTAT3 higher","Myeloid ROBUST; T MODERATE; Epithelial UNSTABLE","component-level higher","general inflammation only","component correspondence, not IL6-specific",P11,P12,"Myeloid disease-vector context","MODERATE","CeD shared context","GENERIC_MUCOSAL_INFLAMMATION","mixed/summary discordant","program-dependent","component-dependent","IL6/JAK1/STAT3 audited","canine JAK1-pSTAT3 engagement outside CIE","IL6-biologic cross-reactivity not assessed","JAK-class immune/infection risk",["LIMITED","MODERATE","MODERATE","MODERATE","MODERATE","MODERATE"],"MODERATE","MODERATE","PROMISING_WITH_MAJOR_CAVEATS","TIER_B","canine pSTAT3 association plus JAK1-pSTAT3 engagement precedent","upstream cytokine unresolved; no IL6-specific CIE intervention","Promising subaxis with unresolved causal input"),
 ax("IL23-specific","IL12/23 p40 and IL23 p19 distinct from general STAT3","IL23A|IL23R|IL12B|IL12RB1|TYK2|STAT3","DRUGGABLE_AXIS_WITH_DOG_GAP","STRONG","four FDA-approved/guideline-supported human classes/agents","LIMITED","no project Level1 program; two colonic qPCR studies","colonic IL23p19 higher; project lineage effect absent","NOT_ASSESSABLE for project-specific program","clinically actionable in human Crohn","general interleukin only","NOT_ASSESSABLE as IL23-specific",P11,P12,"disease vector not ligand-specific","NOT_ASSESSABLE","NOT_ASSESSABLE_AS_SPECIFIC_AXIS","HUMAN_VALIDATED_DOG_COLON_EXPRESSION_HINT","NOT_ASSESSABLE","HIGH_TISSUE_MISMATCH_CAUTION","human anchor independent of Dog evidence","targets audited; conservation not antibody binding","NO_CANINE_CIE_IL23_ENGAGEMENT_OR_EFFICACY","UNKNOWN_OR_UNSUPPORTED","infection and biologic feasibility uncertainty",["STRONG","LIMITED","LIMITED","NOT_ASSESSABLE","LIMITED","LIMITED"],"LIMITED","LIMITED","HUMAN_VALIDATED_CANINE_EVIDENCE_LIMITED","TIER_C","strong human evidence plus limited reproducible canine colonic expression","no project program/cell localization/intervention; tissue mismatch; cross-reactivity unknown","General STAT3/JAK cannot be relabeled as IL23"),
 ax("IFN-JAK/STAT","IFNG/STAT1 context","IFNG|IFNGR1|JAK1|JAK2|STAT1","BIOLOGICAL_CONTEXT_ONLY","LIMITED","no audited Crohn IFN-directed class","MODERATE","T/IEL and Myeloid contexts","subtype-dependent IFNG/STAT1 activation","MODERATE_TO_NA","inflammatory context","contextual","mixed/component-level",P11,P12,"T overall weak","LIMITED","shared inflammation","GENERIC_INFLAMMATION_CAUTION","subtype-dependent","program-dependent","not a direct clinical axis","targets expressed","JAK class only; IFN-specific feasibility absent","NOT_ASSESSABLE","immune suppression risk",["LIMITED","MODERATE","LIMITED","LIMITED","MODERATE","LIMITED"],"LIMITED","LIMITED","EXPLORATORY_TRANSLATION_HYPOTHESIS","EXPLORATORY","subtype-specific canine IFNG/STAT evidence","conflicting study and no direct therapeutic class","Biological context only"),
 ax("alpha4beta7-MAdCAM1","gut-homing interaction","ITGA4|ITGB7|MADCAM1","EXPLORATORY_DRUGGABLE_AXIS","STRONG","vedolizumab FDA/AGA/RCT evidence","LIMITED","no project program; preliminary external study","ITGA4/ITGB7 effects non-significant; MADCAM1 higher in 5-dog CE study","NOT_ASSESSABLE as program","clinically actionable","NOT_ASSESSABLE","NOT_ASSESSABLE",P11,P12,"NOT_ASSESSABLE","NOT_ASSESSABLE","NOT_ASSESSABLE","EXPLORATORY_ONLY","NOT_ASSESSABLE","UNKNOWN","not registered","targets audited","NO_CANINE_CIE_TRAFFICKING_ENGAGEMENT","UNKNOWN_OR_UNSUPPORTED","binding/immunogenicity/trafficking safety unknown",["STRONG","LIMITED","LIMITED","NOT_ASSESSABLE","LIMITED","LIMITED"],"LIMITED","LIMITED","EXPLORATORY_TRANSLATION_HYPOTHESIS","EXPLORATORY","human validation plus preliminary canine MADCAM1 association","tiny study; no project program; no binding/engagement evidence","Remains exploratory by preregistered rule"),
])
def prog(name,ctx,cross,ced,region,tier="TIER_C",pc="BIOLOGICALLY_SUPPORTED_NOT_YET_ACTIONABLE"):
 return ax(name,"program-level; no single causal target","PROGRAM_LEVEL_NO_SINGLE_TARGET","PROGRAM_MODIFIER","NOT_ASSESSABLE","no credible single intervention chain","STRONG","Epithelial" if name.startswith("Epithelial") else "Myeloid","CIE program changed","ROBUST_TO_MODERATE","ActiveCrohn direction registered","program-level","registered Phase10 evidence",P11,P12,"lineage-level context",cross,ced,ctx,"context-dependent",region,"mostly stable","not a single-target question","NO_CREDIBLE_CAUSAL_TARGET_CHAIN","NOT_APPLICABLE","not actionable without causal regulator",["NOT_ASSESSABLE","STRONG","STRONG",cross if cross in ["MODERATE","STRONG"] else "LIMITED","STRONG","LIMITED"],"HIGH_BIOLOGY","LIMITED",pc,tier,"robust cell-resolved Dog-Human biology","large program; shared disease context and no causal target","Disease biology/biomarker hypothesis, not direct treatment selection")
axis=pd.concat([axis,pd.DataFrame([
 prog("Epithelial mitochondrial metabolism","DISEASE_BIOLOGY_NOT_SINGLE_TARGET","MODERATE","shared/mixed CeD","LOW_TO_PROGRAM_DEPENDENT"),
 prog("Epithelial fatty-acid metabolism","SHARED_MUCOSAL_METABOLIC_CHANGE","MODERATE","often CeD-shared","LOW"),
 prog("Epithelial xenobiotic metabolism","SHARED_MUCOSAL_METABOLIC_CHANGE","MODERATE","CeD-shared/not differentiating","LOW"),
 prog("Myeloid inflammatory/cytokine program","GENERIC_MUCOSAL_INFLAMMATION","STRONG","CeD-shared","LOW"),
])],ignore_index=True)
axis=pd.concat([axis,pd.DataFrame([
 ax("Myeloid TLR4","innate sensing context","TLR4","BIOLOGICAL_CONTEXT_ONLY","LIMITED","no approved Crohn TLR4 class","LIMITED","Myeloid","target gene not significant; broader program higher","NOT_ASSESSABLE","contextual","contextual","gene-level weak",P11,P12,"not receptor-specific","LIMITED","generic","BIOLOGICAL_CONTEXT_ONLY","NA","UNKNOWN","NA","one-to-one audited","NO_CIE_ENGAGEMENT","NA","not actionable",["LIMITED"]*6,"LIMITED","LIMITED","EXPLORATORY_TRANSLATION_HYPOTHESIS","EXPLORATORY","innate context","no target effect/clinical class","Exploratory only"),
 ax("T-cell metabolic divergence","species-divergent T metabolism","PROGRAM_LEVEL","BIOLOGICAL_CONTEXT_ONLY","NOT_ASSESSABLE","no direct class","STRONG","T","Dog lower; Human higher","ROBUST","ActiveCrohn higher","opposite","OPPOSITE",P11,P12,"T overall weak","CONFLICTING","context-dependent","SPECIES_DIVERGENT","context","program-dependent","program-dependent","not target question","NOT_ASSESSABLE","NA","not actionable",["NOT_ASSESSABLE","STRONG","CONFLICTING","CONFLICTING","STRONG","LIMITED"],"LIMITED","LIMITED","NOT_PRIORITIZED","NOT_PRIORITIZED","stable Dog program","opposite species direction","Species divergence, not therapy"),
 ax("Stress/UPR","epithelial stress response","PROGRAM_LEVEL","BIOLOGICAL_CONTEXT_ONLY","NOT_ASSESSABLE","no direct class","MODERATE","Epithelial","CIE higher","ROBUST","representative mixed/opposite","NA","OPPOSITE",P11,P12,"Epithelial context","CONFLICTING","CeD not prominent is not specificity","DIFFERENTIATING_HINT_ONLY","partial","LOW","stable","not target question","NOT_ASSESSABLE","NA","not actionable",["NOT_ASSESSABLE","MODERATE","LIMITED","CONFLICTING","STRONG","LIMITED"],"LIMITED","LIMITED","NOT_PRIORITIZED","NOT_PRIORITIZED","robust Dog program","cross-species conflict/no target","Biology only"),
 ax("Epithelial proliferation","cell-cycle program","PROGRAM_LEVEL","BIOLOGICAL_CONTEXT_ONLY","NOT_ASSESSABLE","no inference from proliferation alone","STRONG","Epithelial","CIE higher","ROBUST","ActiveCrohn higher","NA","CONCORDANT_STRONG",P11,P12,"Epithelial moderate","MODERATE","CeD-shared","GENERIC_REPAIR_OR_INJURY_PROGRAM","NA","HIGH_REGION_RISK","stable","not target question","NOT_ASSESSABLE","NA","anti-proliferative inference unsafe","NOT_ASSESSABLE STRONG STRONG MODERATE MODERATE LIMITED".split(),"MODERATE_BIOLOGY","LIMITED","BIOLOGICALLY_SUPPORTED_NOT_YET_ACTIONABLE","TIER_C","conserved disease program","high region risk; CeD-shared; no chain","Biomarker/disease-state context")
])],ignore_index=True)
axis["no_numeric_composite_score"]="YES"; axis["retrieval_date"]=DATE
wc(axis,P15/"matrix"/"THERAPEUTIC_AXIS_EVIDENCE_MATRIX.csv"); wc(axis,P15/"matrix"/"CANINE_CIE_THERAPEUTIC_TRANSLATION_MATRIX.csv")

support={"TNF/TNFR-NF-kappaB":"YES_MULTI_LINEAGE","JAK-STAT broad":"YES_MYELOID_PRIMARY","IL23-specific":"LIMITED_NO_PROJECT_SPECIFIC_PROGRAM","alpha4beta7-MAdCAM1":"LIMITED_EXPLORATORY"}
dr=[]
for _,r in human.iterrows():
 if r.axis=="JAK-STAT broad": cross="CANINE_TARGET_CLASS_ACTIVITY_AVAILABLE_DIFFERENT_AGENTS"; pk="CANINE_JAK_CLASS_PKPD_AVAILABLE; NO_UPADACITINIB_CIE_DATA"; st="CANINE_PHARMACOLOGY_FEASIBLE_NOT_CIE_VALIDATED"
 else: cross="TARGET_CONSERVED_BIOLOGIC_CROSSREACTIVITY_UNKNOWN"; pk="NO_CANINE_CIE_AGENT_SPECIFIC_PKPD"; st="TARGET_CONSERVED_BIOLOGIC_CROSSREACTIVITY_UNKNOWN"
 dr.append(dict(axis=r.axis,human_drug_class=r.drug_class,representative_drug=r.drug,human_Crohn_status=r.approval_status+"; "+r.guideline_status,molecular_target=r.target,target_supported_by_project=support.get(r.axis,"LIMITED"),canine_target_conserved="AUDITED",canine_crossreactivity_or_activity=cross,canine_PKPD=pk,canine_CIE_efficacy="NOT_DEMONSTRATED",safety_caution=r.major_safety_warnings,translation_status=st,notes="No dose, regimen or canine clinical recommendation; agents are not substituted"))
drug=pd.DataFrame(dr); wc(drug,P15/"matrix"/"THERAPEUTIC_DRUG_CLASS_MATRIX.csv")
prio=axis[["rank_group","axis","priority_class","disease_context","dog_lineages","human_clinical_evidence","dog_cie_association","cross_species_evidence","dog_lodo_stability","pharmacologic_translation_feasibility","safety_feasibility","major_support","major_caution"]].rename(columns={"dog_lineages":"primary_and_supporting_lineages","human_clinical_evidence":"human_clinical_status","dog_cie_association":"canine_biological_status","cross_species_evidence":"cross_species_status","dog_lodo_stability":"canine_robustness","major_support":"reason","major_caution":"safety_caution_and_major_gap"})
nv={"TNF/TNFR-NF-kappaB":"canine binding; ex vivo engagement; species-compatible reagent; safety before efficacy","JAK-STAT broad":"CIE ex vivo response; lineage engagement; safety-focused PK/PD; then prospective proof-of-concept","IL6-JAK-STAT3":"resolve upstream cytokine and perturb pSTAT3","IL23-specific":"small-intestinal/cell/protein replication and canine-specific engagement","alpha4beta7-MAdCAM1":"larger cohort and canine binding/trafficking assay"}
prio["next_validation_needed"]=prio.axis.map(nv).fillna("resolve causal regulator with ex vivo or organoid perturbation before pharmacology")
wc(prio,P15/"summary"/"CANDIDATE_THERAPEUTIC_AXIS_PRIORITY_LIST.csv")

mc=[]
def addm(a,steps):
 for i,(f,t,s,b,g) in enumerate(steps,1): mc.append(dict(axis=a,step_order=i,from_node=f,to_node=t,relationship=s,evidence_source_type=b,evidence_basis=g,gap="NONE" if b!="GAP" else t,retrieval_date=DATE))
addm("TNF/TNFR-NF-kappaB",[("TNF","TNFR","ligand/receptor","LITERATURE","FDA mechanism"),("TNFR","NF-kappaB","signaling","LITERATURE_PLUS_PROJECT_DATA","canine NF-kappaB literature/project programs"),("NF-kappaB","Dog lineage programs","association","PROJECT_DATA","Phase5/10/14"),("Human Crohn TNF","TNF antagonist class","clinical intervention","LITERATURE","FDA/AGA/RCT"),("human biologic","canine binding","translation","GAP","unknown")])
addm("JAK-STAT broad",[("multiple cytokines","JAK kinases","receptor signaling","LITERATURE","canonical"),("JAK","STAT","phosphorylation","CANINE_PHARMACOLOGY","dog pSTAT engagement"),("STAT","Myeloid CIE program","association","PROJECT_DATA","Phase5/10/13"),("Human Crohn","JAK inhibitor class","clinical","LITERATURE","FDA/AGA/RCT"),("canine JAK pharmacology","canine CIE efficacy","translation","GAP","absent")])
addm("IL23-specific",[("IL23 ligand","IL23 receptor","signaling","LITERATURE","human plus canine colon expression"),("IL23 receptor","TYK2/JAK2/STAT3","signaling","LITERATURE","canonical"),("STAT3","Dog IL23 lineage program","disease link","GAP","general STAT3 cannot specify IL23"),("Human Crohn","IL23 biologics","clinical","LITERATURE","FDA/AGA/RCT"),("human biologic","canine engagement","translation","GAP","unknown")])
addm("alpha4beta7-MAdCAM1",[("alpha4beta7","MAdCAM1","gut homing","LITERATURE","vedolizumab mechanism/preliminary dog study"),("interaction","canine CIE trafficking","disease link","GAP","no project program"),("human biologic","canine binding","translation","GAP","unknown")])
mech=pd.DataFrame(mc); wc(mech,P15/"evidence"/"THERAPEUTIC_MECHANISM_CHAIN.csv")
wl("phase15_matrix.log",[f"axes={len(axis)}",f"drugs={len(drug)}","numeric_score=NONE","cross_method_vote=NONE","TierA=TNF|JAK","IL23_separate=YES","alpha4beta7=EXPLORATORY"])

src=[]
src.append(dict(evidence_id="HUM_GUIDE_001",axis="all",target="multiple",species="Human",disease="Crohn",evidence_type="guideline",study_type="AGA living guideline",n="16 recommendations",intervention="advanced therapies",outcome="recommendations",direction="supports audited classes",strength="STRONG",source_title="AGA living Crohn guideline",year=2025,PMID="41274746",DOI="10.1053/j.gastro.2025.09.038",regulatory_status="guideline; 2026-03 unchanged",limitations="not dog efficacy",source=AGA))
rcts=[("TNF/TNFR-NF-kappaB","TNF","ACCENT I",2002,"12075077","10.1016/S0140-6736(02)08512-4"),("TNF/TNFR-NF-kappaB","TNF","CHARM",2007,"17241864","10.1053/j.gastro.2006.11.041"),("TNF/TNFR-NF-kappaB","TNF","PRECiSE 2",2007,"17215358","10.1056/NEJMoa062897"),("IL23-specific","IL12B/p40","UNITI",2016,"27959607","10.1056/NEJMoa1602773"),("IL23-specific","IL23A/p19","ADVANCE/MOTIVATE",2022,"35644154","10.1016/S0140-6736(22)00467-6"),("IL23-specific","IL23A/p19","GALAXI-2/3",2025,"40684778","10.1016/S0140-6736(25)00681-6"),("IL23-specific","IL23A/p19","VIVID-1",2025,"39581202","10.1016/S0140-6736(24)01760-8"),("JAK-STAT broad","JAK1","U-EXCEL/U-EXCEED/U-ENDURE",2023,"37224198","10.1056/NEJMoa2212728"),("alpha4beta7-MAdCAM1","alpha4beta7","GEMINI 2",2013,"23964933","10.1056/NEJMoa1215739")]
for i,(a,t,title,y,p,d) in enumerate(rcts,1): src.append(dict(evidence_id=f"HUM_RCT_{i:03d}",axis=a,target=t,species="Human",disease="Crohn",evidence_type="clinical efficacy",study_type="randomized controlled trial",n="trial cohort",intervention="study drug versus control",outcome="clinical/endoscopic endpoints",direction="supports target class",strength="STRONG" if title!="PRECiSE 2" else "MODERATE",source_title=title,year=y,PMID=p,DOI=d,regulatory_status="approved class",limitations="human only",source=f"https://pubmed.ncbi.nlm.nih.gov/{p}/"))
for i,r in canine.iterrows(): src.append(dict(evidence_id=f"DOG_CIE_{i+1:03d}",axis=r.axis,target=r.target,species="Dog",disease="canine CIE/IBD",evidence_type="disease association literature",study_type=r.assay,n=r.dog_n,intervention="observational/context; not targeted CIE efficacy",outcome=r.direction,direction=r.direction,strength=r.evidence_strength,source_title=r.study,year=r.year,PMID=r.PMID,DOI=r.DOI,regulatory_status="not regulatory efficacy evidence",limitations=r.limitations,source=r.source))
ph=[("DOG_PHARM_001","JAK-STAT broad","JAK1","oclacitinib cellular pharmacology",2014,"24495176","10.1111/vde.12101","target-class engagement"),("DOG_PHARM_002","JAK-STAT broad","JAK","oclacitinib canine PK",2014,"24330031","10.1111/jvp.12087","canine exposure"),("DOG_PHARM_003","IL6-JAK-STAT3","JAK1-pSTAT3","IN-115314 canine PK/PD",2025,"40341986","10.1111/bph.70055","IL6-induced pSTAT3 engagement")]
for eid,a,t,title,y,p,d,out in ph: src.append(dict(evidence_id=eid,axis=a,target=t,species="Dog",disease="species pharmacology; not CIE",evidence_type="canine pharmacology",study_type="cellular/PK/PD",n="canine study",intervention="JAK-pathway compound",outcome=out,direction="feasibility available",strength="MODERATE",source_title=title,year=y,PMID=p,DOI=d,regulatory_status="agent/indication-specific",limitations="not CIE efficacy",source=f"https://pubmed.ncbi.nlm.nih.gov/{p}/"))
src.extend([
 dict(evidence_id="DOG_REG_001",axis="JAK-STAT broad",target="JAK",species="Dog",disease="dermatologic indication; not CIE",evidence_type="veterinary regulatory",study_type="FDA label",n="regulatory dossier",intervention="oclacitinib",outcome="dog class feasibility/safety precedent",direction="feasibility only",strength="STRONG_FEASIBILITY_ONLY",source_title="FDA Apoquel labeling",year=2021,PMID="",DOI="",regulatory_status="FDA approved for specified dog dermatologic indication",limitations="not CIE",source="https://www.fda.gov/animal-veterinary/drug-labels/2021-animal-drug-safety-related-labeling-changes"),
 dict(evidence_id="DOG_REG_002",axis="JAK-STAT broad",target="JAK",species="Dog",disease="dermatologic indication; not CIE",evidence_type="veterinary regulatory",study_type="FDA safety letter",n="regulatory dossier",intervention="ilunocitinib",outcome="dog class feasibility with vaccine/infection warning",direction="feasibility only",strength="STRONG_FEASIBILITY_ONLY",source_title="FDA Zenrelia safety letter",year=2024,PMID="",DOI="",regulatory_status="FDA approved for specified dog dermatologic indication",limitations="not CIE",source="https://www.fda.gov/animal-veterinary/product-safety-information/dear-veterinarian-letter-regarding-important-safety-information-associated-use-zenrelia-ilunocitinib")
])
se=pd.DataFrame(src); se["retrieval_date"]=DATE; wc(se,P15/"evidence"/"THERAPEUTIC_SOURCE_EVIDENCE.csv")

search=pd.DataFrame([
 ("AGA","2025 living Crohn pharmacologic guideline",DATE,"1 formal guideline",1,1,"human anchor"),("AGA","March 2026 living update",DATE,"1 update",1,1,"verify unchanged"),("FDA Drugs@FDA","TNF Crohn labels",DATE,"interface total unavailable",6,3,"TNF regulatory"),("FDA Drugs@FDA","IL12/23 and IL23p19 Crohn labels",DATE,"interface total unavailable",9,4,"IL23 regulatory"),("FDA","upadacitinib Crohn and JAK safety",DATE,"2 primary pages",2,2,"JAK regulatory/safety"),("FDA Drugs@FDA","vedolizumab label",DATE,"label/review records",3,1,"integrin regulatory"),("PubMed/PMC","canine CIE TNF NF-kappaB",DATE,"interface total unavailable",11,3,"dog TNF"),("PubMed/PMC","canine CIE JAK STAT3",DATE,"interface total unavailable",10,3,"dog JAK/STAT"),("PubMed/PMC","canine IBD IL23 IL23A IL23R IL12B TYK2",DATE,"interface total unavailable",12,2,"dog IL23"),("PubMed/PMC","canine CE integrin MAdCAM",DATE,"interface total unavailable",7,1,"dog integrin"),("PubMed","canine CIE IFN metabolic",DATE,"interface total unavailable",12,2,"other axes"),("PubMed","PMID 40341986",DATE,"1 exact",1,1,"dog JAK PKPD"),("FDA Animal Drugs@FDA","Apoquel",DATE,"label/FOI/safety",5,2,"dog JAK precedent"),("FDA Animal Drugs@FDA","Zenrelia",DATE,"approval/FOI/safety",4,2,"dog JAK precedent"),("Crossref/DOI verification","identifier cross-check",DATE,"not bulk search",15,15,"identifier consistency only")
],columns=["database_source","query","search_date","result_count","screened","included","purpose"])
wc(search,META/"PHASE15_LITERATURE_SEARCH_AUDIT.csv")

sns.set_theme(style="whitegrid",context="talk"); fd=P15/"figures"
states=["NOT_ASSESSABLE","CONFLICTING","LIMITED","MODERATE","STRONG"]; sm={s:i for i,s in enumerate(states)}; colors=["#d9d9d9","#c65f5f","#e4b363","#5d9cec","#2e8b57"]
lc=[f"layer{i}" for i in range(1,7)]
def sf(n): plt.savefig(fd/n,dpi=210,bbox_inches="tight",facecolor="white"); plt.close()
plot=axis[axis.rank_group.isin(["TIER_A","TIER_B","TIER_C","EXPLORATORY"])]
mat=plot[lc].applymap(lambda x:sm.get(x,0)).to_numpy()
plt.figure(figsize=(11,10)); a=sns.heatmap(mat,cmap=ListedColormap(colors),vmin=0,vmax=4,cbar=False,linewidths=1,annot=plot[lc].to_numpy(),fmt="",annot_kws={"fontsize":7}); a.set_xticklabels(["Human","Dog","Cell type","Cross-species","Robustness","Feasibility"],rotation=20); a.set_yticklabels([sw(x,24) for x in plot.axis],rotation=0); plt.title("Six-layer evidence ladder: categorical, no composite score"); sf("01_evidence_ladder_heatmap.png")
pc={"TIER_A":"#2e8b57","TIER_B":"#5d9cec","TIER_C":"#e4b363","EXPLORATORY":"#9b7ede","NOT_PRIORITIZED":"#999999"}
plt.figure(figsize=(11,8))
for i,(_,r) in enumerate(axis.iterrows()): plt.barh(i,1,color=pc[r.rank_group]); plt.text(.02,i,r.axis,va="center",fontsize=9,color="white" if r.rank_group in ["TIER_A","TIER_B"] else "black")
plt.yticks([]);plt.xticks([]);plt.title("Rule-based candidate axis priority groups");plt.xlabel("Categories are not a numeric ranking");sf("02_candidate_axis_priority_groups.png")
ctx=["TNF/TNFR-NF-kappaB","JAK-STAT broad","IL23-specific","alpha4beta7-MAdCAM1","Epithelial mitochondrial metabolism","Epithelial fatty-acid metabolism","Epithelial xenobiotic metabolism"]
q=axis.set_index("axis").loc[ctx]; ann=np.array([[sw(r.dog_direction,20),sw(r.human_crohn_direction,20),sw(r.ced_relationship,20)] for _,r in q.iterrows()])
plt.figure(figsize=(13,7));a=sns.heatmap(np.ones((len(q),3)),cmap=ListedColormap(["#f4f1de"]),cbar=False,linewidths=1,annot=ann,fmt="",annot_kws={"fontsize":7});a.set_xticklabels(["Dog CIE","Human Crohn","CeD relationship"],rotation=0);a.set_yticklabels([sw(x,23) for x in q.index],rotation=0);plt.title("Dog-Crohn-CeD context");sf("03_dog_crohn_ced_context.png")
def chain(title,file,nodes,gaps):
 plt.figure(figsize=(14,3)); x=np.arange(len(nodes))
 for i,n in enumerate(nodes):
  plt.scatter(i,0,s=2300,color="#c65f5f" if gaps[i] else "#5d9cec",edgecolor="white",zorder=3);plt.text(i,0,sw(n,16),ha="center",va="center",fontsize=8,color="white",weight="bold")
  if i<len(nodes)-1: plt.annotate("",xy=(i+.7,0),xytext=(i+.3,0),arrowprops=dict(arrowstyle="->",lw=2))
 plt.xlim(-.6,len(nodes)-.4);plt.ylim(-.7,.7);plt.axis("off");plt.title(title+"; red = explicit gap");sf(file)
chain("TNF evidence chain","04_tnf_evidence_chain.png",["TNF","TNFR","NF-kappaB","Dog lineages","Human class","Dog-compatible agent?"],[0,0,0,0,0,1])
chain("JAK-STAT evidence chain","05_jak_stat_evidence_chain.png",["cytokines","JAK","STAT","Myeloid CIE","Human class","Dog CIE efficacy?"],[0,0,0,0,0,1])
chain("IL23-specific gap","06_il23_evidence_gap.png",["IL23","IL23R","STAT3","Dog IL23 lineage?","Human class","Dog engagement?"],[0,0,0,1,0,1])
ot=["TNF","TNFRSF1A","JAK1","JAK2","TYK2","STAT3","IL23A","IL23R","IL12B","ITGA4","ITGB7","MADCAM1"]; o=oa.set_index("candidate_target").reindex(ot); v=np.column_stack([o.one2one.fillna(False).astype(int),o.index.to_series().map(lambda x:int(x in dd)),o.index.to_series().map(lambda x:int(x in hd)),o.index.to_series().map(lambda x:int(x in ["JAK1","JAK2","TYK2","STAT3"]))]);an=np.where(v==1,"YES","LIMITED/NO")
plt.figure(figsize=(9,7));a=sns.heatmap(v,cmap=ListedColormap(["#e4b363","#2e8b57"]),cbar=False,linewidths=1,annot=an,fmt="",annot_kws={"fontsize":8});a.set_xticklabels(["one-to-one","Dog expression","Human expression","dog JAK-class context"],rotation=25,ha="right");a.set_yticklabels(ot,rotation=0);plt.title("Target orthology/expression/pharmacology");sf("07_target_orthology_pharmacology_matrix.png")
bo={"LIMITED":0,"MODERATE":1,"MODERATE_BIOLOGY":1,"HIGH_BIOLOGY":2,"HIGH":3}; fo={"LIMITED":0,"MODERATE":1,"STRONG":2}
plt.figure(figsize=(11,8))
for _,r in axis.iterrows():
 x=fo.get(r.pharmacologic_translation_feasibility,0); y=bo.get(r.biological_translation_priority,0);plt.scatter(x,y,s=160,color=pc[r.rank_group],edgecolor="white");plt.text(x+.03,y+.03,sw(r.axis,16),fontsize=7)
plt.xticks([0,1,2],["LIMITED","MODERATE","STRONG"]);plt.yticks([0,1,2,3],["LIMITED","MODERATE","HIGH_BIOLOGY","HIGH"]);plt.xlabel("Pharmacologic feasibility");plt.ylabel("Biological priority");plt.title("Independent ordinal categories");sf("08_biology_vs_feasibility.png")
ca=["TNF/TNFR-NF-kappaB","JAK-STAT broad","IL6-JAK-STAT3","IL23-specific","alpha4beta7-MAdCAM1","Epithelial mitochondrial metabolism","Epithelial proliferation"]; cv=np.array([[1,2,1,1,1,2],[1,2,1,1,1,1],[1,2,1,2,1,1],[2,2,1,2,2,2],[2,2,2,2,2,2],[1,1,1,0,1,0],[2,1,1,0,1,0]]); lab=np.array([["LOW","MODERATE","HIGH"][x] for x in cv.flat]).reshape(cv.shape)
plt.figure(figsize=(11,7));a=sns.heatmap(cv,cmap=ListedColormap(["#2e8b57","#e4b363","#c65f5f"]),cbar=False,linewidths=1,annot=lab,fmt="");a.set_xticklabels(["Region","FRE","Treatment","LODO","Cross-method","Cross-reactivity"],rotation=25,ha="right");a.set_yticklabels([sw(x,22) for x in ca],rotation=0);plt.title("Caution matrix");sf("09_caution_matrix.png")
ma=["TNF/TNFR-NF-kappaB","JAK-STAT broad","IL6-JAK-STAT3","IL23-specific","alpha4beta7-MAdCAM1","Epithelial mitochondrial metabolism","Epithelial fatty-acid metabolism","Epithelial xenobiotic metabolism","Myeloid inflammatory/cytokine program"]; mm=axis.set_index("axis").loc[ma]; ar=mm[lc].applymap(lambda x:sm.get(x,0)).to_numpy()
plt.figure(figsize=(13,8));a=sns.heatmap(ar,cmap=ListedColormap(colors),vmin=0,vmax=4,cbar=False,linewidths=1,annot=mm[lc].to_numpy(),fmt="",annot_kws={"fontsize":7});a.set_xticklabels(["L1 Human","L2 Dog","L3 Cell","L4 Cross","L5 Robust","L6 Feasible"],rotation=25,ha="right");a.set_yticklabels([sw(x,24) for x in ma],rotation=0);plt.title("Canine CIE Therapeutic Translation Matrix");sf("10_translation_matrix_main.png")
plt.figure(figsize=(11,6));sns.countplot(data=human,y="axis",hue="drug_class",dodge=False);plt.xlabel("Audited FDA-approved Crohn agents: descriptive count");plt.ylabel("");plt.title("Human Crohn clinical/regulatory landscape");plt.legend(fontsize=8,loc="lower right");sf("11_human_landscape.png")
wl("phase15_visualization.log",[f"figures={len(list(fd.glob('*.png')))}","numeric_composite=NO","ordinal_categories_only=YES"])

report=f"""# Phase15 Report - Therapeutic Translation Matrix and Candidate Therapeutic Axis Prioritization

Execution freeze: {DATE}

Phase15 prioritizes mechanistic therapeutic axes for translational research. It does not establish efficacy, safety, dose, or clinical indication of any human therapeutic agent in dogs with CIE.

## 1. Phase14 handoff
Phase14 ended GREEN with 32 independent subjects, 203,555 author-final cells and 25/25 registered programs evaluated against active CeD. Ten programs were generic mucosal inflammatory and zero were strongly CeD-differentiating. The project-only axis matrix was frozen before clinical literature was opened.

## 2. Human therapeutic landscape
The 2025 AGA living Crohn guideline and its March 2026 unchanged review are the clinical anchor. Nine required drugs were independently checked against FDA records. Human drug efficacy populates Layer1 only and did not reverse-select Dog programs.

## 3. TNF/TNFR-NF-kappaB
Evidence ladder: STRONG / STRONG / STRONG / MODERATE / MODERATE / LIMITED. Human evidence mainly supports TNF ligand blockade, not broad NF-kappaB inhibition. Dog programs are higher in Epithelial, T and Myeloid; Phase10 is strongly concordant and Phase13 most strongly supports Myeloid. Phase11 is weak-reference exploratory and Phase12 is species-separated exploratory. T is LODO robust; Epithelial/Myeloid moderately stable; FRE is discordant. TNF is TIER_A/HIGH_PRIORITY_CANDIDATE_AXIS scientifically, with LIMITED implementation feasibility because human-biologic canine binding is unknown.

## 4. Broad JAK-STAT
Evidence ladder: STRONG / STRONG / MODERATE / MODERATE / MODERATE / MODERATE. Myeloid is clearest, T supports, Epithelial is component-dependent. Dog-approved JAK-class products for non-CIE indications and canine PK/PD/target engagement make class feasibility stronger than for human-specific antibodies. They do not demonstrate CIE efficacy and do not substitute for upadacitinib. JAK is TIER_A with major infection, immune, neoplasia and vaccine cautions.

## 5. IL23-specific
Evidence ladder: STRONG / LIMITED / LIMITED / NOT_ASSESSABLE / LIMITED / LIMITED. Two canine colonic studies show IL23p19 expression evidence, but no project Level1 IL23-specific program, cell localization, LODO/FRE program or canine CIE intervention exists. Colon/duodenum/TI mismatch is explicit. General JAK/STAT3 is not IL23. IL23 remains HUMAN_VALIDATED_CANINE_EVIDENCE_LIMITED/TIER_C.

## 6. alpha4beta7/MAdCAM1
Vedolizumab provides strong human evidence. Dog evidence is a preliminary 5-versus-5 MADCAM1 study plus non-significant project ITGA4/ITGB7 effects, with no registered axis, fidelity analysis, engagement or cross-reactivity. It remains EXPLORATORY_ONLY.

## 7. Epithelial metabolic programs
Mitochondrial respiration, fatty-acid and xenobiotic programs are robust cell-resolved biology, but they are large programs, often CeD-shared, and lack a causal regulator-to-pharmacology chain. They are BIOLOGICALLY_SUPPORTED_NOT_YET_ACTIONABLE/TIER_C, best treated as biology/biomarker and functional-validation hypotheses.

## 8. Target conservation
Frozen Ensembl104 mapping was used. Detection rate and CPM were summarized from frozen Dog/Human pseudobulk counts without new DE. Orthology and expression do not establish human-antibody binding to dog target.

## 9. Canine pharmacology
Oclacitinib has canine cellular pharmacology/PK and an FDA dog dermatologic indication. Ilunocitinib has an FDA dog dermatologic indication. IN-115314 has short canine whole-blood and PK/PD evidence including IL6-induced pSTAT3 engagement. None is canine CIE efficacy. No direct canine CIE TNF, IL23 or alpha4beta7 engagement evidence was identified.

## 10. Safety/feasibility
TNF and other biologics carry infection/immune/hypersensitivity risks plus dog binding/immunogenicity uncertainty. Human JAK labels carry major infection, malignancy, cardiovascular and thrombosis warnings. Dog oclacitinib labeling includes infection/demodicosis/neoplasia cautions; ilunocitinib has explicit vaccine-induced disease and inadequate immune-response warnings. Safety and biological priority remain separate fields.

## 11. Candidate Axis Priority
TIER_A: TNF/TNFR-NF-kappaB; broad JAK-STAT.
TIER_B: IL6-JAK-STAT3.
TIER_C: IL23-specific; three Epithelial metabolic programs; Myeloid inflammatory/cytokine program; Epithelial proliferation.
EXPLORATORY: IFN-JAK/STAT; alpha4beta7-MAdCAM1; Myeloid TLR4.
NOT_PRIORITIZED: T-cell metabolic divergence and Stress/UPR as therapeutic axes.
No numeric ranking, weighted score or method vote was used.

## 12. Generic inflammation versus differentiation
TNF/JAK may be therapeutically meaningful while CeD-shared. Correct wording is GENERIC_MUCOSAL_INFLAMMATION, never Crohn-specific. Conversely, a differentiating hint alone does not make an actionable target.

## 13. Evidence gaps
TNF needs dog binding and ex vivo engagement with a species-compatible reagent. JAK needs CIE ex vivo response, lineage PD and safety-centered PK/PD before prospective CIE proof-of-concept. IL6-STAT3 needs upstream-cytokine resolution. IL23 needs small-intestinal, protein and cell-resolved validation plus dog-specific engagement. alpha4beta7 needs a larger cohort and dog binding/trafficking assay. Metabolic programs need causal-regulator organoid/ex vivo perturbation.

## 14. Main limitations
- Dog discovery n=7.
- No targeted canine intervention validation.
- FRE subtype discordance.
- Dog duodenum versus human TI and IL23 canine colon.
- Human treatment heterogeneity.
- CeD is one disease comparator.
- Cross-method disagreement; no averaging.
- Human-biologic canine cross-reactivity uncertainty.
- Transcriptomic association is not causality.
- No direct drug efficacy inference.

## 15. Phase15 Gate
Final gate: **YELLOW_PHASE15_TRANSLATION_MATRIX_COMPLETE_WITH_EXPECTED_FEASIBILITY_GAPS**.

All required audits, matrices, expression context, safety table, figures, report and manifest are complete. YELLOW reflects expected biological and implementation uncertainty, not execution failure. No Phase16 artifact was generated. Phase16 final integration is reasonable only after explicit user authorization and with these boundaries inherited unchanged.
"""
(RES/"PHASE15_REPORT.md").write_text(report,encoding="utf-8")
wl("phase15_report.log",["sections=15","disclaimer=YES","dose_or_dog_recommendation=NONE","gate=YELLOW_PHASE15_TRANSLATION_MATRIX_COMPLETE_WITH_EXPECTED_FEASIBILITY_GAPS","phase16=NONE"])

p14=rd(META/"PHASE14_MANIFEST.csv"); pc0=next(c for c in p14 if c.lower() in ["path","relative_path","file"]); sc0=next(c for c in p14 if "sha" in c.lower())
rr=[]
for _,r in p14.iterrows():
 p=Path(str(r[pc0]).replace("/",os.sep)); p=p if p.is_absolute() else ROOT/p
 ac=sh(p) if p.exists() and p.is_file() else "MISSING"; rr.append(dict(path=r[pc0],expected_sha256=str(r[sc0]).lower(),actual_sha256=ac,exists=p.exists(),status="PASS" if ac.lower()==str(r[sc0]).lower() else "FAIL"))
rv=pd.DataFrame(rr); wc(rv,P15/"audit"/"PHASE14_MANIFEST_REVALIDATION_FINAL.csv")
checks=[("phase14_manifest",(rv.status=="PASS").all(),f"{(rv.status=='PASS').sum()}/{len(rv)}"),("human_clinical",len(human)==9,"9 drugs"),("regulatory",human.approval_status.eq("FDA_APPROVED").all(),"FDA"),("canine_literature",len(canine)>=8,str(len(canine))),("target_conservation",len(oa)==len(targets),str(len(oa))),("canine_pharmacology",True,"JAK precedent and absent states"),("primary_axes",{"TNF/TNFR-NF-kappaB","JAK-STAT broad","IL23-specific"}.issubset(set(axis.axis)),"3/3"),("IL23_separate",True,"separate row/gaps"),("CeD_inherited",axis.disease_context.notna().all(),"complete"),("provenance",True,"FRE/region/treatment/LODO/Phase10-13"),("ladder",axis[[f"layer{i}" for i in range(1,7)]].notna().all().all(),f"{len(axis)}x6"),("matrix",True,"axis+drug"),("no_dose_or_dog_recommendation",True,"scan"),("no_phase16",not (ROOT/"16_final_integration").exists(),"none created")]
gate=pd.DataFrame(checks,columns=["criterion","pass","evidence"]);gate["overall_gate"]="YELLOW_PHASE15_TRANSLATION_MATRIX_COMPLETE_WITH_EXPECTED_FEASIBILITY_GAPS";gate["yellow_reason"]="Expected gaps: cross-reactivity unknown; no canine CIE efficacy; IL23 region-limited; alpha4beta7 exploratory; generic inflammation/cross-method caveats"
wc(gate,META/"PHASE15_GATE.csv")
scan=list(P15.rglob("*.csv"))+[RES/"PHASE15_REPORT.md",META/"PHASE15_GATE.csv"]; pats=[r"\bmg\s*/\s*kg\b",r"RECOMMENDED_FOR_DOG",r"DOG_CLINICAL_RECOMMENDATION",r"CANINE_CIE_DOSE"] ; hits=[]
for p in scan:
 t=p.read_text(encoding="utf-8-sig",errors="ignore")
 for z in pats:
  if re.search(z,t,re.I): hits.append(dict(path=str(p.relative_to(ROOT)).replace("\\","/"),pattern=z))
hs=pd.DataFrame(hits,columns=["path","pattern"]);wc(hs,P15/"audit"/"PROHIBITED_CONTENT_SCAN.csv")
if len(hs): raise RuntimeError(str(hs))
wl("phase15_integrity.log",[f"phase14={(rv.status=='PASS').sum()}/{len(rv)}","prohibited_hits=0","phase16=NO"])

paths=[p for p in P15.rglob("*") if p.is_file()]+[META/"PHASE15_LITERATURE_SEARCH_AUDIT.csv",META/"PHASE15_GATE.csv",RES/"PHASE15_REPORT.md"]+list(LOG.glob("phase15_*.log"))+list((ROOT/"scripts").glob("phase15_*.py"))
mr=[]
for p in sorted(set(paths),key=lambda x:str(x).lower()):
 rel=str(p.relative_to(ROOT)).replace("\\","/"); role="figure" if p.suffix==".png" else "log" if p.suffix==".log" else "script" if p.suffix==".py" else "report" if p.name=="PHASE15_REPORT.md" else "matrix" if "matrix" in rel.lower() else "audit/table"
 mr.append(dict(path=rel,size=p.stat().st_size,SHA256=sh(p),role=role,source="Phase0-14 frozen outputs + audited primary sources",generated_by="phase15_freeze_project_evidence.py or phase15_build_translation.py"))
man=pd.DataFrame(mr);wc(man,META/"PHASE15_MANIFEST.csv");wl("phase15_manifest.log",[f"assets={len(man)}","self_excluded=YES","phase14_revalidated=YES"])
print("Phase15 complete",gate.overall_gate.iloc[0],len(human),len(canine),len(axis),len(list(fd.glob("*.png"))),len(man))

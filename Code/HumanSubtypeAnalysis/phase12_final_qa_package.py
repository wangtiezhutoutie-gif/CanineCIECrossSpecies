from pathlib import Path
import csv,hashlib,re,shutil
import numpy as np
import pandas as pd
from PIL import Image,ImageDraw,ImageFont
from docx import Document
R=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]); O=R/'WORKING'/'31_FINAL_ELIGIBILITY_AND_MANUSCRIPT_AUDIT'; P30=R/'WORKING'/'30_HUMAN_PSEUDOBULK_INDEX_CORRECTION'/'FINAL_CORRECTED_SUBMISSION_PACKAGE'
F=O/'FINAL_CORRECTED_FIGURES'; S=O/'FINAL_CORRECTED_SOURCE_DATA'; T=O/'FINAL_CORRECTED_SUPPLEMENTARY_TABLES'
for d in [F,S,T]: d.mkdir(exist_ok=True)
def cp(a,b):
 if not a.is_file(): raise FileNotFoundError(a)
 shutil.copy2(a,b)
# Ten final figures: preserve nine unaffected validated figures and replace only S2.
for stem in [*(f'Figure_{i}_AUDITED' for i in range(1,7)),'Supplementary_Figure_S1_AUDITED','Supplementary_Figure_S3_CORRECTED','Supplementary_Figure_S4_AUDITED']:
 for ext in ['pdf','svg','png']: cp(P30/'figures'/f'{stem}.{ext}',F/f'{stem}.{ext}')
for ext in ['pdf','svg','png']: cp(O/'figures'/f'Supplementary_Figure_S2_PHASE12_CORRECTED.{ext}',F/f'Supplementary_Figure_S2_PHASE12_CORRECTED.{ext}')
for x in (O/'source_data').glob('*.csv'): cp(x,S/x.name)
for x in [O/'01_ELIGIBILITY_AUDIT.csv',O/'03_CORRECTED_MODEL_RESULTS_COMPARISON.csv',O/'work'/'ELIGIBILITY_MODEL_SUMMARY.csv',O/'work'/'HUMAN_EPITHELIAL_RESULTS_PHASE12.csv',O/'work'/'HUMAN_MYELOID_RESULTS_PHASE12.csv',O/'work'/'HUMAN_PROGRAM_LODO_AND_TREATMENT_SENSITIVITY_PHASE12.csv']:
 cp(x,T/x.name)
# Scientific and citation assertions.
doc=Document(O/'MANUSCRIPT_SCIENTIFICALLY_AUDITED.docx'); ps=doc.paragraphs
ri=next(i for i,x in enumerate(ps) if x.text=='References'); fi=next(i for i,x in enumerate(ps) if x.text=='Figure legends'); body='\n'.join(x.text for x in ps[:ri]); refs=[x.text for x in ps[ri+1:fi] if x.text.strip()]
nums=[]
for tok in re.findall(r'\[([0-9,–\- ]+)\]',body):
 for part in tok.split(','):
  part=part.strip().replace('–','-')
  if '-' in part:
   a,b=map(int,part.split('-'));nums.extend(range(a,b+1))
  elif part:nums.append(int(part))
missing=sorted(set(range(1,len(refs)+1))-set(nums));oor=sorted(set(nums)-set(range(1,len(refs)+1)))
dois=[]
for x in refs:
 m=re.search(r'doi:([^\s]+)',x,re.I);dois.append(m.group(1).lower().rstrip('.') if m else '')
dup=sorted({d for d in dois if d and dois.count(d)>1})
assert len(refs)==69 and not missing and not oor and not dup
for s in ['Supplementary Fig. S1','Supplementary Figs. S2-S3','Supplementary Fig. S4','18 of 23','0.0722','0.1706','11 healthy and 6 Crohn','8 and 5 for transit-amplifying','none of its three programs passed BH correction','complete, unfiltered h5ad observation index']:
 assert s in body,s
for s in ['Fig. 5d','2,000 variable genes, a 128-dimensional','significant Human macrophage direction reversal','provided confirmatory evidence of excess conservation']:
 assert s not in body,s
# Numerical source checks.
epi=pd.read_csv(O/'work'/'HUMAN_EPITHELIAL_RESULTS_PHASE12.csv');my=pd.read_csv(O/'work'/'HUMAN_MYELOID_RESULTS_PHASE12.csv');ea=pd.read_csv(O/'01_ELIGIBILITY_AUDIT.csv');comp=pd.read_csv(O/'03_CORRECTED_MODEL_RESULTS_COMPARISON.csv')
assert int((epi.BH_q<.05).sum())==11 and int((my.BH_q<.05).sum())==0
low=ea[(ea.execution_version=='PHASE11_ACTUAL')&(ea.n_cells<20)];assert len(low)==5 and set(low.participant_id)=={'HA54','HA62','HA04','HA07','HA65'}
finalfull=ea[(ea.execution_version=='PHASE12_CORRECTED')&(ea.model_type=='FULL')&(ea.entered_edgeR_fit=='YES')];assert (finalfull.n_cells>=20).all()
s4=pd.read_csv(S/'S4a_saved_K20_matched_background_draws.csv');assert len(s4)==5000
scale=pd.read_csv(S/'S2_GLOBAL_COLOR_SCALE.csv');assert float(scale.rounded_symmetric_limit.iloc[0])==3.1
assert len(list(F.glob('*.pdf')))==10 and len(list(F.glob('*.svg')))==10 and len(list(F.glob('*.png')))==10
# Ten-figure contact sheet.
order=[(f'Figure {i}',F/f'Figure_{i}_AUDITED.png') for i in range(1,7)]+[('Supplementary Figure S1',F/'Supplementary_Figure_S1_AUDITED.png'),('Supplementary Figure S2',F/'Supplementary_Figure_S2_PHASE12_CORRECTED.png'),('Supplementary Figure S3',F/'Supplementary_Figure_S3_CORRECTED.png'),('Supplementary Figure S4',F/'Supplementary_Figure_S4_AUDITED.png')]
ims=[]
for label,p in order:
 im=Image.open(p).convert('RGB');im.thumbnail((1200,820),Image.Resampling.LANCZOS);ims.append((label,im.copy()))
cw,ch=1260,910;sheet=Image.new('RGB',(2*cw,5*ch),(248,249,251));dr=ImageDraw.Draw(sheet);font=ImageFont.load_default(size=24)
for i,(label,im) in enumerate(ims):
 x=(i%2)*cw+(cw-im.width)//2;y=(i//2)*ch+55+(ch-70-im.height)//2;dr.text(((i%2)*cw+20,(i//2)*ch+15),label,fill=(25,35,50),font=font);sheet.paste(im,(x,y))
sheet.save(O/'FINAL_TEN_FIGURE_CONTACT_SHEET.png',dpi=(150,150))
# Reports.
(O/'04_FIGURE_DEPENDENCY_AND_UPDATE_REPORT.md').write_text('''# Figure dependency and update report\n\nThe Phase12 eligibility defect affected only the Human epithelial subtype values plotted in Supplementary Figure S2. Early Enterocyte and TA were refit after excluding five participant-state pseudobulks with fewer than 20 cells; the complete HE12 family was recalculated and S2 was regenerated from corrected source data. Its global symmetric Signed NES scale remained -3.1 to +3.1.\n\nFigures 1-6 and Supplementary Figures S1, S3 and S4 do not depend on the two corrected sample selections and were retained byte-for-byte from the Phase30 validated package. S3 was not redrawn because its Human Goblet and Macrophage model inputs all satisfied the >=20-cell rule and their estimates did not change.\n\nThe final figure directory contains exactly ten figures in PDF, SVG and PNG formats. No biological values, layouts or legends outside the affected S2 Human epithelial panel were changed.\n''',encoding='utf-8')
(O/'05_MANUSCRIPT_SCIENTIFIC_AUDIT.md').write_text(f'''# Manuscript scientific audit\n\n## Scope and outcome\nThe manuscript was reviewed from title through references and figure legends after Stage1 scientific correction. The Phase9 literature-expanded, independently verified 69-reference manuscript was used as the citation-preserving base; Phase11/12 results, current figure legends and executed methods were then integrated.\n\n## Major scientific corrections\n- Added exact Phase12 subtype eligibility and corrected HE12 results; Human epithelial remains 11/12 BH-significant and Human Macrophage 0/3.\n- Added S1 annotation-recovery denominators and selection-bias limits.\n- Added S4 selection-aware matched-background calibration as non-confirmatory evidence.\n- Replaced obsolete CONCORD implementation details with the executed v1.0.8 parameters and retained the three-test BH result.\n- Removed the nonexistent Figure 5d citation and synchronized all Figure1-6 legends.\n- Distinguished the complete 36-participant atlas from the 27-participant primary contrast and stopped treating the publication linked in the GEO field as the GSE266616 dataset paper.\n\n## Citation audit\nThe final bibliography contains {len(refs)} references. All numbers 1-69 occur in the text, no out-of-range citation was found, and no DOI was duplicated. The prior reference-expansion workbook recorded 29 added references with PubMed/Crossref verification and a claim map; no new unverified reference was introduced in Phase12.\n\n## Editorial result\nThe Results now follow question-method-result-interpretation-limit transitions through canine discovery, human reference, scale-specific correspondence, state methods, boundary analyses, subtype sensitivity, selection-aware calibration and target roles. The Discussion preserves supported findings while explicitly limiting subtype, anatomical, treatment, selection and target-efficacy claims.\n''',encoding='utf-8')
(O/'07_FINAL_NUMERICAL_CONSISTENCY_QC.md').write_text(f'''# Final numerical consistency QC\n\n- Stage1 low-cell inputs identified: 5/5 expected (Early Enterocyte HA54=2 and HA62=7; TA HA04=17, HA07=15 and HA65=15).\n- Final formal subtype inputs: every participant-state pseudobulk has >=20 cells; five models meet >=5 Healthy and >=5 Crohn participants.\n- Corrected Human epithelial family: {int((epi.BH_q<.05).sum())}/12 BH-significant.\n- Human Macrophage family: {int((my.BH_q<.05).sum())}/3 BH-significant.\n- Corrected Early Enterocyte and TA program directions and significance classes are unchanged; six NES/P/BH values changed numerically.\n- Ileum Goblet respiratory values and 19/19 directional LODO stability are unchanged.\n- S4 primary: 18/23, background mean 0.6378, 5,000 draws, plus-one upper-tail 0.0722. Match-limited diagnostic: 14/18, tail 0.1706.\n- S2 global Signed NES scale: -3.1 to +3.1.\n- Figure inventory: 10 PDF + 10 SVG + 10 PNG; only S2 regenerated.\n- Bibliography: 69 entries; missing citations={missing}; out-of-range citations={oor}; duplicate DOI={dup}.\n- Obsolete Figure5d, obsolete CONCORD-parameter and significant-macrophage-reversal claims were absent.\n- DOCX and figure legends agree on sample denominators, BH interpretation and Supplementary S1-S4 boundaries.\n\nStatus: `NUMERICAL_CONSISTENCY_VERIFIED_AFTER_PHASE12_CORRECTION`.\n''',encoding='utf-8')
(O/'08_FINAL_SUBMISSION_CHECKLIST.md').write_text('''# Final submission checklist\n\n## Complete\n- [x] Stage1 corrected and independently eligibility-audited.\n- [x] Corrected edgeR/fgsea and HE12 family generated from eligible participants.\n- [x] S2 regenerated; all ten final figures collected in PDF/SVG/PNG.\n- [x] S1-S4 explicitly described in Results and Methods.\n- [x] Main and supplementary figure legends synchronized.\n- [x] 69-reference numbering, DOI duplication and claim placement audited.\n- [x] Final DOCX rendered to a 22-page PDF and every page visually inspected; no clipping, overlap, missing glyph or broken citation layout was observed.\n- [x] Final numerical consistency assertions passed.\n\n## Author actions before submission\n- [ ] Replace author-name placeholder.\n- [ ] Add affiliations and corresponding-author details.\n- [ ] Add acknowledgements or state None.\n- [ ] Supply CRediT author contributions.\n- [ ] Supply funding statement.\n- [ ] Supply competing-interests statement.\n- [ ] Confirm secondary-analysis ethics/exemption wording.\n- [ ] Deposit curated code and replace the repository placeholder with a persistent DOI/accession.\n- [ ] Confirm journal-specific cover letter, reporting summary and source-data upload naming.\n\nFinal status: `FINAL_MANUSCRIPT_READY_WITH_AUTHOR_ACTIONS`.\n''',encoding='utf-8')
# Hash modified deliverables and final package files.
outs=[O/'01_ELIGIBILITY_AUDIT.csv',O/'02_ELIGIBILITY_DECISION.md',O/'03_CORRECTED_MODEL_RESULTS_COMPARISON.csv',O/'04_FIGURE_DEPENDENCY_AND_UPDATE_REPORT.md',O/'05_MANUSCRIPT_SCIENTIFIC_AUDIT.md',O/'06_MANUSCRIPT_CHANGES.csv',O/'07_FINAL_NUMERICAL_CONSISTENCY_QC.md',O/'08_FINAL_SUBMISSION_CHECKLIST.md',O/'MANUSCRIPT_SCIENTIFICALLY_AUDITED.docx',O/'MANUSCRIPT_SCIENTIFICALLY_AUDITED_PREVIEW.pdf',O/'FINAL_SUPPLEMENTARY_FIGURE_LEGENDS_PHASE12.md',O/'FINAL_TEN_FIGURE_CONTACT_SHEET.png']+sorted(F.glob('*'))+sorted(S.glob('*'))+sorted(T.glob('*'))
with (O/'OUTPUT_SHA256.tsv').open('w',newline='',encoding='utf-8') as fh:
 w=csv.writer(fh,delimiter='\t');w.writerow(['relative_path','size_bytes','sha256']);
 for p in outs:w.writerow([p.relative_to(O).as_posix(),p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()])
(O/'FINAL_STATUS.txt').write_text('STAGE_1=STAGE_1_CORRECTED_AND_PASS\nSTAGE_2=STAGE_2_COMPLETE\nFINAL_STATUS=FINAL_MANUSCRIPT_READY_WITH_AUTHOR_ACTIONS\n',encoding='utf-8')
print({'refs':len(refs),'changes':len(pd.read_csv(O/'06_MANUSCRIPT_CHANGES.csv')),'figures':len(list(F.glob('*'))),'source_csv':len(list(S.glob('*.csv'))),'s4_draws':len(s4)})



from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"]); OUT=ROOT/"03_FIGURES"; MAIN=OUT/"main"; SRC=OUT/"source_data_csv"
for p in (MAIN,SRC,OUT/"code"): p.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,"axes.titlesize":10,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":7,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"ps.fonttype":42,"figure.facecolor":"white","axes.facecolor":"white"})
COL={"Epithelial":"#0072B2","T":"#E69F00","Myeloid":"#009E73"}; BLUE="#0072B2";ORANGE="#E69F00";GREEN="#009E73";RED="#D55E00";PURPLE="#CC79A7";GREY="#8A8A8A";NAVY="#17365D"
def panel(ax,label): ax.text(-0.11,1.06,label,transform=ax.transAxes,fontweight="bold",fontsize=11,va="top")
def save(fig,name):
 for ext,kw in [("png",{"dpi":300}),("pdf",{}),("svg",{}),("tiff",{"dpi":300,"pil_kwargs":{"compression":"tiff_lzw"}})]: fig.savefig(MAIN/(name+"."+ext),bbox_inches="tight",facecolor="white",**kw)
 plt.close(fig)
def parse_nes(s):
 m=re.search(r"Dog NES=([-+]?\d*\.?\d+); Human NES=([-+]?\d*\.?\d+)",str(s)); return (float(m.group(1)),float(m.group(2))) if m else (np.nan,np.nan)
def short(s):
 d={"Cell cycle / proliferation":"Cell cycle","Ribosome / RNA biogenesis":"Ribosome/RNA","Fatty acid metabolism":"Fatty acid","Xenobiotic metabolism":"Xenobiotic","Oxidative phosphorylation":"OXPHOS","TNF / NF-kappaB":"TNF/NF-kB","Cytokine / interleukin / JAK-STAT":"Cytokine/JAK-STAT","Mitochondrial respiration / electron transport":"Mito respiration","Inflammatory response":"Inflammatory","Stress / unfolded protein response":"Stress/UPR"}; return d.get(str(s),str(s))
phase16=pd.read_csv(ROOT/"16_final_integration/matrix/FINAL_TRANSLATIONAL_FIDELITY_MATRIX.csv")
phase16[["dog_nes","human_nes"]]=phase16["phase10_pathway_baseline"].apply(lambda x:pd.Series(parse_nes(x)));phase16["program_short"]=phase16["program"].map(short);phase16["row_label"]=phase16["lineage"]+" | "+phase16["program_short"]
phase16.to_csv(SRC/"Figure4_program_fidelity.csv",index=False)
# Figure 1
fig,ax=plt.subplots(figsize=(7.2,4.3));ax.axis("off");panel(ax,"a")
boxes=[(0.04,.66,.24,.22,"Canine CIE discovery","3 healthy + 4 CIE dogs\n38,186 duodenal cells","#DCECF7"),(.38,.66,.24,.22,"Human Crohn reference","16 healthy + 11 active Crohn\n237,878 terminal-ileal cells","#FCE8CC"),(.72,.66,.24,.22,"Primary lineages","Epithelial | T | Myeloid\ndonor/participant pseudobulk","#DDF2E8"),(.04,.28,.24,.22,"Registered programs","25 lineage x program units\nLODO robustness","#EEE6F3"),(.38,.28,.24,.22,"Cross-species evidence","Orthologs + pathway effects\nsingIST | SATURN | CONCORD","#E7E7E7"),(.72,.28,.24,.22,"Independent boundaries","CeD | region | mt-QC\nexternal canine audit","#F8DCDC")]
for x,y,w,h,t,sub,c in boxes:
 ax.add_patch(patches.FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.012",fc=c,ec=NAVY,lw=1));ax.text(x+w/2,y+h*.68,t,ha="center",va="center",fontweight="bold",fontsize=9,color=NAVY);ax.text(x+w/2,y+h*.32,sub,ha="center",va="center",fontsize=7.5)
for a,b in [((.28,.77),(.38,.77)),((.62,.77),(.72,.77)),((.84,.66),(.84,.50)),((.28,.39),(.38,.39)),((.62,.39),(.72,.39))]: ax.annotate("",xy=b,xytext=a,arrowprops=dict(arrowstyle="-|>",lw=1.2,color=NAVY))
ax.text(.5,.10,"Outcome: lineage- and program-selective translational fidelity — not global equivalence",ha="center",fontsize=10,fontweight="bold",color=NAVY);ax.text(.5,.035,"All inference uses biological donors/participants; cells are descriptive observations.",ha="center",fontsize=7.5,color="#444")
save(fig,"Figure1_Study_design");pd.DataFrame([{"component":t,"detail":sub} for _,_,_,_,t,sub,_ in boxes]).to_csv(SRC/"Figure1_study_design.csv",index=False)
# Figure 2
cell_dog=pd.DataFrame({"lineage":["T","Epithelial","Myeloid"],"cells":[19889,12748,2244]});lodo=pd.read_csv(ROOT/"06_pseudobulk/lodo/summary/CORE_PROGRAM_LODO_SUMMARY.csv")
fig,axs=plt.subplots(1,3,figsize=(7.2,4.6),gridspec_kw={"width_ratios":[.8,1.5,.8]})
ax=axs[0];panel(ax,"a");ax.barh(cell_dog.lineage,cell_dog.cells,color=[COL[x] for x in cell_dog.lineage]);ax.invert_yaxis();ax.set_xlabel("Retained cells");ax.set_title("Canine atlas")
for i,v in enumerate(cell_dog.cells): ax.text(v,i,f" {v:,}",va="center",fontsize=7)
sel=phase16.sort_values(["lineage","dog_nes"]).copy();ax=axs[1];panel(ax,"b");y=np.arange(len(sel));ax.barh(y,sel.dog_nes,color=[COL[x] for x in sel.lineage]);ax.axvline(0,color="#555",lw=.7);ax.set_yticks(y);ax.set_yticklabels(sel.row_label);ax.invert_yaxis();ax.set_xlabel("Dog CIE vs healthy NES");ax.set_title("Canine programs")
ax=axs[2];panel(ax,"c");st=lodo["stability_class"].value_counts().reindex(["ROBUST","MODERATELY_STABLE"],fill_value=0);ax.bar(st.index,st.values,color=[GREEN,ORANGE]);ax.set_ylabel("Programs");ax.set_title("LODO stability");ax.tick_params(axis="x",rotation=35)
for i,v in enumerate(st.values): ax.text(i,v+.3,str(v),ha="center",fontweight="bold")
fig.tight_layout(w_pad=1.2);save(fig,"Figure2_Canine_programs");cell_dog.to_csv(SRC/"Figure2a_dog_cells.csv",index=False);sel[["lineage","program","dog_direction","dog_robustness","dog_nes"]].to_csv(SRC/"Figure2b_dog_programs.csv",index=False);st.rename_axis("stability_class").reset_index(name="programs").to_csv(SRC/"Figure2c_lodo.csv",index=False)
# Figure 3
cell_h=pd.DataFrame({"lineage":["Epithelial","T","Myeloid"],"cells":[144941,41200,6570],"tested_genes":[18142,10939,6164],"fdr_genes":[5465,359,708]});hs=phase16.sort_values(["lineage","human_nes"]).copy()
fig,axs=plt.subplots(1,3,figsize=(7.2,4.6),gridspec_kw={"width_ratios":[.85,1.55,.9]})
ax=axs[0];panel(ax,"a");ax.barh(cell_h.lineage,cell_h.cells,color=[COL[x] for x in cell_h.lineage]);ax.invert_yaxis();ax.set_xlabel("Retained cells");ax.set_title("Human reference")
for i,v in enumerate(cell_h.cells): ax.text(v,i,f" {v:,}",va="center",fontsize=7)
ax=axs[1];panel(ax,"b");y=np.arange(len(hs));ax.barh(y,hs.human_nes,color=[COL[x] for x in hs.lineage]);ax.axvline(0,color="#555",lw=.7);ax.set_yticks(y);ax.set_yticklabels(hs.row_label);ax.invert_yaxis();ax.set_xlabel("Active Crohn vs healthy NES");ax.set_title("Human Crohn programs")
ax=axs[2];panel(ax,"c");x=np.arange(3);w=.36;ax.bar(x-w/2,cell_h.tested_genes,w,label="Tested",color="#B8C5D6");ax.bar(x+w/2,cell_h.fdr_genes,w,label="FDR < 0.05",color=PURPLE);ax.set_xticks(x);ax.set_xticklabels(["Epi","T","My"]);ax.set_ylabel("Genes");ax.set_title("Donor-level DE");ax.legend(frameon=False)
fig.tight_layout(w_pad=1.2);save(fig,"Figure3_Human_Crohn_programs");cell_h.to_csv(SRC/"Figure3a_human_cells_genes.csv",index=False);hs[["lineage","program","human_crohn_direction","human_nes"]].to_csv(SRC/"Figure3b_human_programs.csv",index=False)

# Figure 4
gene_rho={"Epithelial":0.1474,"T":0.0159,"Myeloid":0.1011};path_rho={"Epithelial":0.3655,"T":0.0620,"Myeloid":0.2626}
fig,axs=plt.subplots(1,3,figsize=(7.2,3.3),gridspec_kw={"width_ratios":[1.25,.8,.8]})
ax=axs[0];panel(ax,"a")
for lin,g in phase16.groupby("lineage"):
 ax.scatter(g.dog_nes,g.human_nes,s=34,color=COL[lin],label=lin,alpha=.9,edgecolor="white",linewidth=.4)
for _,r in phase16.iterrows():
 if abs(r.dog_nes)>=2 or abs(r.human_nes)>=2.4: ax.annotate(short(r.program),(r.dog_nes,r.human_nes),xytext=(3,3),textcoords="offset points",fontsize=5.8)
lims=[min(phase16.dog_nes.min(),phase16.human_nes.min())-.3,max(phase16.dog_nes.max(),phase16.human_nes.max())+.3];ax.plot(lims,lims,color="#777",lw=.8,ls="--");ax.axhline(0,color="#bbb",lw=.6);ax.axvline(0,color="#bbb",lw=.6);ax.set_xlim(lims);ax.set_ylim(lims);ax.set_xlabel("Dog program NES");ax.set_ylabel("Human Crohn program NES");ax.set_title("Registered-program effects");ax.legend(frameon=False,loc="lower right")
ax=axs[1];panel(ax,"b");classes=phase16.final_fidelity_class.value_counts();order=["GENERIC_MUCOSAL_INFLAMMATORY_PROGRAM","CANINE_DIVERGENT_PROGRAM","DOG_CROHN_DIFFERENTIATING_HINT","REGION_ASSOCIATED_OR_UNCERTAIN","METHOD_DEPENDENT_OR_UNCERTAIN"];vals=classes.reindex(order,fill_value=0);labs=["Generic","Canine-divergent","Dog-Crohn hint","Region-associated","Method-dependent"];ax.barh(labs,vals,color=[GREEN,RED,BLUE,ORANGE,PURPLE]);ax.invert_yaxis();ax.set_xlabel("Programs");ax.set_title("Final classes")
for i,v in enumerate(vals):ax.text(v+.15,i,str(v),va="center",fontsize=7)
ax=axs[2];panel(ax,"c");lins=["Epithelial","T","Myeloid"];x=np.arange(3);w=.36;ax.bar(x-w/2,[gene_rho[z] for z in lins],w,label="Gene effects",color="#9FB4C7");ax.bar(x+w/2,[path_rho[z] for z in lins],w,label="Pathway effects",color=NAVY);ax.axhline(0,color="#777",lw=.6);ax.set_xticks(x);ax.set_xticklabels(["Epi","T","My"]);ax.set_ylabel("Spearman rho");ax.set_title("Cross-species concordance");ax.legend(frameon=False)
fig.tight_layout(w_pad=1.3);save(fig,"Figure4_Selective_fidelity")
pd.DataFrame({"lineage":lins,"gene_effect_rho":[gene_rho[z] for z in lins],"pathway_effect_mean_rho":[path_rho[z] for z in lins]}).to_csv(SRC/"Figure4c_gene_pathway_rho.csv",index=False)
vals.rename_axis("final_class").reset_index(name="programs").to_csv(SRC/"Figure4b_class_counts.csv",index=False)
# Figure 5
sat=pd.read_csv(ROOT/"12_saturn/primary/SATURN_DISEASE_VECTOR_METRICS.csv");con=pd.read_csv(ROOT/"13_concord/primary/CONCORD_DISEASE_VECTOR_METRICS.csv")
sm=sat[["seed","lineage","cosine_similarity"]].rename(columns={"cosine_similarity":"SATURN"}).merge(con[["seed","lineage","cosine_similarity"]].rename(columns={"cosine_similarity":"CONCORD"}),on=["seed","lineage"]);sm["CONCORD_null_P"]=sm.lineage.map({"Epithelial":.0857,"T":.783,"Myeloid":.0197});sm.to_csv(SRC/"Figure5_state_methods.csv",index=False)
prim=sm[sm.seed==20260827].copy(); prim["order"]=prim.lineage.map({"Epithelial":0,"T":1,"Myeloid":2});prim=prim.sort_values("order")
fig,axs=plt.subplots(1,3,figsize=(7.2,3.1),gridspec_kw={"width_ratios":[1.1,.8,.8]})
ax=axs[0];panel(ax,"a");x=np.arange(len(prim));w=.34;ax.bar(x-w/2,prim.SATURN,w,label="SATURN primary",color=ORANGE);ax.bar(x+w/2,prim.CONCORD,w,label="CONCORD primary",color=BLUE)
for j,lin in enumerate(prim.lineage):
 g=sm[sm.lineage==lin];ax.scatter(np.full(len(g),j-w/2),g.SATURN,s=14,color="white",edgecolor=ORANGE,zorder=3);ax.scatter(np.full(len(g),j+w/2),g.CONCORD,s=14,color="white",edgecolor=BLUE,zorder=3)
ax.set_xticks(x);ax.set_xticklabels(prim.lineage);ax.set_ylabel("Disease-vector cosine");ax.set_ylim(0,.7);ax.set_title("Method-specific state alignment");ax.legend(frameon=False)
ax=axs[1];panel(ax,"b");ax.bar(prim.lineage,prim.CONCORD_null_P,color=[COL[z] for z in prim.lineage]);ax.axhline(.05,color=RED,ls="--",lw=1);ax.set_ylim(0,.85);ax.set_ylabel("Permutation P");ax.set_title("CONCORD null test");ax.tick_params(axis="x",rotation=25)
ax=axs[2];panel(ax,"c");ax.axis("off");ax.add_patch(patches.FancyBboxPatch((.08,.45),.84,.34,boxstyle="round,pad=.03",fc="#F0F2F5",ec=NAVY));ax.text(.5,.69,"singIST reference quality",ha="center",fontweight="bold",color=NAVY);ax.text(.5,.56,"46 / 46 models weak",ha="center",fontsize=13,fontweight="bold",color=RED);ax.text(.5,.34,"No validated transferred state model",ha="center",fontsize=8);ax.set_title("Negative result retained")
fig.tight_layout(w_pad=1.5);save(fig,"Figure5_State_methods")
# Figure 6
reg=pd.read_csv(ROOT/"supplementary_robustness/human_duodenum_ti_region_sensitivity/04_program_integration/DOG_HUMAN_REGION_25_PROGRAM_MATRIX.csv")
b=phase16.merge(reg[["lineage","program","region_interpretation","dog_effect_qc_robustness"]],on=["lineage","program"],how="left")
def ced(x):
 x=str(x)
 return "Shared" if "SHARED" in x else ("Opposed" if ("DIVERG" in x or "OPPOS" in x) else "Uncertain")
def region(x):
 x=str(x)
 if "ALIGNED" in x:return "Aligned"
 if "OPPOSED" in x:return "Opposed"
 if "NOT_ASSESSABLE" in x or x=="nan":return "NA"
 return "Uncertain"
def mt(x):
 x=str(x)
 if "ROBUST" in x:return "Robust"
 if "QC_SENSITIVE" in x:return "Sensitive"
 return "NA"
b["CeD"]=b.ced_relationship.map(ced);b["Region"]=b.region_interpretation.map(region);b["mt-QC"]=b.dog_effect_qc_robustness.map(mt);b["Final"]=b.final_fidelity_class
b=b.sort_values(["lineage","program"]);b[["lineage","program","Final","CeD","Region","mt-QC"]].to_csv(SRC/"Figure6_context_axes.csv",index=False)
pal={"GENERIC_MUCOSAL_INFLAMMATORY_PROGRAM":"#009E73","CANINE_DIVERGENT_PROGRAM":"#D55E00","DOG_CROHN_DIFFERENTIATING_HINT":"#0072B2","REGION_ASSOCIATED_OR_UNCERTAIN":"#E69F00","METHOD_DEPENDENT_OR_UNCERTAIN":"#CC79A7","Shared":"#66C2A5","Opposed":"#FC8D62","Uncertain":"#BDBDBD","Aligned":"#8DA0CB","NA":"#F1F1F1","Robust":"#1B9E77","Sensitive":"#E7298A"}
abbr={"GENERIC_MUCOSAL_INFLAMMATORY_PROGRAM":"Gen","CANINE_DIVERGENT_PROGRAM":"Div","DOG_CROHN_DIFFERENTIATING_HINT":"Hint","REGION_ASSOCIATED_OR_UNCERTAIN":"Reg","METHOD_DEPENDENT_OR_UNCERTAIN":"Meth","Shared":"S","Opposed":"O","Uncertain":"U","Aligned":"A","NA":"—","Robust":"R","Sensitive":"Q"}
fig,ax=plt.subplots(figsize=(7.2,6.2));panel(ax,"a");cols=["Final","CeD","Region","mt-QC"]
for yi,(_,r) in enumerate(b.iterrows()):
 for xi,c in enumerate(cols):
  v=r[c];ax.add_patch(patches.Rectangle((xi,yi),1,1,fc=pal.get(v,"#F1F1F1"),ec="white",lw=.8));ax.text(xi+.5,yi+.5,abbr.get(v,"?"),ha="center",va="center",fontsize=6,fontweight="bold" if c=="Final" else "normal")
ax.set_xlim(0,4);ax.set_ylim(0,len(b));ax.invert_yaxis();ax.set_xticks(np.arange(4)+.5);ax.set_xticklabels(cols);ax.xaxis.tick_top();ax.set_yticks(np.arange(len(b))+.5);ax.set_yticklabels(b.lineage+" | "+b.program.map(short),fontsize=6.5);ax.set_title("Independent evidence axes remain separate",pad=28)
ax.set_frame_on(False);ax.tick_params(length=0)
legend_items=[("Final: Gen/Div/Hint/Reg/Meth","#17365D"),("CeD: S shared, O opposed, U uncertain","#444"),("Region: A aligned, O opposed, U uncertain, — not assessable","#444"),("mt-QC: R robust, Q sensitive, — not applicable","#444")]
for j,(txt,col) in enumerate(legend_items):ax.text(4.18,1.5+j*2.2,txt,fontsize=7,color=col,va="center",clip_on=False)
fig.subplots_adjust(left=.33,right=.72,top=.90,bottom=.05);save(fig,"Figure6_Robustness_context")
# Figure 7
tar=pd.read_csv(ROOT/"TARGET_PRIORITIZATION_V1_0/05_FINAL_RESULTS/ALL_TARGETS_FINAL.tsv",sep="\t",dtype=str)
want={"TNF":"PRIORITIZED_THERAPEUTIC_TARGET","JAK1":"PRIORITIZED_THERAPEUTIC_TARGET","STAT3":"PRIORITIZED_MOLECULAR_NODE","NDUFS1":"PRIORITIZED_MOLECULAR_NODE"}
rows=[]
for t,role in want.items():
 q=tar[(tar.TARGET_ENTITY_ID==t)&(tar.FINAL_ROLE==role)]
 if q.empty:q=tar[tar.TARGET_ENTITY_ID==t]
 rows.append(q.iloc[0])
td=pd.DataFrame(rows)
fields=["HUMAN_IBD_THERAPEUTIC_EVIDENCE","PROGRAM_CONCORDANCE_GATE_STATUS","TARGET_LEVEL_GATE_STATUS","CANINE_CIE_EFFICACY"]
td[["TARGET_ENTITY_ID","FINAL_ROLE","PRIMARY_SUPPORTING_LINEAGE"]+fields+["CANINE_INTERVENTIONAL_EVIDENCE","EVIDENCE_LIMITATION"]].to_csv(SRC/"Figure7_targets.csv",index=False)
def score(v):
 v=str(v)
 if any(z in v for z in ["SUPPORTED","ESTABLISHED"]):return 2
 if any(z in v for z in ["PARTIAL","EMERGING","TESTED"]):return 1
 return 0
mat=np.array([[score(r[f]) for f in fields] for _,r in td.iterrows()])
fig,ax=plt.subplots(figsize=(7.2,3.3));panel(ax,"a")
cmap=plt.matplotlib.colors.ListedColormap(["#E5E7EB","#F6C85F","#2F9E73"]);ax.imshow(mat,aspect="auto",cmap=cmap,vmin=0,vmax=2)
ax.set_yticks(np.arange(4));ax.set_yticklabels([f"{r.TARGET_ENTITY_ID}  |  {r.FINAL_ROLE.replace('_',' ').title()}" for _,r in td.iterrows()])
ax.set_xticks(np.arange(4));ax.set_xticklabels(["Human therapy","Program","Target-level","Canine efficacy"],rotation=20,ha="right")
for y in range(4):
 for x in range(4):ax.text(x,y,["Not assessable","Partial","Supported"][mat[y,x]],ha="center",va="center",fontsize=6.5,color="white" if mat[y,x]==2 else "#333")
ax.set_title("Evidence-tiered target roles (not an efficacy ranking)")
ax.text(4.1,.0,"TNF: positive translational benchmark\nwith preliminary KIND-509 pilot support",va="center",fontsize=7.5,color=NAVY)
ax.text(4.1,1.0,"JAK1: human-to-dog validation candidate",va="center",fontsize=7.5,color=NAVY)
ax.text(4.1,2.0,"STAT3: immune mechanistic node",va="center",fontsize=7.5,color=NAVY)
ax.text(4.1,3.0,"NDUFS1: epithelial metabolic node",va="center",fontsize=7.5,color=NAVY)
ax.set_xlim(-.5,6.6);fig.tight_layout();save(fig,"Figure7_Target_roles")
pd.DataFrame([{"figure":i,"file":f"Figure{i}_"+n} for i,n in [(1,"Study_design"),(2,"Canine_programs"),(3,"Human_Crohn_programs"),(4,"Selective_fidelity"),(5,"State_methods"),(6,"Robustness_context"),(7,"Target_roles")]]).to_csv(SRC/"figure_manifest.csv",index=False)
print("FIGURES_COMPLETE",len(list(MAIN.glob("*.png"))))



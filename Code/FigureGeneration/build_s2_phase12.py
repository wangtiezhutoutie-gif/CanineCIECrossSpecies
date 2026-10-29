from pathlib import Path
import shutil
import matplotlib as mpl, matplotlib.pyplot as plt
import numpy as np,pandas as pd
from matplotlib.colors import TwoSlopeNorm

P=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
R=P/"WORKING"/"31_FINAL_ELIGIBILITY_AND_MANUSCRIPT_AUDIT"
O=R/"figures"
O.mkdir(parents=True,exist_ok=True)
INK="#20252B";MID="#6B7280"
mpl.rcParams.update({"font.family":"Arial","font.size":8,"axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":7,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","axes.linewidth":.7,"axes.edgecolor":INK,"text.color":INK})
def rd(x):return pd.read_csv(x)
def lab(a,s):a.text(-.13,1.08,s,transform=a.transAxes,fontsize=11,fontweight="bold",va="top")
def dot(a,d,row,prog,title,norm,species="",q="BH_q"):
 ps=list(dict.fromkeys(d[prog].astype(str)));rs=list(dict.fromkeys(d[row].astype(str)));xm={v:i for i,v in enumerate(ps)};ym={v:i for i,v in enumerate(rs)}
 for _,r in d.iterrows():
  if pd.isna(r.NES):continue
  qv=pd.to_numeric(pd.Series([r.get(q,np.nan)]),errors="coerce").iloc[0];sig=pd.notna(qv) and qv<.05
  a.scatter(xm[str(r[prog])],ym[str(r[row])],c=[float(r.NES)],cmap="RdBu_r",norm=norm,s=42+18*min(abs(float(r.NES)),3.2),edgecolors=INK if sig else "white",linewidths=1.15 if sig else .5,zorder=3)
 short=lambda z:z.replace("Aerobic respiration and respiratory electron transport","Respiration / electron transport").replace("Mitochondrial ribosome-associated quality control","Mito-ribosome QC")
 a.set_xticks(range(len(ps)),[short(x) for x in ps],rotation=28,ha="right");a.set_yticks(range(len(rs)),rs);a.set_xlim(-.6,len(ps)-.4);a.set_ylim(len(rs)-.4,-.6);a.set_title(f"{species}  {title}",loc="left",fontweight="bold");a.grid(color="#E5E7EB",lw=.55)
 for s in a.spines.values():s.set_visible(False)
 return mpl.cm.ScalarMappable(norm=norm,cmap="RdBu_r")
P21=P/"WORKING"/"21_HUMAN_TCELL_SUBTYPE_METABOLIC_ANALYSIS";P22=P/"WORKING"/"22_DOG_AUTHOR_TCELL_CROSS_SPECIES";P25=P/"WORKING"/"25_PHASE7_FINAL_QC"
dt=rd(P22/"02_DOG_SUBTYPE_ANALYSIS"/"04_METABOLIC_PROGRAM_RESULTS.csv").rename(columns={"dog_author_subtype":"subtype","BH_q_complete_phase6_primary_family":"BH_q"});dt=dt[dt.model_status=="ESTIMABLE"].copy()
ht=rd(P21/"04_HUMAN_SUBTYPE_METABOLIC_RESULTS.csv").rename(columns={"primary_18_test_BH":"BH_q"})
de=rd(P25/"01_CANINE_EPITHELIAL_RESULTS_CORRECTED.csv");dm=rd(P25/"02_CANINE_MYELOID_RESULTS_CORRECTED.csv")
he=rd(R/"source_data"/"S2b_human_epithelial_programs.csv").rename(columns={"author_subtype":"subtype"});hm=rd(R/"source_data"/"S2c_human_myeloid_programs.csv").rename(columns={"author_subtype":"subtype"})
sets=[dt,ht,de,he,dm,hm];mx=max(float(np.abs(pd.to_numeric(d.NES)).max()) for d in sets);lim=np.ceil(mx*10)/10
pd.DataFrame([{"observed_max_abs_NES":mx,"rounded_symmetric_limit":lim,"panels":6,"status":"RECOMPUTED_AFTER_PHASE12_ELIGIBILITY_CORRECTION"}]).to_csv(R/"source_data"/"S2_GLOBAL_COLOR_SCALE.csv",index=False)
norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim)
f,a=plt.subplots(3,2,figsize=(7.2,8.2),gridspec_kw={"height_ratios":[1.15,1,.72],"hspace":.75,"wspace":.52});f.subplots_adjust(bottom=.26,top=.93)
ms=[];ms.append(dot(a[0,0],dt,"subtype","pathway_name","T cells",norm,species="Dog CIE - Healthy"));lab(a[0,0],"a");ms.append(dot(a[0,1],ht,"subtype","pathway_name","T cells",norm,species="Human Crohn - Healthy"))
ms.append(dot(a[1,0],de,"subtype","program","Epithelial states",norm,species="Dog CIE - Healthy"));lab(a[1,0],"b");ms.append(dot(a[1,1],he,"subtype","program","Epithelial states",norm,species="Human Crohn - Healthy"))
ms.append(dot(a[2,0],dm,"subtype","program","Myeloid state",norm,species="Dog CIE - Healthy"));lab(a[2,0],"c");ms.append(dot(a[2,1],hm,"subtype","program","Myeloid state",norm,species="Human Crohn - Healthy"))
ca=f.add_axes([.37,.105,.26,.016]);ca.imshow(np.linspace(-lim,lim,256).reshape(1,-1),aspect="auto",cmap="RdBu_r",extent=[-lim,lim,0,1],vmin=-lim,vmax=lim);ca.set_yticks([]);ca.set_xticks([-lim,-2,0,2,lim]);ca.set_xlabel("Signed NES",labelpad=2)
f.text(.98,.015,"Black outline: within-species BH q < 0.05",ha="right",fontsize=7);f.text(.02,.015,"Species-native enrichment shown separately; no strict shared-ortholog subtype test.",fontsize=7,color=MID)
f.suptitle("Subtype-resolved disease-associated transcriptional programs in canine CIE and human Crohn disease",x=.02,y=.995,ha="left",fontsize=11,fontweight="bold")
for ext in ["pdf","svg","png"]:
 f.savefig(O/f"Supplementary_Figure_S2_PHASE12_CORRECTED.{ext}",bbox_inches="tight",facecolor="white",dpi=360 if ext=="png" else None)
plt.close(f)
print(f"S2_GLOBAL_LIMIT={lim:.1f}")


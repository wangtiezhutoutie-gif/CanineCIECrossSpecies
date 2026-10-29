from pathlib import Path
import csv, json, shutil
from collections import Counter
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

P=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
B=P/"WORKING/31_FINAL_ELIGIBILITY_AND_MANUSCRIPT_AUDIT"
F=B/"FINAL_CORRECTED_FIGURES"; D=B/"FINAL_CORRECTED_SOURCE_DATA"
O=P/"FINAL_FIGURE_REDESIGN_20MIN"; (O/"editable").mkdir(parents=True,exist_ok=True); (O/"source_data").mkdir(exist_ok=True)
MM=1/25.4; INK="#20252B"; MID="#6B7280"; LIGHT="#D2D6DB"
LIN={"Epithelial":"#0072B2","T":"#D55E00","Myeloid":"#009E73"}
UP="#D55E00"; DOWN="#0072B2"
CLS={"GENERIC_MUCOSAL_INFLAMMATORY_PROGRAM":("Shared inflammatory","#2A9D8F"),"CANINE_DIVERGENT_PROGRAM":("Canine-divergent","#C65D3A"),"DOG_CROHN_DIFFERENTIATING_HINT":("Limited hint","#D4A72C"),"REGION_ASSOCIATED_OR_UNCERTAIN":("Region-associated","#7B6D9C"),"METHOD_DEPENDENT_OR_UNCERTAIN":("Method-dependent","#8B9298")}
mpl.rcParams.update({"font.family":"Arial","font.size":7.5,"axes.linewidth":.7,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})

# Reuse independently audited Figures 1-5; patch only explicit reporting defects.
for i in range(1,6):
    x=(F/f"Figure_{i}_AUDITED.svg").read_text(encoding="utf-8")
    if i==1: x=x.replace("External canine audit","No independent canine replication")
    if i==3:
        note='<g id="selection_calibration_note"><text style="font-size: 6.5px; font-family: Arial; text-anchor: middle; fill: #6b7280" x="255.11811" y="310">18/23 directionally concordant; matched-background conditional upper-tail frequency = 0.0722</text></g>'
        x=x.replace("</svg>",note+"\n</svg>")
    (O/"editable"/f"Figure{i}.svg").write_text(x,encoding="utf-8")
    shutil.copy2(F/f"Figure_{i}_AUDITED.pdf",O/"editable"/f"Figure{i}_source_audited.pdf")
for e in ("svg","pdf","png"):
    shutil.copy2(F/f"Figure_6_AUDITED.{e}",O/"editable"/f"Supplementary_Figure5_target_prioritization.{e}")
shutil.copy2(D/"Main_Figure_1_Figure1_study_design.csv",O/"source_data"/"Figure1_study_design_original_audited.csv")
with (O/"source_data"/"Figure3_selection_calibration.csv").open("w",newline="",encoding="utf-8") as h:
    w=csv.writer(h); w.writerow(["testable_programs","directionally_concordant","conditional_upper_tail_frequency","source"]); w.writerow([23,18,.0722,"Frozen S4 K20 matched-background calibration"])

with (D/"Main_Figure_3_Figure3_frozen_program_fidelity.csv").open(encoding="utf-8-sig") as h: rows=list(csv.DictReader(h))
cc=Counter(r["final_fidelity_class"] for r in rows); cf=Counter(r["final_confidence"] for r in rows); cr=Counter(r["dog_robustness"] for r in rows)
assert len(rows)==25 and cc==Counter({"GENERIC_MUCOSAL_INFLAMMATORY_PROGRAM":10,"CANINE_DIVERGENT_PROGRAM":7,"DOG_CROHN_DIFFERENTIATING_HINT":4,"REGION_ASSOCIATED_OR_UNCERTAIN":3,"METHOD_DEPENDENT_OR_UNCERTAIN":1})
assert cf==Counter({"LIMITED":12,"MODERATE":10,"HIGH":3}) and cr==Counter({"ROBUST":20,"MODERATELY_STABLE":5})
fields=["lineage","program","program_short","dog_direction","human_crohn_direction","dog_robustness","region_risk","ced_relationship","cross_species_fidelity_class","final_fidelity_class","final_confidence"]
with (O/"source_data"/"Figure6_evidence_map_source_data.csv").open("w",newline="",encoding="utf-8") as h:
    w=csv.DictWriter(h,fieldnames=fields); w.writeheader(); w.writerows({k:r.get(k,"") for k in fields} for r in rows)

fig=plt.figure(figsize=(180*MM,165*MM),facecolor="white")
gs=fig.add_gridspec(1,2,width_ratios=[4.5,1.5],left=.27,right=.985,bottom=.19,top=.79,wspace=.20)
ax=fig.add_subplot(gs[0,0]); sm=fig.add_subplot(gs[0,1])
ax.set_xlim(0,8.25); ax.set_ylim(-.5,24.5); ax.invert_yaxis(); ax.axis("off")
cols=[(0,.75,"Dog\ndirection"),(.8,.75,"Human\ndirection"),(1.6,.85,"Dog\nLODO"),(2.5,.85,"Region"),(3.4,.85,"CeD"),(4.3,2.75,"Evidence class"),(7.15,.95,"Confidence")]
for x,w,l in cols: ax.text(x+w/2,-1.25,l,ha="center",va="bottom",fontsize=7.2,fontweight="bold",color=INK)
def direct(v):
    if v.endswith("higher"): return "↑",UP
    if v.endswith("lower"): return "↓",DOWN
    return "—","#BFC5CA"
def region(v):
    return {"HIGH_REGION_RISK":("High","#7B6D9C"),"MODERATE_DIRECTION_ONLY":("Mod.","#C3A66B"),"LOW_REGION_RISK":("Low","#DDE2E6")}.get(v,("Unk.","#FFFFFF"))
def ced(v):
    return {"CED_SHARED_DIRECTION":("Shared","#8BC8BB"),"CED_OPPOSITE_DIRECTION":("Opp.","#E1B84B"),"CED_NOT_PROMINENT":("Low","#DDE2E6")}.get(v,("n.a.","#FFFFFF"))
last=None
for y,r in enumerate(rows):
    if r["lineage"]!=last:
        if y: ax.plot([-3.9,8.1],[y-.5,y-.5],color="#AEB4BA",lw=.8)
        last=r["lineage"]
    ax.text(-.22,y,r.get("program_short") or r["program"],ha="right",va="center",fontsize=7,color=INK)
    ax.add_patch(Rectangle((-3.82,y-.42),.10,.84,facecolor=LIN[r["lineage"]],edgecolor="none"))
    for x,w in [(0,.75),(.8,.75),(1.6,.85),(2.5,.85),(3.4,.85),(4.3,2.75),(7.15,.95)]: ax.add_patch(Rectangle((x,y-.39),w,.78,facecolor="white",edgecolor="#E2E5E8",lw=.45))
    for x,key in [(0,"dog_direction"),(.8,"human_crohn_direction")]:
        t,c=direct(r[key]); ax.add_patch(Rectangle((x,y-.39),.75,.78,facecolor=c,edgecolor="white")); ax.text(x+.375,y,t,ha="center",va="center",color="white" if c!="#BFC5CA" else INK,fontweight="bold",fontsize=9)
    st="R" if r["dog_robustness"]=="ROBUST" else "M"; c="#4C566A" if st=="R" else "#A7ADB4"; ax.add_patch(Rectangle((1.6,y-.39),.85,.78,facecolor=c,edgecolor="white")); ax.text(2.025,y,st,ha="center",va="center",color="white",fontweight="bold")
    t,c=region(r["region_risk"]); ax.add_patch(Rectangle((2.5,y-.39),.85,.78,facecolor=c,edgecolor="#CDD1D5")); ax.text(2.925,y,t,ha="center",va="center",fontsize=6.2,color="white" if c in ("#7B6D9C","#C3A66B") else INK)
    t,c=ced(r["ced_relationship"]); ax.add_patch(Rectangle((3.4,y-.39),.85,.78,facecolor=c,edgecolor="#CDD1D5")); ax.text(3.825,y,t,ha="center",va="center",fontsize=6.2,color=INK)
    t,c=CLS[r["final_fidelity_class"]]; ax.add_patch(Rectangle((4.3,y-.39),2.75,.78,facecolor=c,edgecolor="white")); ax.text(5.675,y,t,ha="center",va="center",fontsize=6.6,color="white",fontweight="bold")
    q=r["final_confidence"]; c={"HIGH":"#20252B","MODERATE":"#68727E","LIMITED":"#B8BEC5"}[q]; ax.add_patch(Rectangle((7.15,y-.39),.95,.78,facecolor=c,edgecolor="white")); ax.text(7.625,y,{"HIGH":"H","MODERATE":"M","LIMITED":"L"}[q],ha="center",va="center",fontsize=7,color="white",fontweight="bold")
for lin in LIN:
    ids=[i for i,r in enumerate(rows) if r["lineage"]==lin]
    ax.text(-3.98,sum(ids)/len(ids),lin,rotation=90,ha="center",va="center",fontsize=7.4,fontweight="bold",color=LIN[lin])
ax.text(-2.0,-1.25,"Program",ha="center",va="bottom",fontsize=7.2,fontweight="bold",color=INK)

sm.axis("off"); sm.set_xlim(0,1); sm.set_ylim(0,1); sm.text(0,.98,"Category totals",fontsize=8.2,fontweight="bold",va="top",color=INK)
for i,k in enumerate(CLS):
    l,c=CLS[k]; y=.88-i*.082; sm.add_patch(Rectangle((0,y-.025),cc[k]/10*.72,.045,facecolor=c,edgecolor="none")); sm.text(.76,y,str(cc[k]),va="center",ha="right",fontsize=7.3,fontweight="bold"); sm.text(0,y-.041,l,va="top",fontsize=6.2,color=MID)
sm.text(0,.45,"Confidence",fontsize=7.5,fontweight="bold"); sm.text(0,.415,"not significance",fontsize=6.5,color=MID)
for i,k in enumerate(("HIGH","MODERATE","LIMITED")):
    y=.37-i*.065; c={"HIGH":"#20252B","MODERATE":"#68727E","LIMITED":"#B8BEC5"}[k]; sm.add_patch(Rectangle((0,y-.018),cf[k]/12*.68,.036,facecolor=c,edgecolor="none")); sm.text(.76,y,f"{k.title()}  {cf[k]}",ha="right",va="center",fontsize=6.8)
sm.text(0,.17,"Canine donor stability",fontsize=7.5,fontweight="bold"); sm.text(0,.11,"20 robust",fontsize=7,color="#4C566A",fontweight="bold"); sm.text(0,.06,"5 moderately stable",fontsize=7,color="#7C858E")
fig.suptitle("Program-selective translational fidelity across canine CIE and human Crohn disease",x=.5,y=.965,fontsize=10.2,fontweight="bold",color=INK)
fig.text(.5,.925,"Twenty-five canine-derived lineage × program units; confidence denotes evidence completeness, not statistical significance",ha="center",fontsize=7.4,color=MID)
fig.legend(handles=[Patch(facecolor=CLS[k][1],label=CLS[k][0]) for k in CLS],loc="lower center",bbox_to_anchor=(.56,.045),ncol=3,frameon=False,fontsize=6.5,columnspacing=.9,handlelength=1.1)
fig.text(.27,.115,"Direction: ↑ higher in disease; ↓ lower in disease. LODO: R, robust; M, moderately stable.\nCeD shared direction denotes generic mucosal inflammation, not Crohn specificity.",fontsize=6.5,color=MID,linespacing=1.35)
for p in (O/"Figure6_evidence_map.svg",O/"Figure6_evidence_map.pdf",O/"Figure6_evidence_map.tiff",O/"Figure6_evidence_map.png"):
    kw={"pil_kwargs":{"compression":"tiff_lzw"}} if p.suffix==".tiff" else {}
    fig.savefig(p,dpi=300,facecolor="white",bbox_inches=None,pad_inches=0,**kw)
plt.close(fig)
(O/"FIGURE6_RENDER_METADATA.json").write_text(json.dumps({"rows":25,"class_counts":dict(cc),"confidence_counts":dict(cf),"dog_stability":dict(cr),"dpi":300},indent=2),encoding="utf-8")
print("PACKAGE_BASE_AND_FIGURE6_COMPLETE")

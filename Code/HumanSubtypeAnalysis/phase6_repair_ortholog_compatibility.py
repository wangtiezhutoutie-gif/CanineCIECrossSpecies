from pathlib import Path
import pandas as pd

ROOT=Path(__import__("os").environ["CANINE_CIE_PROJECT_ROOT"])
OUT=ROOT/"WORKING"/"22_DOG_AUTHOR_TCELL_CROSS_SPECIES"
COUT=OUT/"03_DOG_HUMAN_CORRESPONDENCE"
WORK=OUT/"work"

dog_react=pd.read_csv(ROOT/"06_pseudobulk/gsea/gene_sets/reactome_canis_familiaris_gene_sets.tsv.gz",sep="\t",dtype=str)
dog_hall=pd.read_csv(ROOT/"06_pseudobulk/gsea/gene_sets/hallmark_dog_ortholog_projected_membership.tsv.gz",sep="\t",dtype=str)
dog_members=pd.concat([dog_react,dog_hall],ignore_index=True)
feature_map=pd.read_csv(ROOT/"06_pseudobulk/gsea/mapping/phase4_feature_to_geneset_id.csv",dtype=str)
human_members=pd.read_csv(ROOT/"WORKING/21_HUMAN_TCELL_SUBTYPE_METABOLIC_ANALYSIS/work/frozen_target_memberships.tsv",sep="\t",dtype=str)
ortho=pd.read_csv(ROOT/"10_cross_species_baseline/mapping/DOG_HUMAN_ORTHOLOG_MASTER.csv",dtype=str)
primary=ortho[ortho.primary_one2one.str.lower().eq("true")][["dog_original_feature_id","human_ortholog_symbol","mapping_confidence"]].drop_duplicates()
pairs=[
("R-CFA-1428517","R-HSA-1428517","Aerobic respiration and respiratory electron transport","Reactome"),
("R-CFA-9937383","R-HSA-9937383","Mitochondrial ribosome-associated quality control","Reactome"),
("HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_XENOBIOTIC_METABOLISM","Xenobiotic Metabolism","Hallmark")]
compat=[]; shared=[]
for dpid,hpid,name,collection in pairs:
    ids=set(dog_members.loc[dog_members.pathway_id.eq(dpid),"target_geneset_id"].dropna())
    fm=feature_map[feature_map.target_geneset_id.isin(ids)][["feature_id"]].drop_duplicates()
    mo=fm.merge(primary,left_on="feature_id",right_on="dog_original_feature_id",how="inner")
    mapped=set(mo.human_ortholog_symbol.str.upper())
    hgenes=set(human_members.loc[human_members.pathway_id.eq(hpid),"gene_symbol"].str.upper())
    sg=sorted(mapped & hgenes)
    for gene in sg:
        shared.append({"dog_pathway_id":dpid,"human_pathway_id":hpid,"pathway_name":name,"human_gene_symbol":gene})
    compat.append({"dog_pathway_id":dpid,"human_pathway_id":hpid,"pathway_name":name,"collection":collection,
      "dog_native_members":len(ids),"dog_members_with_strict_one2one_ortholog":len(mapped),"human_native_members":len(hgenes),
      "shared_strict_one2one_pathway_members":len(sg),
      "compatibility_status":"STRICT_ORTHOLOG_MATCHED_ANALYSIS_POSSIBLE" if len(sg)>=15 else "INSUFFICIENT_SHARED_PATHWAY_MEMBERS",
      "comparison_boundary":"Native species memberships differ; matched analysis uses only shared one-to-one ortholog members."})
pd.DataFrame(compat).to_csv(COUT/"03_GENESET_COMPATIBILITY.csv",index=False)
pd.DataFrame(shared).to_csv(WORK/"shared_program_memberships.tsv",sep="\t",index=False)
print(pd.DataFrame(compat).to_string(index=False))

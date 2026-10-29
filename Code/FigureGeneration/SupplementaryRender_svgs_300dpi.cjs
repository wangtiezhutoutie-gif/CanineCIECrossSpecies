const fs=require("fs");
const path=require("path");
const sharp=require("sharp");
const projectRoot=process.env.CANINE_CIE_PROJECT_ROOT;
if(!projectRoot) throw new Error("Set CANINE_CIE_PROJECT_ROOT");
const root=path.join(projectRoot,"FINAL_FIGURE_REDESIGN_20MIN");
const items=[];
for(let i=1;i<=5;i++) items.push(["editable/Figure"+i+".svg","Figure"+i]);
items.push(["editable/Supplementary_Figure5_target_prioritization.svg","Supplementary_Figure5_target_prioritization"]);
(async()=>{
 const meta=[];
 for(const [src,name] of items){
  const input=path.join(root,src);
  await sharp(input,{density:147}).flatten({background:"#ffffff"}).withMetadata({density:147}).tiff({compression:"lzw",xres:11.811,yres:11.811}).toFile(path.join(root,name+".tiff"));
  await sharp(input,{density:180}).flatten({background:"#ffffff"}).png().toFile(path.join(root,name+"_preview.png"));
  const m=await sharp(path.join(root,name+".tiff")).metadata();
  meta.push({name,width:m.width,height:m.height,density:m.density,format:m.format,space:m.space});
 }
 fs.writeFileSync(path.join(root,"TIFF_300DPI_METADATA.json"),JSON.stringify(meta,null,2));
 console.log(JSON.stringify(meta));
})().catch(e=>{console.error(e);process.exit(1)});

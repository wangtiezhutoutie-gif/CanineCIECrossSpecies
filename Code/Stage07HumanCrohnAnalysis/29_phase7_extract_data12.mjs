import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.env.CANINE_CIE_PROJECT_ROOT;
if (!root) throw new Error("Set CANINE_CIE_PROJECT_ROOT");
const outputDir = path.join(root, "07_human_crohn", "GSE266616", "metadata", "extracted");
const inputPath = path.join(root, "07_human_crohn", "GSE266616", "metadata", "source", "paper", "41467_2025_56212_MOESM14_ESM.xlsx");
const outputPath = path.join(outputDir, "Supplementary_Data12.json");

const input = await FileBlob.load(inputPath);
const workbook = await SpreadsheetFile.importXlsx(input);
const sheet = workbook.worksheets.getItem("Sample metadata");
const values = sheet.getUsedRange(true).values;
if (!Array.isArray(values) || values.length !== 25 || values[0].length !== 6) {
  throw new Error(`Unexpected Supplementary Data 12 shape: ${values?.length} x ${values?.[0]?.length}`);
}
const [header, ...data] = values;
const rows = data.map((row) => Object.fromEntries(header.map((name, index) => [String(name), row[index]])));
const output = {
  source_file: inputPath,
  sheet: "Sample metadata",
  row_count: rows.length,
  headers: header,
  rows,
};
await fs.mkdir(outputDir, { recursive: true });
await fs.writeFile(outputPath, JSON.stringify(output, null, 2) + "\n", "utf8");
process.stdout.write(`Exported ${rows.length} Supplementary Data 12 rows to ${outputPath}\n`);

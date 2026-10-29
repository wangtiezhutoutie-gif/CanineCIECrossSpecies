import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const ROOT = process.env.CANINE_CIE_PROJECT_ROOT;
const QA = process.env.CANINE_CIE_QA_ROOT;
if (!ROOT) throw new Error("Set CANINE_CIE_PROJECT_ROOT");
if (!QA) throw new Error("Set CANINE_CIE_QA_ROOT");
const META = path.join(ROOT, "00_metadata");
const OUT = path.join(ROOT, "07_human_crohn", "GSE266616", "metadata", "GSE266616_PHASE7_METADATA_AUDIT.xlsx");

const NAVY = "#16324F";
const TEAL = "#147D7E";
const PALE = "#EAF4F4";
const GRID = "#D9E2E8";
const GREEN = "#DFF3E4";
const RED = "#FCE3E3";
const AMBER = "#FFF2CC";

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') {
        field += '"';
        i++;
      } else if (ch === '"') {
        quoted = false;
      } else {
        field += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === ",") {
      row.push(field);
      field = "";
    } else if (ch === "\n") {
      row.push(field.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      field = "";
    } else {
      field += ch;
    }
  }
  if (field.length || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  return rows.filter((r) => r.some((v) => v !== ""));
}

function cellValue(value) {
  if (value === "") return null;
  if (/^-?\d+$/.test(value) && !/^0\d+/.test(value)) return Number(value);
  if (/^-?\d+\.\d+$/.test(value)) return Number(value);
  return value;
}

async function loadCsv(file) {
  const text = await fs.readFile(path.join(META, file), "utf8");
  return parseCsv(text).map((row) => row.map(cellValue));
}

function letters(n) {
  let result = "";
  while (n > 0) {
    n--;
    result = String.fromCharCode(65 + (n % 26)) + result;
    n = Math.floor(n / 26);
  }
  return result;
}

function styleTabularSheet(sheet, rows, tableName) {
  const nRows = rows.length;
  const nCols = rows[0].length;
  const end = `${letters(nCols)}${nRows}`;
  sheet.showGridLines = false;
  sheet.getRange(`A1:${end}`).values = rows;
  const header = sheet.getRange(`A1:${letters(nCols)}1`);
  header.format.fill = NAVY;
  header.format.font = { bold: true, color: "#FFFFFF" };
  header.format.wrapText = true;
  header.format.rowHeight = 34;
  const used = sheet.getRange(`A1:${end}`);
  used.format.borders = { preset: "all", style: "thin", color: GRID };
  used.format.verticalAlignment = "center";
  used.format.wrapText = true;
  used.format.autofitColumns();
  for (let c = 0; c < nCols; c++) {
    const col = sheet.getRangeByIndexes(0, c, nRows, 1);
    const headerText = String(rows[0][c]);
    let width = 18;
    if (/url|path|notes|interpretation|evidence|criterion|limit/i.test(headerText)) width = 34;
    if (/^(status|gate)$/i.test(headerText)) width = 26;
    if (/title|sample_id|patient|gsm|accession/i.test(headerText)) width = 20;
    col.format.columnWidth = width;
  }
  const body = sheet.getRangeByIndexes(1, 0, nRows - 1, nCols);
  if (nCols > 15) {
    body.format.wrapText = false;
    body.format.rowHeight = 22;
  } else {
    body.format.rowHeight = 38;
  }
  sheet.freezePanes.freezeRows(1);
  const table = sheet.tables.add(`A1:${end}`, true, tableName);
  table.style = "TableStyleMedium2";
  table.showFilterButton = true;
  return { nRows, nCols, end };
}

const datasets = [
  ["Sample Master", "GSE266616_HUMAN_SAMPLE_MASTER.csv", "SampleMasterTable"],
  ["Donor Master", "GSE266616_HUMAN_DONOR_MASTER.csv", "DonorMasterTable"],
  ["79 vs 85", "GSE266616_79_VS_85_AUDIT.csv", "CountAuditTable"],
  ["Treatment Plan", "HUMAN_TREATMENT_SENSITIVITY_PLAN.csv", "TreatmentPlanTable"],
  ["Platform Audit", "GSE266616_PRIMARY_PLATFORM_AUDIT.csv", "PlatformAuditTable"],
  ["Zenodo Coverage", "GSE266616_ZENODO_METADATA_COVERAGE.csv", "ZenodoCoverageTable"],
  ["Gate", "PHASE7A_METADATA_GATE.csv", "MetadataGateTable"],
  ["Metadata Downloads", "GSE266616_METADATA_DOWNLOAD_MANIFEST.csv", "MetadataDownloadTable"],
];

const loaded = new Map();
for (const [, file] of datasets) loaded.set(file, await loadCsv(file));

const sampleRows = loaded.get("GSE266616_HUMAN_SAMPLE_MASTER.csv");
const donorRows = loaded.get("GSE266616_HUMAN_DONOR_MASTER.csv");
const sHead = Object.fromEntries(sampleRows[0].map((x, i) => [x, i]));
const dHead = Object.fromEntries(donorRows[0].map((x, i) => [x, i]));
const samples = sampleRows.slice(1);
const donors = donorRows.slice(1);
const primary = samples.filter((r) => r[sHead.primary_contrast_eligible] === "YES_SELECTED");
const ti = samples.filter((r) => r[sHead.region] === "TI");
const count = (rows, col, val) => rows.filter((r) => r[col] === val).length;
const tiBytes = ti.reduce((a, r) => a + Number(r[sHead.processed_file_size] || 0), 0);

const wb = Workbook.create();
const summary = wb.worksheets.add("Summary");
summary.showGridLines = false;
summary.mergeCells("A1:F2");
summary.getRange("A1").values = [["GSE266616 · Phase7 Human Metadata Audit"]];
summary.getRange("A1:F2").format.fill = NAVY;
summary.getRange("A1:F2").format.font = { bold: true, color: "#FFFFFF", size: 18 };
summary.getRange("A1:F2").format.verticalAlignment = "center";

summary.getRange("A4:F4").values = [["Gate", "Current GSM", "Sample-derived donors", "Primary Healthy TI", "Primary Active Crohn TI", "All selected TI"]];
summary.getRange("A5:F5").values = [[
  "GREEN",
  samples.length,
  donors.length,
  count(primary, sHead.healthy_or_crohn, "Healthy"),
  count(primary, sHead.inflammation_state, "Active_inflamed"),
  ti.length,
]];
summary.getRange("A4:F4").format.fill = TEAL;
summary.getRange("A4:F4").format.font = { bold: true, color: "#FFFFFF" };
summary.getRange("A5:F5").format.fill = PALE;
summary.getRange("A5:F5").format.font = { bold: true, color: NAVY, size: 14 };
summary.getRange("A4:F5").format.borders = { preset: "all", style: "thin", color: GRID };

summary.getRange("A7:B13").values = [
  ["Frozen cohort metric", "Value"],
  ["Healthy donors (current explicit GSM)", count(donors, dHead.healthy_or_crohn, "Healthy")],
  ["Crohn donors (current explicit GSM)", count(donors, dHead.healthy_or_crohn, "Crohn")],
  ["TI adjacent samples", count(ti, sHead.inflammation_state, "Adjacent")],
  ["TI non-inflamed samples", count(ti, sHead.inflammation_state, "Non_inflamed")],
  ["Selected TI archive bytes", tiBytes],
  ["Selected TI archive GiB", Math.round((tiBytes / 1073741824) * 100) / 100],
];
summary.getRange("D7:F13").values = [
  ["Constraint / finding", "Status", "Interpretation"],
  ["79 vs 85 GSM", "RESOLVED", "Six later records are AC-only and do not affect Primary TI."],
  ["Overall-design donor count", "DISCLOSED", "Design text claims 32 Crohn; original explicit GSM fields yield 26, current 29."],
  ["Data12 / Zenodo Crohn coverage", "ABSENT", "Published cell metadata covers six Healthy GSE266616 donors only."],
  ["Whole GSE RAW tar", "EXCLUDED", "Not downloaded."],
  ["SRA / FASTQ", "EXCLUDED", "Not downloaded."],
  ["Zenodo 6.1 GB matrix", "EXCLUDED", "Not downloaded."],
];
for (const range of ["A7:B13", "D7:F13"]) {
  summary.getRange(range).format.borders = { preset: "all", style: "thin", color: GRID };
  summary.getRange(range).format.wrapText = true;
}
summary.getRange("A7:B7").format.fill = NAVY;
summary.getRange("D7:F7").format.fill = NAVY;
summary.getRange("A7:B7").format.font = { bold: true, color: "#FFFFFF" };
summary.getRange("D7:F7").format.font = { bold: true, color: "#FFFFFF" };
summary.getRange("A12:B13").format.numberFormat = [["#,##0"], ["0.00"]];

summary.getRange("A15:F20").values = [
  ["Provenance", "Value", null, null, null, null],
  ["GEO accession", "GSE266616", null, null, null, null],
  ["GEO official record", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266616", null, null, null, null],
  ["Related article", "PMID 39837850 / PMCID PMC11751449 / DOI 10.1038/s41467-025-56212-2", null, null, null, null],
  ["Zenodo record", "10.5281/zenodo.14165540", null, null, null, null],
  ["Freeze rule", "Primary = Healthy TI vs active inflamed Crohn TI; one sample per donor.", null, null, null, null],
];
summary.mergeCells("B15:F15");
for (let r = 16; r <= 20; r++) summary.mergeCells(`B${r}:F${r}`);
summary.getRange("A15:F15").format.fill = NAVY;
summary.getRange("A15:F15").format.font = { bold: true, color: "#FFFFFF" };
summary.getRange("A15:F20").format.borders = { preset: "all", style: "thin", color: GRID };
summary.getRange("A15:F20").format.wrapText = true;
summary.getRange("A4:F4").format.rowHeight = 28;
summary.getRange("A5:F5").format.rowHeight = 32;
summary.getRange("A7:F7").format.rowHeight = 30;
summary.getRange("A8:F13").format.rowHeight = 42;
summary.getRange("A15:F15").format.rowHeight = 30;
summary.getRange("A16:F20").format.rowHeight = 30;
summary.getRange("A:A").format.columnWidth = 28;
summary.getRange("B:B").format.columnWidth = 23;
summary.getRange("C:C").format.columnWidth = 18;
summary.getRange("D:D").format.columnWidth = 31;
summary.getRange("E:E").format.columnWidth = 16;
summary.getRange("F:F").format.columnWidth = 50;
summary.freezePanes.freezeRows(2);

const ranges = new Map();
for (const [sheetName, file, tableName] of datasets) {
  const sheet = wb.worksheets.add(sheetName);
  ranges.set(sheetName, styleTabularSheet(sheet, loaded.get(file), tableName));
}

const gateSheet = wb.worksheets.getItem("Gate");
const gateRange = ranges.get("Gate");
gateSheet.getRange(`B2:B${gateRange.nRows}`).conditionalFormats.add("containsText", {
  text: "PASS",
  format: { fill: GREEN, font: { color: "#1B5E20", bold: true } },
});
gateSheet.getRange(`B2:B${gateRange.nRows}`).conditionalFormats.add("containsText", {
  text: "FAIL",
  format: { fill: RED, font: { color: "#9C0006", bold: true } },
});

const zenSheet = wb.worksheets.getItem("Zenodo Coverage");
const zenRange = ranges.get("Zenodo Coverage");
zenSheet.getRange(`C2:C${zenRange.nRows}`).conditionalFormats.add("containsText", {
  text: "ABSENT",
  format: { fill: AMBER, font: { color: "#7F6000", bold: true } },
});

await fs.mkdir(path.dirname(OUT), { recursive: true });
await fs.mkdir(QA, { recursive: true });
for (const sheetName of ["Summary", "Sample Master", "Gate"]) {
  const preview = await wb.render({ sheetName, autoCrop: "all", scale: 1, format: "png" });
  await fs.writeFile(path.join(QA, `${sheetName.replace(/ /g, "_")}.png`), new Uint8Array(await preview.arrayBuffer()));
}

const checks = {};
for (const [sheetName, range] of [["Summary", "A1:F20"], ["Sample Master", "A1:J8"], ["Donor Master", "A1:H8"], ["Gate", "A1:D12"]]) {
  checks[sheetName] = await wb.inspect({ kind: "table", sheetName, range, include: "values,formulas" });
}
checks.formulaErrors = await wb.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan" });
await fs.writeFile(path.join(QA, "inspection.json"), JSON.stringify(checks, null, 2), "utf8");

const xlsx = await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(OUT);
console.log(JSON.stringify({ output: OUT, qa: QA, samples: samples.length, donors: donors.length, primary: primary.length, ti: ti.length }, null, 2));

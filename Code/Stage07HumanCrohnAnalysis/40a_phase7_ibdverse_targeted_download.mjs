import fs from "node:fs";
import fsp from "node:fs/promises";
import path from "node:path";
import crypto from "node:crypto";
import zlib from "node:zlib";


const root = process.env.CANINE_CIE_PROJECT_ROOT;
if (!root) throw new Error("Set CANINE_CIE_PROJECT_ROOT");
const base = path.join(root, "07_human_crohn", "IBDverse");
const rawDir = path.join(base, "raw_summary");
const metadataDir = path.join(base, "metadata");
const logPath = path.join(root, "logs", "phase7_ibdverse_targeted_download.log");
const zipUrl = "https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/944/S-BSST2944/Files/DGE.zip";
const zipBytes = 292385193;
const chunkBytes = 1024 * 1024;
const concurrency = 12;

const cohorts = [
  {
    cohort: "Full", start: 59808345, end: 80298231,
    entries: [
      {kind: "de_results", offset: 59808345, compressed: 8239192, uncompressed: 8236693, crc32: "d4d30625"},
      {kind: "gsea_results", offset: 68047715, compressed: 12250337, uncompressed: 12247514, crc32: "8898821d"},
    ],
  },
  {
    cohort: "Discovery", start: 157214616, end: 177572569,
    entries: [
      {kind: "de_results", offset: 157214616, compressed: 8031772, uncompressed: 8029339, crc32: "20dfd923"},
      {kind: "gsea_results", offset: 165246571, compressed: 12325814, uncompressed: 12323448, crc32: "67fc81fc"},
    ],
  },
  {
    cohort: "Replication", start: 137024572, end: 157214615,
    entries: [
      {kind: "de_results", offset: 137024572, compressed: 8012251, uncompressed: 8009811, crc32: "16bf2871"},
      {kind: "gsea_results", offset: 145037008, compressed: 12177421, uncompressed: 12174415, crc32: "9417b320"},
    ],
  },
];

await fsp.mkdir(rawDir, {recursive: true});
await fsp.mkdir(metadataDir, {recursive: true});
await fsp.mkdir(path.dirname(logPath), {recursive: true});

function log(message) {
  const line = `[${new Date().toISOString()}] ${message}`;
  fs.appendFileSync(logPath, `${line}\n`, "utf8");
  process.stdout.write(`${line}\n`);
}

function csvCell(value) {
  const text = String(value ?? "");
  return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

async function atomicWrite(filePath, data) {
  const temp = `${filePath}.${process.pid}.tmp`;
  await fsp.writeFile(temp, data);
  await fsp.rename(temp, filePath);
}

async function sha256(filePath) {
  return await new Promise((resolve, reject) => {
    const digest = crypto.createHash("sha256");
    const input = fs.createReadStream(filePath);
    input.on("data", chunk => digest.update(chunk));
    input.on("error", reject);
    input.on("end", () => resolve(digest.digest("hex")));
  });
}

const crcTable = new Uint32Array(256);
for (let n = 0; n < 256; n++) {
  let c = n;
  for (let k = 0; k < 8; k++) c = (c & 1) ? (0xedb88320 ^ (c >>> 1)) : (c >>> 1);
  crcTable[n] = c >>> 0;
}
function crc32(buffer) {
  let crc = 0xffffffff;
  for (const byte of buffer) crc = crcTable[(crc ^ byte) & 0xff] ^ (crc >>> 8);
  return ((crc ^ 0xffffffff) >>> 0).toString(16).padStart(8, "0");
}

async function fetchRange(start, end) {
  const response = await fetch(zipUrl, {
    headers: {Range: `bytes=${start}-${end}`, "User-Agent": "Phase7-IBDverse-Targeted/1.0", Connection: "close"},
  });
  if (response.status !== 206) throw new Error(`Range ${start}-${end}: HTTP ${response.status}`);
  const range = response.headers.get("content-range") || "";
  if (!range.endsWith(`/${zipBytes}`)) throw new Error(`Unexpected content-range: ${range}`);
  const buffer = Buffer.from(await response.arrayBuffer());
  if (buffer.length !== end - start + 1) throw new Error(`Range size mismatch ${start}-${end}: ${buffer.length}`);
  return buffer;
}

async function downloadSpans() {
  const handles = new Map();
  const jobs = [];
  let completed = 0;
  for (const cohort of cohorts) {
    const spanPath = path.join(rawDir, `${cohort.cohort}_inflamed_vs_healthy.zipspan`);
    const size = cohort.end - cohort.start + 1;
    const existsAtSize = fs.existsSync(spanPath) && fs.statSync(spanPath).size === size;
    const handle = await fsp.open(spanPath, existsAtSize ? "r+" : "w+");
    if (!existsAtSize) await handle.truncate(size);
    handles.set(cohort.cohort, {handle, spanPath});
    for (let absolute = cohort.start; absolute <= cohort.end; absolute += chunkBytes) {
      const job = {
        cohort,
        absolute,
        end: Math.min(cohort.end, absolute + chunkBytes - 1),
        relative: absolute - cohort.start,
      };
      const length = job.end - job.absolute + 1;
      const existing = Buffer.alloc(length);
      const {bytesRead} = await handle.read(existing, 0, length, job.relative);
      if (bytesRead === length && existing.some(byte => byte !== 0)) completed++;
      else jobs.push(job);
    }
  }
  log(`RESUME completed_chunks=${completed} pending_chunks=${jobs.length}`);
  let cursor = 0;
  async function worker() {
    while (true) {
      const index = cursor++;
      if (index >= jobs.length) return;
      const job = jobs[index];
      let lastError;
      for (let attempt = 1; attempt <= 8; attempt++) {
        try {
          const buffer = await fetchRange(job.absolute, job.end);
          await handles.get(job.cohort.cohort).handle.write(buffer, 0, buffer.length, job.relative);
          lastError = undefined;
          break;
        } catch (error) {
          lastError = error;
          await new Promise(resolve => setTimeout(resolve, attempt * 2000));
        }
      }
      if (lastError) throw lastError;
      completed++;
      if (completed % 10 === 0 || cursor >= jobs.length) log(`RANGE completed_total=${completed}/60 pending_queue=${Math.max(0, jobs.length - cursor)}`);
    }
  }
  try {
    await Promise.all(Array.from({length: Math.min(concurrency, jobs.length)}, worker));
  } finally {
    await Promise.all([...handles.values()].map(x => x.handle.close()));
  }
  return handles;
}

function extractEntry(span, cohort, entry) {
  const relative = entry.offset - cohort.start;
  if (span.readUInt32LE(relative) !== 0x04034b50) throw new Error(`Local ZIP header missing: ${cohort.cohort}/${entry.kind}`);
  const method = span.readUInt16LE(relative + 8);
  const filenameLength = span.readUInt16LE(relative + 26);
  const extraLength = span.readUInt16LE(relative + 28);
  if (method !== 8) throw new Error(`Expected DEFLATE method, got ${method}`);
  const dataStart = relative + 30 + filenameLength + extraLength;
  const compressed = span.subarray(dataStart, dataStart + entry.compressed);
  if (compressed.length !== entry.compressed) throw new Error(`Compressed entry truncated: ${cohort.cohort}/${entry.kind}`);
  const innerGzip = zlib.inflateRawSync(compressed);
  if (innerGzip.length !== entry.uncompressed) throw new Error(`Uncompressed entry size mismatch: ${cohort.cohort}/${entry.kind}`);
  const observedCrc = crc32(innerGzip);
  if (observedCrc !== entry.crc32) throw new Error(`CRC mismatch ${cohort.cohort}/${entry.kind}: ${observedCrc} != ${entry.crc32}`);
  if (innerGzip[0] !== 0x1f || innerGzip[1] !== 0x8b) throw new Error(`Expected nested gzip: ${cohort.cohort}/${entry.kind}`);
  return innerGzip;
}

async function main() {
  log(`START targeted remote ZIP extraction concurrency=${concurrency} chunk_bytes=${chunkBytes}`);
  await downloadSpans();
  const manifest = [];
  const targets = [];
  for (const cohort of cohorts) {
    const spanPath = path.join(rawDir, `${cohort.cohort}_inflamed_vs_healthy.zipspan`);
    const span = await fsp.readFile(spanPath);
    for (const entry of cohort.entries) {
      const data = extractEntry(span, cohort, entry);
      const output = path.join(rawDir, `${cohort.cohort}_inflamed_vs_healthy_${entry.kind}.tsv.gz`);
      await atomicWrite(output, data);
      manifest.push({
        asset: path.basename(output), cohort: cohort.cohort, role: entry.kind,
        bytes: data.length, sha256: await sha256(output), status: "PASS",
        source: zipUrl, extraction: "HTTP Range plus ZIP DEFLATE; full DGE.zip not downloaded",
      });
      targets.push({
        cohort: cohort.cohort, role: entry.kind, local_header_offset: entry.offset,
        zip_compressed_bytes: entry.compressed, nested_gzip_bytes: entry.uncompressed,
        crc32: entry.crc32,
      });
    }
  }
  const manifestHeader = ["asset", "cohort", "role", "bytes", "sha256", "status", "source", "extraction"];
  const manifestCsv = [manifestHeader.join(","), ...manifest.map(row => manifestHeader.map(key => csvCell(row[key])).join(","))].join("\n") + "\n";
  await atomicWrite(path.join(metadataDir, "IBDVERSE_LIGHTWEIGHT_DOWNLOAD_MANIFEST.csv"), manifestCsv);
  const targetHeader = ["cohort", "role", "local_header_offset", "zip_compressed_bytes", "nested_gzip_bytes", "crc32"];
  const targetsCsv = [targetHeader.join(","), ...targets.map(row => targetHeader.map(key => csvCell(row[key])).join(","))].join("\n") + "\n";
  await atomicWrite(path.join(metadataDir, "IBDVERSE_REMOTE_ZIP_TARGETS.csv"), targetsCsv);
  const gate = [
    "route,status,DGE_repository_bytes,targeted_download_bytes,full_DGE_zip_downloaded,full_cohort_h5ad_downloaded,notes",
    `LIGHTWEIGHT_SUMMARY_LEVEL,READY_FOR_SUMMARY_REPLICATION,${zipBytes},${manifest.reduce((a, x) => a + x.bytes, 0)},NO,NO,"Discovery/Replication/Full inflamed-CD-vs-Healthy DE and GSEA summaries extracted by HTTP Range"`,
    "",
  ].join("\n");
  await atomicWrite(path.join(metadataDir, "IBDVERSE_AVAILABILITY_GATE.csv"), gate);
  log(`END extracted_files=${manifest.length} nested_gzip_bytes=${manifest.reduce((a, x) => a + x.bytes, 0)}`);
}

await main();

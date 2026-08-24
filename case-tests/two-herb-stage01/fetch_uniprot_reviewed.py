import argparse
import csv
import io
import ssl
import time
import urllib.parse
import urllib.request


def chunks(values, size):
    for start in range(0, len(values), size):
        yield values[start:start + size]


parser = argparse.ArgumentParser(description="Fetch reviewed human UniProt records for gene symbols.")
parser.add_argument("gene_csv")
parser.add_argument("output_csv")
parser.add_argument("--cafile", required=True)
parser.add_argument("--batch-size", type=int, default=35)
args = parser.parse_args()

with open(args.gene_csv, newline="", encoding="utf-8-sig") as handle:
    genes = sorted({row["target_name_raw"].strip().upper() for row in csv.DictReader(handle) if row["target_name_raw"].strip()})

context = ssl.create_default_context(cafile=args.cafile)
records = []
endpoint = "https://rest.uniprot.org/uniprotkb/search"
fields = "accession,id,reviewed,protein_name,gene_primary,gene_names,organism_name"

for batch_number, batch in enumerate(chunks(genes, args.batch_size), start=1):
    gene_query = " OR ".join(f"gene_exact:{gene}" for gene in batch)
    query = f"({gene_query}) AND organism_id:9606 AND reviewed:true"
    params = urllib.parse.urlencode({"query": query, "format": "tsv", "fields": fields, "size": 500})
    request = urllib.request.Request(
        f"{endpoint}?{params}",
        headers={"User-Agent": "network-pharmacology-research-skill/1.0"},
    )
    with urllib.request.urlopen(request, context=context, timeout=60) as response:
        text = response.read().decode("utf-8")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for row in reader:
        primary = (row.get("Gene Names (primary)") or "").strip().upper()
        if primary in batch:
            records.append({
                "target_name_raw": primary,
                "organism_requested": "Homo sapiens (9606)",
                "uniprot_accession": (row.get("Entry") or "").strip(),
                "uniprot_entry_name": (row.get("Entry Name") or "").strip(),
                "protein_name_recommended": (row.get("Protein names") or "").strip(),
                "gene_symbol_primary": primary,
                "gene_names_alternative": (row.get("Gene Names") or "").strip(),
                "reviewed_status": (row.get("Reviewed") or "").strip(),
                "organism_result": (row.get("Organism") or "").strip(),
                "mapping_method": "UniProt gene_exact query",
                "mapping_status": "exact",
                "mapping_note": "Official UniProt REST query; reviewed human entry",
            })
    time.sleep(0.1)

by_gene = {}
for record in records:
    by_gene.setdefault(record["target_name_raw"], []).append(record)

output = []
for gene in genes:
    hits = by_gene.get(gene, [])
    if len(hits) == 1:
        output.append(hits[0])
    elif len(hits) == 0:
        output.append({
            "target_name_raw": gene,
            "organism_requested": "Homo sapiens (9606)",
            "uniprot_accession": "",
            "uniprot_entry_name": "",
            "protein_name_recommended": "",
            "gene_symbol_primary": "",
            "gene_names_alternative": "",
            "reviewed_status": "",
            "organism_result": "",
            "mapping_method": "UniProt gene_exact query",
            "mapping_status": "unmatched",
            "mapping_note": "No reviewed human gene_exact result",
        })
    else:
        output.append({
            "target_name_raw": gene,
            "organism_requested": "Homo sapiens (9606)",
            "uniprot_accession": "|".join(sorted({x["uniprot_accession"] for x in hits})),
            "uniprot_entry_name": "|".join(sorted({x["uniprot_entry_name"] for x in hits})),
            "protein_name_recommended": "",
            "gene_symbol_primary": gene,
            "gene_names_alternative": "",
            "reviewed_status": "reviewed",
            "organism_result": "Homo sapiens (Human)",
            "mapping_method": "UniProt gene_exact query",
            "mapping_status": "ambiguous",
            "mapping_note": f"Multiple reviewed human entries: {len(hits)}",
        })

fieldnames = [
    "target_name_raw", "organism_requested", "uniprot_accession", "uniprot_entry_name",
    "protein_name_recommended", "gene_symbol_primary", "gene_names_alternative",
    "reviewed_status", "organism_result", "mapping_method", "mapping_status", "mapping_note",
]
with open(args.output_csv, "w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(output)

print(f"genes={len(genes)} exact={sum(x['mapping_status'] == 'exact' for x in output)} ambiguous={sum(x['mapping_status'] == 'ambiguous' for x in output)} unmatched={sum(x['mapping_status'] == 'unmatched' for x in output)}")

import argparse
import csv
import json
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request


parser = argparse.ArgumentParser(description="Fetch PubChem identity candidates for unresolved source compounds.")
parser.add_argument("screening_csv")
parser.add_argument("output_csv")
parser.add_argument("--cafile", required=True)
args = parser.parse_args()

with open(args.screening_csv, newline="", encoding="utf-8-sig") as handle:
    rows = [row for row in csv.DictReader(handle) if row.get("identity_match_status") != "exact"]

context = ssl.create_default_context(cafile=args.cafile)
seen = set()
output = []
for row in rows:
    name = (row.get("compound_name_raw") or "").strip()
    lookup_key = name.casefold()
    if not name or lookup_key in seen:
        continue
    seen.add(lookup_key)
    encoded = urllib.parse.quote(name, safe="")
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded}/property/Title,IUPACName,InChIKey,CanonicalSMILES,IsomericSMILES/JSON"
    )
    record = {
        "query_name": name,
        "pubchem_cid_candidate": "",
        "title": "",
        "iupac_name": "",
        "inchi_key": "",
        "canonical_smiles": "",
        "isomeric_smiles": "",
        "lookup_status": "unmatched",
        "source_url": url,
        "review_status": "manual review required",
    }
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "network-pharmacology-research-skill/1.0"})
        with urllib.request.urlopen(request, context=context, timeout=45) as response:
            payload = json.loads(response.read().decode("utf-8"))
        props = payload.get("PropertyTable", {}).get("Properties", [])
        if props:
            prop = props[0]
            record.update({
                "pubchem_cid_candidate": str(prop.get("CID", "")),
                "title": prop.get("Title", ""),
                "iupac_name": prop.get("IUPACName", ""),
                "inchi_key": prop.get("ConnectivitySMILES", "") if False else prop.get("InChIKey", ""),
                "canonical_smiles": prop.get("ConnectivitySMILES", prop.get("CanonicalSMILES", "")),
                "isomeric_smiles": prop.get("SMILES", prop.get("IsomericSMILES", "")),
                "lookup_status": "candidate found",
            })
    except urllib.error.HTTPError as error:
        record["lookup_status"] = f"HTTP {error.code}"
    except Exception as error:
        record["lookup_status"] = f"error: {type(error).__name__}"
    output.append(record)
    time.sleep(0.15)

fieldnames = list(output[0].keys()) if output else [
    "query_name", "pubchem_cid_candidate", "title", "iupac_name", "inchi_key",
    "canonical_smiles", "isomeric_smiles", "lookup_status", "source_url", "review_status",
]
with open(args.output_csv, "w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(output)

print(f"queries={len(output)} candidates={sum(row['lookup_status'] == 'candidate found' for row in output)}")

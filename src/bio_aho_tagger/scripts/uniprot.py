try:
    from .utils import download_file
except Exception:
    import importlib
    download_file = importlib.import_module("utils").download_file
import ahocorasick
import pickle
import gzip
import re


def parse_uniprot_dat(filepath, targets):
    """
    Parse the UniProt .dat.gz file and extract proteins for the requested organisms.

    Args:
        filepath: Path to the UniProt Swiss-Prot file (e.g., uniprot_sprot.dat.gz)
        targets: dict mapping a short key to a target organism definition, e.g.:
            {
              "human": {"os": "Homo sapiens", "taxid": "9606"},
              "mouse": {"os": "Mus musculus", "taxid": "10090"},
              "rat":   {"os": "Rattus norvegicus", "taxid": "10116"}
            }

    Returns:
        Dict[str, list[dict]]: mapping target key to list of protein dicts
    """
    results = {k: [] for k in targets.keys()}

    # precompile regexes
    recname_re = re.compile(r"DE\s+RecName:\s+Full=(.+);")
    altname_re = re.compile(r"DE\s+AltName:\s+Full=(.+);")
    ox_re = re.compile(r"OX\s+NCBI_TaxID=(\d+);")

    with gzip.open(filepath, "rt") as file:
        protein = {}
        matched_key = None  # which target this entry belongs to

        for line in file:
            # get the uniprot accession
            if line.startswith("AC"):
                # primary accession (first one before the semicolon)
                accessions = line.split()[1].strip(";")
                # flush previous protein if matched
                if protein and matched_key is not None:
                    results[matched_key].append(protein)

                # start a new protein record
                protein = {
                    "accession": accessions,
                    "preferred_name": "",
                    "synonyms": [],
                }
                matched_key = None

            # organism (OS): try to match by organism string
            if line.startswith("OS"):
                for key, spec in targets.items():
                    if spec.get("os") and spec["os"] in line:
                        matched_key = key
                        break

            # taxonomy id (OX): match by numeric taxid
            if line.startswith("OX"):
                m = ox_re.match(line)
                if m:
                    taxid = m.group(1)
                    for key, spec in targets.items():
                        if spec.get("taxid") == taxid:
                            matched_key = key
                            break

            # preferred name (RecName) - strip references in curly braces
            recname_match = recname_re.match(line)
            if recname_match:
                protein["preferred_name"] = re.sub(
                    r"\{.*?\}", "", recname_match.group(1)
                ).strip()

            # synonyms (AltName) - strip references in curly braces
            altname_match = altname_re.match(line)
            if altname_match:
                protein["synonyms"].append(
                    re.sub(r"\{.*?\}", "", altname_match.group(1)).strip()
                )

        # flush last
        if protein and matched_key is not None:
            results[matched_key].append(protein)

    return results


def add_proteins_to_automaton(automaton, proteins):
    """Add a list of protein dicts into an existing Aho-Corasick automaton."""
    for prot in proteins:
        acc = prot["accession"]
        pref_name = prot["preferred_name"]
        for s in prot["synonyms"]:
            syn = s.lower()
            # If a synonym appears in multiple organisms, the last added wins.
            automaton.add_word(syn, (syn, (pref_name, "Protein", f"uniprot:{acc}")))
    return automaton


def main():
    url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"
    filename = "uniprot_sprot.dat.gz"

    download_file(url, filename)

    targets = {
        "human": {"os": "Homo sapiens", "taxid": "9606"},
        "mouse": {"os": "Mus musculus", "taxid": "10090"},
        "rat": {"os": "Rattus norvegicus", "taxid": "10116"},
    }

    proteins_by_org = parse_uniprot_dat(filename, targets)

    # Build a single automaton with a defined indexing order:
    # 1) rat, 2) mouse, 3) human (later additions override earlier ones on duplicate synonyms)
    order = ["rat", "mouse", "human"]
    automaton = ahocorasick.Automaton()
    counts = {}
    for org in order:
        prots = proteins_by_org.get(org, [])
        counts[org] = len(prots)
        add_proteins_to_automaton(automaton, prots)

    automaton.make_automaton()

    out_file = "swissprot_rat_mouse_human.pkl"
    with open(out_file, "wb") as fh:
        pickle.dump(automaton, fh)
    print(
        f"Saved combined automaton to {out_file} (rat={counts['rat']}, mouse={counts['mouse']}, human={counts['human']})."
    )


if __name__ == "__main__":
    main()

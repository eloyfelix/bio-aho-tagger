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
    # OX taxonomy id
    ox_re = re.compile(r"OX\s+NCBI_TaxID=(\d+);")
    # Detect DE section starts
    de_recname_start = re.compile(r"DE\s+RecName:")
    de_altname_start = re.compile(r"DE\s+AltName:")
    # Extract tokens like Full=...; or Short=...; within a DE line
    de_token_re = re.compile(r"(Full|Short)=(.+?);")

    with gzip.open(filepath, "rt") as file:
        protein = {}
        matched_key = None  # which target this entry belongs to
        de_section = None  # current DE section context: 'RecName' | 'AltName' | None

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
                de_section = None

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

            # Parse DE lines capturing RecName Full (preferred) and Short (synonyms),
            # and AltName Full/Short as synonyms. Remove any references in {curly braces}.
            if line.startswith("DE"):
                # update section on explicit RecName/AltName starts
                if de_recname_start.match(line):
                    de_section = "RecName"
                elif de_altname_start.match(line):
                    de_section = "AltName"

                # extract tokens present in this line (works for both start and continuation lines)
                for tok, val in de_token_re.findall(line):
                    cleaned = re.sub(r"\{.*?\}", "", val).strip()
                    if tok == "Full" and de_section == "RecName":
                        # set preferred name (first Full encountered wins)
                        if not protein["preferred_name"]:
                            protein["preferred_name"] = cleaned
                    elif tok in ("Short", "Full") and de_section in ("RecName", "AltName"):
                        # treat RecName Short and all AltName names as synonyms
                        protein["synonyms"].append(cleaned)

        # flush last
        if protein and matched_key is not None:
            results[matched_key].append(protein)

    return results


def add_proteins_to_automaton(automaton, proteins):
    """Add a list of protein dicts into an existing Aho-Corasick automaton."""
    for prot in proteins:
        acc = prot["accession"]
        pref_name = prot.get("preferred_name", "")

        # Build list of terms to index: preferred name + synonyms, de-duplicated (case-insensitive)
        terms = []
        if pref_name:
            terms.append(pref_name)
        terms.extend(prot.get("synonyms", []))

        seen = set()
        for name in terms:
            if not name:
                continue
            key = name.lower()
            if key in seen:
                continue
            seen.add(key)
            # If a term appears in multiple organisms, the last added wins (order controlled by caller)
            automaton.add_word(key, (key, (pref_name, "Protein", f"uniprot:{acc}")))
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

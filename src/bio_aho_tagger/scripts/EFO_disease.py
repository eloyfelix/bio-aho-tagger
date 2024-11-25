from collections import defaultdict
from nltk.corpus import stopwords
import ahocorasick
import requests
import pickle
import nltk

nltk.download("stopwords")


sparql_query = """
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>

SELECT DISTINCT ?term ?label ?exactSynonym ?exactSynonymType ?narrowSynonym ?narrowSynonymType
WHERE {
    VALUES ?parentTerm { 
        <http://www.ebi.ac.uk/efo/EFO_0000408> 
        <http://purl.obolibrary.org/obo/HP_0000118> 
    }
    ?term rdfs:subClassOf* ?parentTerm .
    ?term rdfs:label ?label .
    FILTER (?term != ?parentTerm) .

    # Optional block for exact synonyms
    OPTIONAL {
        ?term oboInOwl:hasExactSynonym ?exactSynonym .
        OPTIONAL {
            ?exactSynonymAxiom
                owl:annotatedSource ?term ;
                owl:annotatedProperty oboInOwl:hasExactSynonym ;
                owl:annotatedTarget ?exactSynonym .
            OPTIONAL {
                ?exactSynonymAxiom oboInOwl:hasSynonymType ?exactSynonymType .
            }
        }
        FILTER NOT EXISTS {
            ?exactSynonymAxiom oboInOwl:hasSynonymType ?exactSynonymType .
            FILTER(?exactSynonymType IN (
                <http://purl.obolibrary.org/obo/mondo#ABBREVIATION>,
                <http://purl.obolibrary.org/obo/mondo/mondo-base#ABBREVIATION>,
                <http://purl.obolibrary.org/obo/hp#abbreviation>
            ))
        }
    }

    # Optional block for narrow synonyms
    OPTIONAL {
        ?term oboInOwl:hasNarrowSynonym ?narrowSynonym .
        OPTIONAL {
            ?narrowSynonymAxiom
                owl:annotatedSource ?term ;
                owl:annotatedProperty oboInOwl:hasNarrowSynonym ;
                owl:annotatedTarget ?narrowSynonym .
            OPTIONAL {
                ?narrowSynonymAxiom oboInOwl:hasSynonymType ?narrowSynonymType .
            }
        }
        FILTER NOT EXISTS {
            ?narrowSynonymAxiom oboInOwl:hasSynonymType ?narrowSynonymType .
            FILTER(?narrowSynonymType IN (
                <http://purl.obolibrary.org/obo/mondo#ABBREVIATION>,
                <http://purl.obolibrary.org/obo/mondo/mondo-base#ABBREVIATION>,
                <http://purl.obolibrary.org/obo/hp#abbreviation>
            ))
        }
    }
}
"""


def query_sparql(url, sparql_query):
    response = requests.get(url, params={"query": sparql_query, "format": "json"})
    data = response.json()
    return data


def format_data(data):
    diseases = defaultdict(lambda: {"label": None, "synonyms": defaultdict(set)})

    for row in data["results"]["bindings"]:
        ontology_id = ":".join(row["term"]["value"].split("/")[-1].split("_"))
        label = row["label"]["value"]

        # add the label
        diseases[ontology_id]["label"] = label

        # add exact synonyms
        exact_synonym = row.get("exactSynonym", {}).get("value")
        if exact_synonym:
            diseases[ontology_id]["synonyms"]["exact"].add(exact_synonym)

        # add narrow synonyms
        narrow_synonym = row.get("narrowSynonym", {}).get("value")
        if narrow_synonym:
            diseases[ontology_id]["synonyms"]["narrow"].add(narrow_synonym)
    return diseases


def main():

    stop_words = set(stopwords.words("english"))

    data = query_sparql("http://localhost:3030/efo/query", sparql_query)
    diseases = format_data(data)

    # 0 for label
    # 1 for exact synonym
    # 2 for narrow synonym
    automaton = ahocorasick.Automaton()

    for ontology_id, data in diseases.items():
        if not isinstance(ontology_id, str):
            continue
        label = data["label"]

        # add exact synonyms
        for syn in data["synonyms"]["exact"]:
            syn = syn.lower().replace("’", "'")
            if syn not in stop_words:
                automaton.add_word(syn, (syn, 1, (label, "Disease", ontology_id)))

        # add narrow synonyms
        for syn in data["synonyms"]["narrow"]:
            syn = syn.lower().replace("’", "'")
            if syn not in stop_words:
                automaton.add_word(syn, (syn, 2, (label, "Disease", ontology_id)))

        # add main term
        ll = label.lower().replace("’", "'")
        automaton.add_word(ll, (ll, 0, (label, "Disease", ontology_id)))

    automaton.make_automaton()

    with open("efo_disease.pkl", "wb") as file:
        pickle.dump(automaton, file)


if __name__ == "__main__":
    main()

import importlib.resources
import pickle


built_in_dicts = {
    "chembl_smiles": "chembl_smiles.pkl",
    "efo_disease": "efo_disease.pkl",
}


def merge_results(*lists):
    """
    Merge results from different automatons, keeping the longest match when there are overlaps,
    while allowing exact matches with different entity types to be retained.

    Example:
    disease_matches = [
        (0, 6, ('cancer', 'disease', 'MONDO:0004992')),
        (0, 11, ('lung cancer', 'disease', 'MONDO:0008903')),
        (20, 25, ('heart', 'disease', 'EFO:0003777')) # EFO:0003777 doesn't really include 'heart' synonym
    ]
    organ_matches = [
        (0, 4, ('lung', 'organ', 'UBERON:0002048')),
        (20, 25, ('heart', 'organ', 'UBERON:0000948'))
    ]
    merge_results(disease_matches, organ_matches)

    Expected Output:
    [
        (0, 11, ('lung cancer', 'disease', 'MONDO:0005070')),  # "lung" removed as it's a substring
        (20, 25, ('heart', 'disease', 'EFO:0003777')),      # "heart" kept as both "disease" and "organ"
        (20, 25, ('heart', 'organ', 'UBERON:0000948'))
    ]
    """
    all_matches = []
    for automaton_matches in lists:
        for start, end, (term, entity_type, entity_id) in automaton_matches:
            all_matches.append((start, end, term, entity_type, entity_id))

    # Sort matches by start position, prioritizing longer matches in case of overlap
    all_matches.sort(key=lambda x: (x[0], -(x[1] - x[0])))

    filtered_matches = []
    for start, end, term, entity_type, entity_id in all_matches:
        # Check if the current match is a true substring of any previously added match in filtered_matches
        if any(
            f_start <= start and end <= f_end and (f_start != start or f_end != end)
            for f_start, f_end, _ in filtered_matches
        ):
            # Skip this match if it is a true substring
            continue

        # Otherwise, add this match to the filtered list
        filtered_matches.append((start, end, (term, entity_type, entity_id)))

    return filtered_matches


class BioAhoTagger:

    stop_chars = [" ", ",", ".", "\n"]

    def __init__(self, file_path=None):
        self.automaton = self.load_automaton(file_path)

    def load_automaton(self, automaton):
        if not automaton:
            print(
                f"Use one of the built-in dictionaries: {[k for k in built_in_dicts.keys()]}\n"
                "e.g.: bta = BioAhoTagger('chembl_smiles')\n"
                "or use your own .pkl automaton file:\n"
                "e.g.: bta = BioAhoTagger('my_path/my_own_automaton.pkl')"
            )
            return None
        if automaton in built_in_dicts:
            # load from package resource
            with importlib.resources.open_binary(
                "bio_aho_tagger.data", built_in_dicts[automaton]
            ) as file:
                return pickle.load(file)
        else:
            # load from external file
            with open(automaton, "rb") as file:
                return pickle.load(file)

    def get(self, name):
        return self.automaton.get(name.lower(), None)

    def extract_entities(self, text):
        text = text.lower()
        entities = []
        for end_index, original_value in self.automaton.iter_long(text):
            start_index = end_index - len(original_value[0]) + 1
            if end_index + 1 <= len(text):
                if (end_index + 1 == len(text)) or (
                    text[end_index + 1] in self.stop_chars
                ):
                    if start_index == 0 or text[start_index - 1] in self.stop_chars:
                        entities.append((start_index, end_index + 1, (original_value)))
        return entities

import importlib.resources
import pickle


built_in_dicts = {
    "chembl_smiles": "chembl_smiles.pkl",
    "swissprot_rat_mouse_human": "swissprot_rat_mouse_human.pkl",
}


def merge_results(*lists):
    """
    Merge results from different automatons, keeping the longest match when there are overlaps,
    while allowing exact matches with different entity types to be retained.
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

    stop_chars = [" ", ",", ".", "\n", "\t", "<", ">", "(", ")", "/"]

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
        return self.automaton.get(name.lower().replace("’", "'"), None)

    def extract_entities(self, text):
        text = text.lower().replace("’", "'")
        text_length = len(text)
        entities = []

        for end_index, original_value in self.automaton.iter_long(text):
            start_index = end_index - len(original_value[0]) + 1
            next_char = text[end_index + 1] if end_index + 1 < text_length else None
            prev_char = text[start_index - 1] if start_index > 0 else None
            if (
                (next_char is None or next_char in self.stop_chars or next_char == ":")
                and not (next_char in {">", "("} and end_index != 0)
                and (prev_char is None or prev_char in self.stop_chars)
                and not (prev_char in {"<", ")"} and start_index != text_length - 1)
            ):
                entities.append((start_index, end_index + 1, original_value))

        return entities

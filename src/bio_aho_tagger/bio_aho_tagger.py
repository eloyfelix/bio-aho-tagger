import importlib.resources
import pickle


built_in_dicts = {
    "chembl_smiles": "chembl_smiles.pkl",
    "efo_disease": "efo_disease.pkl",
}


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

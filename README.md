# BioAhoTagger

Baseline dictionary based package to tag biological entities. Using [pyahocorasick](https://github.com/WojciechMula/pyahocorasick).

## To use it

Extract all entities found in a text (can use one of the built-in demo automatons or your own .pkl file)

```python
from bio_aho_tagger import BioAhoTagger

bt = BioAhoTagger("efo_disease")
disease = bt.get("lung cancer")

entities = bt.extract_entities(
    "Diabetes mellitus is a chronic condition characterized by elevated blood glucose levels, "
    "which can lead to complications such as cardiovascular disease and kidney failure. "
    "Hypertension, often called the \"silent killer,\" increases the risk of stroke, heart attack, "
    "and vision loss if left untreated. Early detection and management of diseases like asthma, "
    "rheumatoid arthritis, and chronic obstructive pulmonary disease can significantly improve "
    "quality of life."
)
```

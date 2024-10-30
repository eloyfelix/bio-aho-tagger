# BioAhoTagger

Baseline dictionary based package to tag biological entities. Using [pyahocorasick](https://github.com/WojciechMula/pyahocorasick).

## To use it

Extract all entities found in a text (can use one of the built-in demo automatons or your own .pkl file)

```python
from bio_aho_tagger import BioAhoTagger

bt = BioAhoTagger("efo_disease")
disease = bt.get("lung cancer")

entities = bt.extract_entities("The doctor prescribed metformin for managing diabetes and suggested amoxicillin to treat the bacterial infection.")
```

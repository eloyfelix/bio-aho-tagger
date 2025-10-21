# BioAhoTagger

Baseline dictionary based package to tag biological entities. Using [pyahocorasick](https://github.com/WojciechMula/pyahocorasick).

## To use it

Extract all entities found in a text (can use one of the built-in demo automatons or your own .pkl file)

```python
from bio_aho_tagger import BioAhoTagger

bt = BioAhoTagger("swissport_rat_mouse_human")
target = bt.get(target)

entities = bt.extract_entities(text)
```

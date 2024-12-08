from jpype import startJVM, getDefaultJVMPath, JPackage, shutdownJVM, isJVMStarted
from concurrent.futures import ThreadPoolExecutor
from importlib import resources


if not isJVMStarted():
    startJVM(
        getDefaultJVMPath(),
        f"-Djava.class.path={resources.files('bio_aho_tagger') / 'data/opsin-core-2.8.0-jar-with-dependencies.jar'}",
        convertStrings=False,
    )

opsin = JPackage("uk").ac.cam.ch.wwmm.opsin
NameToStructure = opsin.NameToStructure


def name_to_structure(name):
    nts = NameToStructure.getInstance()
    try:
        smiles = nts.parseToSmiles(name)
        return name, smiles or None
    except Exception as e:
        return name, None


def batch_process(names, num_threads=8):
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(name_to_structure, names))

    smiles_dict = {name: smiles for name, smiles in results if smiles}
    return smiles_dict

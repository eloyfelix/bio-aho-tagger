import regex
import unicodedata
from collections import defaultdict
from .bio_aho_tagger import BioAhoTagger
from .opsin import batch_process


def merge_results(*lists):
    all_matches = [match for automaton_matches in lists for match in automaton_matches]
    
    # Sort matches by start position, then by length (descending)
    all_matches.sort(key=lambda x: (x[0], -(x[1] - x[0])))

    filtered_matches = []
    for start, end, match_info in all_matches:
        # Skip this match if it's fully contained by any existing match
        if any(f_start <= start and end <= f_end for f_start, f_end, _ in filtered_matches):
            continue
        
        # Remove any previously added matches that are fully contained within this match
        filtered_matches = [
            existing_match for existing_match in filtered_matches
            if not (start <= existing_match[0] and existing_match[1] <= end)
        ]

        # Add the current match
        filtered_matches.append((start, end, match_info))
    return filtered_matches


class IUPACNameExtractor:
    def __init__(self):
        # organic chemistry terms and suffixes
        self.organic_terms = [
            'ane', 'ene', 'yne', 'ol', 'al', 'one', 'oic\s*acid', 
            'ide', 'ester', 'amide', 
            'acid', 'aldehyde', 'ketone', 'nitrile', 'isocyanate'
        ]
        
        # chemical prefixes and modifiers
        self.chemical_prefixes = [
            'di', 'tri', 'tetra', 'penta', 'hexa', 'hepta', 'octa', 'nona', 'deca',
            'iso', 'neo', 'cyclo', 'ortho', 'meta', 'para',
            'methyl', 'ethyl', 'propyl', 'butyl', 'pentyl', 'hexyl',
            'hydroxy', 'chloro', 'fluoro', 'bromo', 'iodo', 'amino', 
            'nitro', 'sulfono', 'phosphono', 'isocyanato', 'diisocyanato'
        ]
        
        # common root name fragments for organic compounds
        self.root_fragments = [
            'meth', 'eth', 'prop', 'but', 'pent', 'hex', 'hept', 'oct', 'non', 'dec',
            'undec', 'dodec', 'tridec', 'tetradec',
            'cyclo', 'benz', 'tolu', 'phenyl', 'acetyl', 'carb', 'alkyl', 
            'isocyanato', 'diisocyanato'
        ]
        
        # common non-chemical word endings to filter out
        self.non_chemical_endings = [
            'al', 'ial', 'ual', 'onal', 'ical', 'itive', 
            'ive', 'ible', 'able', 'ment', 'guide', 'zone', 
            'sis', 'ous', 'ity', 'age', 'ure', 'equal'
        ]
        self.iupac_pattern = self._build_iupac_pattern()
    
    def _build_iupac_pattern(self):
        """
        Build an advanced regex pattern for IUPAC name extraction
        with multiple validation checks
        """
        # pattern that ensures chemical validity
        pattern = regex.compile(
            r'''
            (?<!\w)  # not part of another word
            (?P<name>
                (?:
                    (?:\d+(?:,\d+)*-)+  # multi-position matching like "1,4-" or "1-3-"
                )?
                (?:''' + '|'.join(map(regex.escape, self.chemical_prefixes)) + r''')?   # optional chemical prefixes
                (?:''' + '|'.join(map(regex.escape, self.root_fragments)) + r''')?     # optional root fragments
                [A-Za-z]+                     # main root name
                (?:-?\d+)?                    # additional position indicators
                (?:-?(?:''' + '|'.join(map(regex.escape, self.chemical_prefixes)) + r'''))*  # possible modifiers
                (?:''' + '|'.join(map(regex.escape, self.organic_terms)) + r''')       # required suffixes
            )
            \b  # word boundary to prevent partial matches
            ''', 
            regex.VERBOSE | regex.IGNORECASE | regex.UNICODE
        )
        
        return pattern
    
    def normalize_name(self, name):
        """
        Normalize chemical name to remove diacritics and standardize
        """
        # Convert to lowercase and normalize unicode
        normalized = unicodedata.normalize('NFKD', name.lower())
        return ''.join(
            char for char in normalized 
            if not unicodedata.combining(char)
        )
    
    def is_valid_chemical_name(self, name):
        """
        Validation for chemical names with more robust checks
        """
        # Normalize for validation checks
        name_lower = self.normalize_name(name)
        
        # Minimum and maximum length requirements
        if len(name) < 4 or len(name) > 50:
            return False
        
        # Filter out non-chemical word endings
        if any(name_lower.endswith(ending) for ending in self.non_chemical_endings):
            return False
        
        # Validate position indicators more rigorously
        position_parts = name.split('-')
        
        # Strict check for position indicators
        first_part = position_parts[0]
        if ',' in first_part:
            # Ensure the first part is a valid multi-position indicator
            if not regex.match(r'^\d+(?:,\d+)+$', first_part):
                return False
        
        # Check for chemical-specific elements
        has_chemical_root = any(
            root in name_lower 
            for root in self.root_fragments + self.chemical_prefixes
        )
        
        # Check if name ends with a valid organic term
        has_valid_suffix = any(
            regex.search(term + r'\b', name_lower) 
            for term in self.organic_terms
        )
        
        # Additional chemical name characteristics
        chemical_indicators = [
            'cyclo', 'di', 'tri', 'tetra', 
            'hydroxy', 'chloro', 'fluoro', 
            'amino', 'nitro', 'isocyanato'
        ]
        
        has_chemical_indicator = any(
            indicator in name_lower 
            for indicator in chemical_indicators
        )
        
        # Ensure alphanumeric content is reasonable
        alpha_ratio = sum(c.isalpha() for c in name) / len(name)
        
        # Final validation
        return (
            (has_chemical_root or has_valid_suffix or has_chemical_indicator) 
            and alpha_ratio > 0.6
        )
    
    def extract_iupac_names(self, text):
        """
        Extract potential IUPAC names with their positions in the original text
        
        Returns:
            Dict of {name: [(start, end), ...]}
        """
        # Find potential matches with positions in original text
        matches = list(self.iupac_pattern.finditer(text, overlapped=True))
        
        # Group matches efficiently
        name_positions = defaultdict(list)
        for match in matches:
            name = match.group('name')
            
            # Validate the potential chemical name
            if self.is_valid_chemical_name(name):
                name_positions[name].append((
                    match.start('name'),  # Start position in original text
                    match.end('name')     # End position in original text
                ))
        
        # Convert defaultdict to regular dict and sort positions
        return {
            name: sorted(set(positions))  # Remove potential duplicates
            for name, positions in name_positions.items()
        }


def filter_and_format_chemicals(chemical_dict, chemical_data):
    result = []

    for chemical_name, positions in chemical_dict.items():
        if chemical_name in chemical_data.keys():
            smiles = chemical_data.get(chemical_name, '')
            for start, end in positions:
                result.append((start, end, (chemical_name, (chemical_name, 'Chemical', smiles))))

    return result


def extract_chem_entities(text):
    # extract potential IUPAC names
    extractor = IUPACNameExtractor()
    entities = extractor.extract_iupac_names(text)

    # run the potential IUPAC names through OPSIN
    opsin_dict = batch_process(entities)

    # keep only the ones that OPSIN managed to resolve
    filtered = filter_and_format_chemicals(entities, opsin_dict)

    # get trivial names using chembl synonyms
    bat = BioAhoTagger("chembl_smiles")
    bat_entities = bat.extract_entities(text)

    # merge both sets of compounds
    merged = merge_results(bat_entities, filtered)
    return merged

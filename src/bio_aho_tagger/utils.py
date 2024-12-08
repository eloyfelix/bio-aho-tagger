
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


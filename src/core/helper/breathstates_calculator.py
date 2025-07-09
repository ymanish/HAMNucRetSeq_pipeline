from typing import List, Tuple

def total_states_index(length: int) -> List[Tuple[int, int]]:
    """
    Generate all possible binding states for 'b_index' or 'ph_index' styles. 
    Meaning, it gives the states in terms of index positions from left and right.
    So left=0, right=13 for the "b_index" style means the nucleosme is completely wrapped,
    and left=0, right=27 for the "ph_index" style means the nucleosome is completely wrapped as well.
    Returns pairs (left, right) where 0 <= left <= right < length, representing
    nucleosome breathing states based on index positions.
    
    Args:
        length (int): The maximum index position (exclusive).
    
    Returns:
        List[Tuple[int, int]]: List of (left, right) state pairs.
    """
    
    states = []
    for left in range(length):
        for right in range(left, length):
            states.append((left, right))
    return states

def total_open_states(length: int) -> List[Tuple[int, int]]:
    """
    Generate all possible open states for 'open_sites' style. Meaning, it gives the states in terms of open sites from left and right.
    Returns pairs (left, right) where left >= 0, right >= 0, and left + right <= length,
    representing open sites in nucleosome breathing.
    
    Args:
        length (int): The maximum sum of left and right positions (inclusive).
    
    Returns:
        List[Tuple[int, int]]: List of (left, right) state pairs.
    """

    open_states = []
    for left in range(length + 1):
        for right in range(length + 1):
            if left + right <= length:
                open_states.append((left, right))
    return open_states


if __name__ == "__main__":
        
    STYLE = "ph_index"  # or "open_sites" or "ph_index"

    if STYLE == "open_sites":
        states = total_open_states()
    elif STYLE == "b_index":
        states = total_states_index(length=14)
    else:
        states = total_states_index(length=28)
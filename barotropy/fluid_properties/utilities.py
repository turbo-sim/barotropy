
import numpy as np


def is_float(element: any) -> bool:
    """
    Check if the given element can be converted to a float.

    Parameters
    ----------
    element : any
        The element to be checked.

    Returns
    -------
    bool
        True if the element can be converted to a float, False otherwise.
    """

    if element is None:
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False
    
    
def _generate_coolprop_input_table():
    """Create table of input pairs as string to be copy-pasted in Sphinx documentation"""
    inputs_table = ".. list-table:: CoolProp input mappings\n"
    inputs_table += "   :widths: 50 30\n"
    inputs_table += "   :header-rows: 1\n\n"
    inputs_table += "   * - Input pair name\n"
    inputs_table += "     - Input pair mapping\n"
    from .high_level import INPUT_PAIRS
    for name, value in INPUT_PAIRS:
        inputs_table += f"   * - {name}\n"
        inputs_table += f"     - {value}\n"

    return inputs_table


def states_to_dict(states):
    """
    Convert a list of state objects into a dictionary.
    Each key is a field name of the state objects, and each value is a NumPy array of all the values for that field.
    """
    state_dict = {}
    for field in states[0].keys():
        state_dict[field] = np.array([getattr(state, field) for state in states])
    return state_dict


def states_to_dict_2d(states):
    """
    Convert a 2D list (grid) of state objects into a dictionary.
    Each key is a field name of the state objects, and each value is a 2D NumPy array of all the values for that field.

    Parameters
    ----------
    states_grid : list of list of objects
        A 2D grid where each element is a state object with the same keys.

    Returns
    -------
    dict
        A dictionary where keys are field names and values are 2D arrays of field values.
    """
    state_dict_2d = {}
    for i, row in enumerate(states):
        for j, state in enumerate(row):
            for field in state.keys():
                if field not in state_dict_2d:
                    m, n = len(states), len(row)
                    state_dict_2d[field] = np.empty((m, n), dtype=object)
                state_dict_2d[field][i, j] = state[field]

    return state_dict_2d
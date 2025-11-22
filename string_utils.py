def split_before_uppercases(formula):
    start = 0
    end = 1
    elements_lst = []

    if not formula:
        return elements_lst

    while end < len(formula):
        if formula[end].isupper():
            elements_lst.append(formula[start:end])
            start = end
        end+=1

    elements_lst.append(formula[start:])

    return elements_lst


def split_at_digit(formula):
    digit_location = 0

    while digit_location < len(formula) and not formula[digit_location].isdigit():
        digit_location += 1

    if digit_location == len(formula):
        return formula, 1

    num_end = digit_location
    while num_end < len(formula) and formula[num_end].isdigit():
        num_end += 1

    prefix = formula[:digit_location]
    numeric_part_str = formula[digit_location:num_end]
    numeric_part_int = int(numeric_part_str)

    return prefix, numeric_part_int



def count_atoms_in_molecule(molecular_formula):
    """Takes a molecular formula (string) and returns a dictionary of atom counts.
    Example: 'H2O' → {'H': 2, 'O': 1}"""

    # Step 1: Initialize an empty dictionary to store atom counts
    atoms_count_dict = {}

    for atom in split_before_uppercases(molecular_formula):
        atom_name, atom_count = split_at_digit(atom)
        # Step 2: Update the dictionary with the atom name and count
        atoms_count_dict[atom_name] = atoms_count_dict.get(atom_name, 0) + atom_count

    # Step 3: Return the completed dictionary
    return atoms_count_dict
    
def parse_chemical_reaction(reaction_equation):
    """Takes a reaction equation (string) and returns reactants and products as lists.  
    Example: 'H2 + O2 -> H2O' → (['H2', 'O2'], ['H2O'])"""
    reaction_equation = reaction_equation.replace(" ", "")  # Remove spaces for easier parsing
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")
    
def count_atoms_in_reaction(molecules_list):
    """Takes a list of molecular formulas and returns a list of atom count dictionaries.  
    Example: ['H2', 'O2'] → [{'H': 2}, {'O': 2}]"""
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count


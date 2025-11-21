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
    for char_index, char in enumerate(formula):
        if char.isdigit():
            return formula[:char_index], int(formula[char_index:])
    return formula, 1

def count_atoms_in_molecule(molecular_formula):
    atoms_count_dict = {}
    for atom in split_before_uppercases(molecular_formula):
        atom_name, atom_count = split_at_digit(atom)
        atoms_count_dict[atom_name] = atoms_count_dict.get(atom_name, 0) + atom_count    
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

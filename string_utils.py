def split_at_digit(formula):
"""
Splits an element string into (element_name, count).
Example: 'H2' → ('H', 2), 'O' → ('O', 1)
"""
if not formula:
return "", 1

for idx, char in enumerate(formula):
if char.isdigit():
start = idx
end = idx
while end < len(formula) and formula[end].isdigit():
end += 1
return formula[:start], int(formula[start:end])

return formula, 1 # No digit → count = 1


def split_before_uppercases(formula):
"""
Splits a chemical formula based on uppercase letters.
Example: 'C6H12O6' → ['C6', 'H12', 'O6']
"""
if not formula:
return []

split_formula = []
start = 0

for i in range(1, len(formula)):
if formula[i].isupper():
split_formula.append(formula[start:i])
start = i

split_formula.append(formula[start:])
return split_formula


def count_atoms_in_molecule(molecular_formula):
"""
Takes a molecular formula (string) and returns a dictionary of atom counts.
Example: 'H2O' → {'H': 2, 'O': 1}
"""

# Step 1: Initialize a dictionary
atom_counts = {}

# Step 2: Split the molecule into pieces like ['H2', 'O']
for atom in split_before_uppercases(molecular_formula):
atom_name, atom_count = split_at_digit(atom)

# Step 3: Update counts
if atom_name in atom_counts:
atom_counts[atom_name] += atom_count
else:
atom_counts[atom_name] = atom_count

# Step 4: Return final counts
return atom_counts


# -----------------------------
# Reaction utility functions
# -----------------------------

def parse_chemical_reaction(reaction_equation):
"""
Example: 'H2 + O2 -> H2O'
Returns: (['H2', 'O2'], ['H2O'])
"""
reaction_equation = reaction_equation.replace(" ", "")
reactants, products = reaction_equation.split("->")
return reactants.split("+"), products.split("+")


def count_atoms_in_reaction(molecules_list):
"""
Takes a list like ['H2', 'O2'] and returns:
[{'H': 2}, {'O': 2}]
"""
molecules_atoms_count = []
for molecule in molecules_list:
molecules_atoms_count.append(count_atoms_in_molecule(molecule))
return molecules_atoms_count



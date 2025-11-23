from sympy import symbols, Eq, solve as sympy_solve

ELEMENTS = [
'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
'Sb', 'I', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
'Rg', 'Cn', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo'
]


def generate_equation_for_element(compounds, coefficients, element):
"""Generate symbolic equation representing atom balance for one element."""
equation = 0
for i, compound in enumerate(compounds):
    if element in compound:
equation += coefficients[i] * compound[element]
return equation


def build_equations(reactant_atoms, product_atoms):
"""Build symbolic equations for each element to balance the reaction."""

# Create coefficients: a0...a(n-1) for reactants, b0...b(m-1) for products
reactant_coefficients = list(symbols(f'a0:{len(reactant_atoms)}'))
product_coefficients = list(symbols(f'b0:{len(product_atoms)}'))

# Set last product coefficient to 1 (fixing the scale)
product_coefficients = product_coefficients[:-1] + [1]

equations = []

# Build equation for every element appearing in either side
for element in ELEMENTS:
lhs = generate_equation_for_element(reactant_atoms, reactant_coefficients, element)
rhs = generate_equation_for_element(product_atoms, product_coefficients, element)
    if lhs != 0 or rhs != 0:
equations.append(Eq(lhs, rhs))

# return equations + all coefficients except the fixed b coefficient
return equations, reactant_coefficients + product_coefficients[:-1]


def my_solve(equations, coefficients):
"""Solve the system of linear equations for the reaction coefficients."""

solution = sympy_solve(equations, coefficients)

# If number of solved vars does NOT match → incomplete solution
    if len(solution) != len(coefficients):
return None 

result = []
for c in coefficients:
result.append(float(solution[c]))

return result

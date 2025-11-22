import math
import sympy

def parse_chemical_reaction(reaction_equation):
    """Takes a reaction equation (string) and returns reactants and products as lists.  
    Example: 'H2 + O2 -> H2O' → (['H2', 'O2'], ['H2O'])"""
    reaction_equation = reaction_equation.replace(" ", "")  # Remove spaces for easier parsing
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")

def count_atoms_in_reaction(molecules_list):
    """Takes a list of molecular formulas and returns a dictionary of total atom counts across all molecules."""
    total_atoms = {}
    for molecule_formula in molecules_list:
        atoms_in_mol = count_atoms_in_molecule(molecule_formula) # Use the updated count_atoms_in_molecule
        for element, count in atoms_in_mol.items():
            total_atoms[element] = total_atoms.get(element, 0) + count
    return total_atoms
    
def balance_reaction(equation_str):
    """Balances a chemical reaction equation and returns the coefficients."""
    # Use user's parse_chemical_reaction
    reactants_str_list, products_str_list = parse_chemical_reaction(equation_str)

    # Convert string lists to atom count dictionaries using the updated count_atoms_in_molecule
    reactants = [count_atoms_in_molecule(c) for c in reactants_str_list]
    products = [count_atoms_in_molecule(c) for c in products_str_list]

    num_reactants = len(reactants)
    num_products = len(products)

    # Create symbolic variables for coefficients
    r_coeffs = sympy.symbols(f'r0:{num_reactants}')
    p_coeffs = sympy.symbols(f'p0:{num_products}')

    all_elements = set()
    for compound in reactants + products:
        all_elements.update(compound.keys())

    equations = []
    for element in all_elements:
        # Sum of element atoms on reactant side
        reactant_sum = generate_equation_for_element(reactants, r_coeffs, element)
        # Sum of element atoms on product side
        product_sum = generate_equation_for_element(products, p_coeffs, element)
        equations.append(sympy.Eq(reactant_sum, product_sum))

    # Add a normalization equation: set the first reactant coefficient to 1 for initial solving
    if r_coeffs:
        equations.append(sympy.Eq(r_coeffs[0], 1))
    else:
        raise ValueError("No reactants found, cannot normalize equation.")

    # Solve the system of equations
    all_coeffs_symbols = list(r_coeffs) + list(p_coeffs)
    solution = sympy.solve(equations, all_coeffs_symbols)

    if not solution:
        raise ValueError("Could not find a solution to balance the equation.")

    # Extract numerical coefficients as SymPy Rational objects
    coeffs_values = [solution[sym] for sym in all_coeffs_symbols]

    # Find the least common multiple (LCM) of all denominators to scale all coefficients to integers,
    # then simplify by GCD.

    rational_coeffs = [c for c in coeffs_values if isinstance(c, sympy.Rational)]
    denominators = [c.q for c in rational_coeffs if c.q != 1]

    lcm_val = 1
    if denominators:
        lcm_val = denominators[0]
        for denom in denominators[1:]:
            lcm_val = math.lcm(lcm_val, denom)
    
    scaled_coeffs = [c * lcm_val for c in coeffs_values]
    
    integer_numerators = [int(c) for c in scaled_coeffs]

    # Find GCD of all integer numerators to get the smallest integer coefficients
    gcd_val = 1
    if integer_numerators: # Ensure list is not empty before computing GCD
        non_zero_numerators = [abs(n) for n in integer_numerators if n != 0]
        if non_zero_numerators:
            gcd_val = non_zero_numerators[0]
            for n in non_zero_numerators[1:]:
                gcd_val = math.gcd(gcd_val, n)
    
    simplified_integer_coeffs = [n // gcd_val for n in integer_numerators]

    # Normalize by the coefficient of the LAST PRODUCT to match test expectations.
    # If there are products, find the coefficient corresponding to the last product.
    if num_products > 0:
        normalizing_coeff_index = num_reactants + num_products - 1
        normalizer = simplified_integer_coeffs[normalizing_coeff_index]

        if normalizer != 0:
            final_coeffs = [sympy.Rational(n, normalizer) for n in simplified_integer_coeffs]
        else:
            # Fallback if the normalizing coefficient is zero (should not happen for a balanced equation)
            final_coeffs = [sympy.Rational(n, 1) for n in simplified_integer_coeffs]
    elif num_reactants > 0:
        # If no products, normalize by the first reactant's coefficient (already set to 1 initially)
        # This branch ensures it handles cases where only reactants are present, though not typical for balancing.
        normalizer = simplified_integer_coeffs[0]
        if normalizer != 0:
            final_coeffs = [sympy.Rational(n, normalizer) for n in simplified_integer_coeffs]
        else:
            final_coeffs = [sympy.Rational(n, 1) for n in simplified_integer_coeffs]
    else:
        # No reactants or products (should have been caught by earlier ValueError)
        final_coeffs = [] 

    return final_coeffs

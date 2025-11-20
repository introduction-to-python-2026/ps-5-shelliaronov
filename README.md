[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/I65RxAd9)
# Chemical Reaction Balancing Assignment

Welcome to the chemical reaction balancing assignment! Your task is to complete the provided Python project by filling in the missing code and writing the necessary imports. Follow the instructions carefully to ensure all components work together seamlessly.

---

## Instructions

### 1. Complete the Missing Functions in `string_utils.py`
In the file `string_utils.py`, you will find placeholders for two functions:  
- `split_before_uppercases`  
- `split_at_digit`  

These functions were part of last week's assignment. Use your previous implementation or write new ones based on the requirements:  
- **`split_before_uppercases`**: Splits a chemical formula (string) into a list of elements based on uppercase letters.  
- **`split_at_digit`**: Splits an element string into a tuple of (element name, count).  

---

### 2. Fill in Missing Lines in `count_atoms_in_molecule`
In the function `count_atoms_in_molecule` inside `string_utils.py`, fill in the missing lines. Here's the partially implemented function:

```python
def count_atoms_in_molecule(molecular_formula):
    """Takes a molecular formula (string) and returns a dictionary of atom counts.  
    Example: 'H2O' → {'H': 2, 'O': 1}"""

    # Step 1: Initialize an empty dictionary to store atom counts

    for atom in split_by_capitals(molecular_formula):
        atom_name, atom_count = split_at_number(atom)
        
        # Step 2: Update the dictionary with the atom name and count

    # Step 3: Return the completed dictionary
```

### Your tasks are:

1. Add a line to initialize an empty dictionary to store atom counts.  
2. Add a line to update the dictionary with the atom name and count.  
3. Add a line to return the completed dictionary.  

---

### 3. Write All Import Statements
In `main.py`, import all the necessary functions from `string_utils.py` and `equation_utils.py` so that `balance_reaction` works correctly.

- Identify which functions are needed from each file.  
- Write the appropriate import statements at the top of `main.py`.

In `equation_utils.py`, ensure that the necessary functions from `sympy` are imported.  
Note: The function `my_solve` uses the `solve` function from `sympy` but is referenced as `sympy_solve`.  


- Identify which functions are needed from each file.  
- Write the appropriate import statements at the top of `main.py`.  


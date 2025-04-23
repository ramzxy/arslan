import re
import numpy as np
from fractions import Fraction

def parse_formula(formula):
    """Extract elements and their counts from a chemical formula."""
    formula = formula.strip()
    
    pattern = r'([A-Z][a-z]*)(\d*)'
    elements = {}
    for element, count in re.findall(pattern, formula):
        count = int(count) if count else 1
        elements[element] = elements.get(element, 0) + count
    
    return elements

def extract_charge(compound):
    """Extract charge from an ionic formula (Fe+3, SO4-2, etc.)."""
    charge = 0
    base_formula = compound
    
    if '+' in compound:
        parts = compound.split('+')
        base_formula = parts[0]
        if len(parts) > 1 and parts[1].isdigit():
            charge = int(parts[1])
        else:
            charge = 1
    elif '-' in compound:
        parts = compound.split('-')
        base_formula = parts[0]
        if len(parts) > 1 and parts[1].isdigit():
            charge = -int(parts[1])
        else:
            charge = -1
    
    return base_formula, charge

def parse_equation(equation):
    """Split a chemical equation into reactants and products."""
    equation = re.sub(r'(\w+)\^(\d*)([+-])', r'\1\3\2', equation)
    
    sides = equation.split('->')
    if len(sides) != 2:
        raise ValueError("Invalid equation format. Use format 'A + B -> C + D'")
    
    reactants = [comp.strip() for comp in sides[0].split(' + ')]
    products = [comp.strip() for comp in sides[1].split(' + ')]
    
    # Remove any empty strings
    reactants = [r for r in reactants if r]
    products = [p for p in products if p]
    
    return reactants, products

def build_compound_data(compounds):
    """Create a dictionary with element and charge data for each compound."""
    compound_data = {}
    
    for compound in compounds:
        base_formula, charge = extract_charge(compound)
        elements = parse_formula(base_formula)
        compound_data[compound] = {
            'elements': elements,
            'charge': charge,
            'base_formula': base_formula
        }
    
    return compound_data

def build_matrix(reactants, products, all_elements, compound_data):
    """Create coefficient matrix for the chemical equation system."""
    n_compounds = len(reactants) + len(products)
    matrix = np.zeros((len(all_elements), n_compounds))
    
    # Fill in coefficients for reactants (positive values)
    for i, reactant in enumerate(reactants):
        data = compound_data[reactant]
        for element, count in data['elements'].items():
            if element in all_elements:
                row = all_elements.index(element)
                matrix[row, i] = count
        
        # Add charge if present
        if 'charge' in all_elements and data['charge'] != 0:
            row = all_elements.index('charge')
            matrix[row, i] = data['charge']
    
    # Fill in coefficients for products (negative values)
    for i, product in enumerate(products):
        data = compound_data[product]
        for element, count in data['elements'].items():
            if element in all_elements:
                row = all_elements.index(element)
                matrix[row, i + len(reactants)] = -count
        
        # Add charge if present
        if 'charge' in all_elements and data['charge'] != 0:
            row = all_elements.index('charge')
            matrix[row, i + len(reactants)] = -data['charge']
    
    return matrix

def gaussian_elimination(matrix):
    """Reduce matrix to row echelon form using Gaussian elimination."""
    m, n = matrix.shape
    r = 0  # Current row
    
    for c in range(n):
        # Find the pivot row
        pivot_row = None
        for i in range(r, m):
            if matrix[i, c] != 0:
                pivot_row = i
                break
        
        if pivot_row is None:
            continue  # No pivot in this column
        
        # Swap rows to bring the pivot to the current row
        if pivot_row != r:
            matrix[[r, pivot_row]] = matrix[[pivot_row, r]]
        
        # Normalize the pivot row
        pivot = matrix[r, c]
        matrix[r] = matrix[r] / pivot
        
        # Eliminate other rows
        for i in range(m):
            if i != r and matrix[i, c] != 0:
                factor = matrix[i, c]
                matrix[i] = matrix[i] - factor * matrix[r]
        
        r += 1
        if r == m:
            break
    
    return matrix

def solve_coefficients(matrix):
    """Calculate integer coefficients from the reduced row echelon form."""
    m, n = matrix.shape
    
    # Set last variable to 1
    solution = np.zeros(n)
    solution[-1] = 1
    
    # Back substitution
    for i in range(m-1, -1, -1):
        # Find the pivot column for this row
        pivot_col = -1
        for j in range(n-1):
            if matrix[i, j] != 0:
                pivot_col = j
                break
        
        if pivot_col != -1:
            # Calculate the variable value
            rhs = 0
            for j in range(pivot_col+1, n):
                rhs += matrix[i, j] * solution[j]
            solution[pivot_col] = -rhs / matrix[i, pivot_col]
    
    # Convert to integers using LCM
    fractions = [Fraction(float(s)).limit_denominator() for s in solution]
    lcm = 1
    for f in fractions:
        lcm = lcm * f.denominator // np.gcd(lcm, f.denominator)
    
    coefficients = [int(f * lcm) for f in fractions]
    
    # Simplify by dividing by GCD if possible
    gcd = coefficients[0]
    for c in coefficients[1:]:
        gcd = np.gcd(gcd, c)
    
    if gcd > 1:
        coefficients = [c // gcd for c in coefficients]
    
    return coefficients

def format_balanced_equation(reactants, products, coefficients):
    """Format the balanced equation as a string."""
    n_reactants = len(reactants)
    
    # Format reactants
    reactant_terms = []
    for i, reactant in enumerate(reactants):
        coef = coefficients[i]
        if coef == 0:
            continue
        if coef == 1:
            reactant_terms.append(reactant)
        elif coef > 0:
            reactant_terms.append(f"{coef}{reactant}")
    
    # Format products
    product_terms = []
    for i, product in enumerate(products):
        coef = coefficients[i + n_reactants]
        if coef == 0:
            continue
        if coef == 1:
            product_terms.append(product)
        elif coef > 0:
            product_terms.append(f"{coef}{product}")
    
    return " + ".join(reactant_terms) + " -> " + " + ".join(product_terms)

def display_compound_data(compound_data):
    """Show parsed compound information."""
    print("\nCompound analysis:")
    for compound, data in compound_data.items():
        elements_with_charge = data['elements'].copy()
        if data['charge'] != 0:
            elements_with_charge['charge'] = data['charge']
        print(f"{compound}: {elements_with_charge}")

def check_balance(reactants, products, all_elements, coefficients, compound_data):
    """Verify if equation is balanced by comparing element counts on both sides."""
    n_reactants = len(reactants)
    
    # Calculate vectors for each side
    reactant_vectors = []
    product_vectors = []
    
    # Calculate element counts for reactants
    total_reactants = [0] * len(all_elements)
    for i, reactant in enumerate(reactants):
        data = compound_data[reactant]
        vector = []
        for elem in all_elements:
            if elem == 'charge':
                vector.append(data['charge'])
            else:
                vector.append(data['elements'].get(elem, 0))
        
        scaled_vector = [v * coefficients[i] for v in vector]
        reactant_vectors.append((reactant, scaled_vector, coefficients[i]))
        total_reactants = [total_reactants[j] + scaled_vector[j] for j in range(len(all_elements))]
    
    # Calculate element counts for products
    total_products = [0] * len(all_elements)
    for i, product in enumerate(products):
        data = compound_data[product]
        vector = []
        for elem in all_elements:
            if elem == 'charge':
                vector.append(data['charge'])
            else:
                vector.append(data['elements'].get(elem, 0))
        
        scaled_vector = [v * coefficients[i + n_reactants] for v in vector]
        product_vectors.append((product, scaled_vector, coefficients[i + n_reactants]))
        total_products = [total_products[j] + scaled_vector[j] for j in range(len(all_elements))]
    
    # Print LHS and RHS
    print("\nLHS:", end=" ")
    for i, (compound, _, coef) in enumerate(reactant_vectors):
        if coef == 0:
            continue
        if i > 0 and coef > 0:
            print(" +", end=" ")
        if coef == 1:
            print(f"{compound}", end="")
        else:
            print(f"{coef}{compound}", end="")
    print(f" = {total_reactants}")
    
    print("RHS:", end=" ")
    for i, (compound, _, coef) in enumerate(product_vectors):
        if coef == 0:
            continue
        if i > 0 and coef > 0:
            print(" +", end=" ")
        if coef == 1:
            print(f"{compound}", end="")
        else:
            print(f"{coef}{compound}", end="")
    print(f" = {total_products}")
    
    # Check if balanced
    is_balanced = total_reactants == total_products
    if is_balanced:
        print("\n=> Vectors are equal => Equation is balanced")
    else:
        print("\n=> Vectors are not equal => Equation is not balanced")
    
    return is_balanced

def display_element_vectors(compounds, all_elements, compound_data):
    """Display element vectors for each compound."""
    print("\nElement vectors for each compound:")
    
    for compound in compounds:
        data = compound_data[compound]
        vector = []
        for elem in all_elements:
            if elem == 'charge':
                vector.append(data['charge'])
            else:
                vector.append(data['elements'].get(elem, 0))
        
        print(f"{compound} -> {vector}")

def balance_equation(equation):
    """Balance chemical equations including those with ionic charges."""
    print(f"Balancing equation: {equation}")
    
    # Parse equation
    reactants, products = parse_equation(equation)
    
    # Extract elements and charges
    compound_data = build_compound_data(reactants + products)
    display_compound_data(compound_data)
    
    # Collect all elements
    all_elements = set()
    has_charge = False
    
    for data in compound_data.values():
        all_elements.update(data['elements'].keys())
        if data['charge'] != 0:
            has_charge = True
    
    # Add charge if needed
    all_elements = sorted(list(all_elements))
    if has_charge:
        all_elements.append('charge')
    
    print(f"\nElements tracked: {all_elements}")
    
    # Build and solve matrix
    matrix = build_matrix(reactants, products, all_elements, compound_data)
    
    print("\nMatrix representation of the equation:")
    print(matrix)
    
    row_echelon = gaussian_elimination(matrix.copy())
    
    print("\nRow Echelon Form:")
    print(row_echelon)
    
    coefficients = solve_coefficients(row_echelon)
    
    # Display results
    all_compounds = reactants + products
    display_element_vectors(all_compounds, all_elements, compound_data)
    check_balance(reactants, products, all_elements, coefficients, compound_data)
    
    balanced_equation = format_balanced_equation(reactants, products, coefficients)
    
    return balanced_equation

# Example usage
if __name__ == "__main__":
    print("\n=== Example 1: Ethane Combustion ===")
    equation = "C2H6 + O2 -> CO2 + H2O"
    balanced = balance_equation(equation)
    print(f"\nBalanced equation: {balanced}")
    
    print("\n" + "="*50)
    print("=== Example 2: Iron Oxidation ===")
    equation = "Fe + O2 -> Fe2O3"
    balanced = balance_equation(equation)
    print(f"\nBalanced equation: {balanced}")
    
    print("\n" + "="*50)
    print("=== Example 3: Ammonia Oxidation ===")
    equation = "NH3 + O2 -> NO + H2O"
    balanced = balance_equation(equation)
    print(f"\nBalanced equation: {balanced}")
    
    print("\n" + "="*50)
    print("=== Example 4: Ionic Reaction (Permanganate + Oxalate) ===")
    equation = "MnO4- + C2O4-2 + H+ -> Mn+2 + CO2 + H2O"
    balanced = balance_equation(equation)
    print(f"\nBalanced equation: {balanced}")
    
    print("\n" + "="*50)
    print("=== Example 5: Another Ionic Reaction ===")
    equation = "Fe+3 + I- -> Fe+2 + I2"
    balanced = balance_equation(equation)
    print(f"\nBalanced equation: {balanced}")
    
    print("\n" + "="*50)
    print("=== Example 6: Dichromate Reaction ===")
    equation = "Cu+2 + Cr2O7-2 + H+ -> Cu+3 + Cr+3 + H2O"
    balanced = balance_equation(equation)
    print(f"\nBalanced equation: {balanced}")
    
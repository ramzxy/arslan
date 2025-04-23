# Chemical Equation Balancer

A Python tool for balancing chemical equations, including those with ionic charges.

## About
This project was developed as part of coursework for the Twente Pathway College. It implements a mathematical approach to balance chemical equations using matrix operations and Gaussian elimination.

## Features

- Balances both regular and ionic chemical equations
- Handles charges in ionic compounds (Fe+3, SO4-2, etc.)
- Uses matrix-based approach with Gaussian elimination
- Provides step-by-step explanation of the balancing process
- Verifies balance by comparing element and charge counts

## Usage

Run the program directly:

```bash
python chemical_balancer.py
```

Or import the balance_equation function in your own code:

```python
from chemical_balancer import balance_equation

# Example: Balance a combustion reaction
equation = "C2H6 + O2 -> CO2 + H2O"
balanced_eq = balance_equation(equation)
print(balanced_eq)  # 2C2H6 + 7O2 -> 4CO2 + 6H2O

# Example: Balance an ionic reaction
equation = "Fe+3 + I- -> Fe+2 + I2"
balanced_eq = balance_equation(equation)
print(balanced_eq)  # 2Fe+3 + 2I- -> 2Fe+2 + I2
```

## Input Format

- Use standard chemical formulas with "+" separating compounds
- For ionic compounds, use the format Fe+3, SO4-2, etc.
- Separate reactants and products with "->"

Examples:
- C2H6 + O2 -> CO2 + H2O
- Fe+3 + I- -> Fe+2 + I2
- MnO4- + C2O4-2 + H+ -> Mn+2 + CO2 + H2O

## How It Works

1. Parses the chemical equation to identify compounds
2. Extracts elements and charges from each compound
3. Builds a matrix representing element/charge conservation
4. Solves the matrix using Gaussian elimination
5. Converts the solution to integer coefficients
6. Verifies balance by comparing element and charge counts on both sides

## Requirements

- Python 3.x
- NumPy
- re (standard library)
- fractions (standard library)

## Academic Context
This program is part of a chemistry project for the Twente Pathway College curriculum. It demonstrates the application of linear algebra to solve chemical balancing problems. 
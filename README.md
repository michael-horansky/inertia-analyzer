# inertia-analyzer
Analyzer that calculates tensors of inertia of discrete systems using symbolic language

# Usage
The program asks you to input a list of pointlike masses ("objects") in your system, each one as a line in the format [mass, displacement x component, displacement y component, displacement z component], four expressions separated by commas. After you input the last object, input an empty string by pressing ENTER and the program does the rest. You don't have to input the displacements from the center of mass; the program automatically calculates it and shifts the coordinates there, so you can pick displacements which are the simplest to input.

The expressions you input should be in the standard 'python-like' format, e.g. '**' as exponentiation, sqrt(...) as square root etc. Feel free to use any parametric symbols you like; the code uses the 'sympy' library to interpret the parameters and work with them symbolically, in a wolframalpha fashion. For more precise instructions refer to the sympy documentation.

# Outputs
The program calculates and outputs the following information (as of 22 Oct 2022):
1. The **total mass** of the system.
2. The **displacement of the center of mass** from the input origin.
3. The **displacements relative to the center of mass** of all objects.
4. The **tensor of inertia** for the center of mass.
5. The **principal axes** and their respective eigenvalue moments of inertia. If evaluable, the program also decides whether rotation about each principal axis is possible and stable.

All outputs use symbolic language, contain only the parameters mentioned in the inputs and are simplified through the standard sympy routine.

# Dependencies
Python 3, numpy, sympy

# Plans for future
Add arbitrary coordinate transformation with full Euler angles support; add a continuous density spatial function analyzer.

# Contact
michael.horansky@gmail.com

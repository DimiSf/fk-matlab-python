import numpy as np
import sympy as sp  # For symbolic math

# Define symbolic variables for theta values
th1, th2, th3, th4, th5, th6, th7 = sp.symbols('th1 th2 th3 th4 th5 th6 th7')

# DH parameter matrix using degrees
dh_matrix = [[0, 0.333, 0, th1],         # theta1
             [0, 0, -90, th2],          # theta2
             [0, 0.316, 90, th3],       # theta3
             [0.0825, 0, 90, th4],      # theta4
             [-0.0825, 0.384, -90, th5],# theta5
             [0, 0, 90, th6],           # theta6
             [0.088, 0, 90, th7]]        # theta7

# DH transformation function, using SymPy
def dh_transform(a, d, alpha, theta):
    # Convert angles from degrees to radians
    alpha_rad = sp.rad(alpha)  # Convert alpha from degrees to radians
    theta_rad = theta  # Keep theta as a symbol (it will be in degrees)

    A = sp.Matrix([
        [sp.cos(theta_rad), -sp.sin(theta_rad) * sp.cos(alpha_rad), sp.sin(theta_rad) * sp.sin(alpha_rad), a * sp.cos(theta_rad)],
        [sp.sin(theta_rad), sp.cos(theta_rad) * sp.cos(alpha_rad), -sp.cos(theta_rad) * sp.sin(alpha_rad), a * sp.sin(theta_rad)],
        [0, sp.sin(alpha_rad), sp.cos(alpha_rad), d],
        [0, 0, 0, 1]
    ])
    return A

# Compute individual transformation matrices
A_matrices = []
for i, (a, d, alpha, theta) in enumerate(dh_matrix):
    A_i = dh_transform(a, d, alpha, theta)
    A_matrices.append(A_i)

# Print transformation matrices
for i, A in enumerate(A_matrices):
    print(f"Transformation Matrix A_{i+1} (from frame {i+1} to frame {i}):")
    sp.pprint(A, use_unicode=True)
    print("\n")
print("---------------------------------------------------------\n")
T_matrices = []  # List to hold transformation matrices
T_matrices.append(A_matrices[0])  # T0_1 = A_1

for j in range(1, len(A_matrices)):
    T_new = T_matrices[-1] @ A_matrices[j]  # Multiply the last T matrix with the current A matrix
    T_matrices.append(T_new)  # Append the new transformation matrix

# Print the resulting transformation matrices
for i, T in enumerate(T_matrices):
    print(f"Transformation Matrix T_0_{i+1} (from frame {i+1} to base):")
    sp.pprint(T, use_unicode=True)
    print("\n")

T_endeff = A_matrices[0] @A_matrices[1] @A_matrices[2]@ A_matrices[3]@ A_matrices[4]@A_matrices[5]@ A_matrices[6]
sp.pprint(T_endeff, use_unicode=True)
if T_endeff == T_matrices[6]:
    print("\n")
    print("Correct! Check!")
else:
    print("Fail")
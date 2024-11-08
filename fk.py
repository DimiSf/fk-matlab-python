import numpy as np
import sympy as sp  # For symbolic math

# Define symbolic variables for theta values
th1, th2, th3, th4, th5, th6, th7 = sp.symbols('th1 th2 th3 th4 th5 th6 th7')

# DH parameter matrix using degrees
dh_matrix = [[0, 0.333, 0, th1],         # theta1
             [0, 0, -90, th2],           # theta2
             [0, 0.316, 90, th3],        # theta3
             [0.0825, 0, 90, th4],       # theta4
             [-0.0825, 0.384, -90, th5], # theta5
             [0, 0, 90, th6],            # theta6
             [0.088, 0, 90, th7]]        # theta7

# DH transformation function using Craig's method
def dh_transform(a, d, alpha, theta):
    # Convert alpha from degrees to radians for computation
    alpha_rad = sp.rad(alpha)  
    
    # Define each transformation matrix according to Craig's method
    Rx_alpha = sp.Matrix([
        [1, 0, 0, 0],
        [0, sp.cos(alpha_rad), -sp.sin(alpha_rad), 0],
        [0, sp.sin(alpha_rad), sp.cos(alpha_rad), 0],
        [0, 0, 0, 1]
    ])
    
    Dx_a = sp.Matrix([
        [1, 0, 0, a],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    
    Rz_theta = sp.Matrix([
        [sp.cos(theta), -sp.sin(theta), 0, 0],
        [sp.sin(theta), sp.cos(theta), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    
    Dz_d = sp.Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, d],
        [0, 0, 0, 1]
    ])
    
    # Combine transformations in the order specified by Craig's method
    A = Rx_alpha * Dx_a * Rz_theta * Dz_d
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

# Compute cumulative transformation matrices T_0_i
T_matrices = []  # List to hold cumulative transformation matrices
T_matrices.append(A_matrices[0])  # T0_1 = A_1

for j in range(1, len(A_matrices)):
    T_new = T_matrices[-1] @ A_matrices[j]  # Multiply the last T matrix with the current A matrix
    T_matrices.append(T_new)  # Append the new transformation matrix

# Print the resulting cumulative transformation matrices
for i, T in enumerate(T_matrices):
    print(f"Cumulative Transformation Matrix T_0_{i+1} (from base to frame {i+1}):")
    sp.pprint(T, use_unicode=True)
    print("\n")

# Compute the end-effector transformation matrix directly by chaining all A matrices
T_endeff = A_matrices[0] @ A_matrices[1] @ A_matrices[2] @ A_matrices[3] @ A_matrices[4] @ A_matrices[5] @ A_matrices[6]
print("End-Effector Transformation Matrix T_0_7:")
sp.pprint(T_endeff, use_unicode=True)

# Verify if the end-effector transformation matches the cumulative transformation T_0_7
if T_endeff == T_matrices[6]:
    print("\nCorrect! The transformation matches.")
else:
    print("Fail! The transformation does not match.")

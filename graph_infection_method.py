import numpy as np
import math
from utils import (
    adj_mat_ground_truth,
    angle_dot,
    get_abs_angle,
)

# example arbitrary non-negative real matrix A
# power itertaion will NOT converge for this A
A = np.array([
    [0, 0, 0, 1],
    [1, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 2, 1, 0]])

######################################
# ground truth of the matrix or graph
######################################
_, _, true_dom_eigVec = adj_mat_ground_truth(A, input_is_matrix=True)


"""
Our Method:
"""
Delta_t = 1
num_steps = 10

# initial x and I
x_0 = np.ones((A.shape[0], 1))
I_0 = np.sum(x_0)

x_old = x_0
I_old = I_0

# whether to normalize to unit vector
normalize = True

for i in range(1, 1+num_steps):
    x_new = x_old + A @ x_old * Delta_t
    I_new = np.sum(x_new)

    m = (math.log(I_new) - math.log(I_old)) / Delta_t
    dominant_eigenval = (math.exp(m * Delta_t) - 1) / Delta_t

    if normalize:
        # normalize by unit vector
        x_new /= np.sqrt(np.sum(x_new ** 2))
        I_new = np.sum(x_new)

    x_old = x_new
    I_old = I_new

    _, angle = angle_dot(true_dom_eigVec, x_old)
    abs_angle = get_abs_angle(angle)

print(f"dominant_eigenvalue = {dominant_eigenval}")
print(f"dominant eigenvector:\n{x_old}")
print(f"dot_angles = {abs_angle}")

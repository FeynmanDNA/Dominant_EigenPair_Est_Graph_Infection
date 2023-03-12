Const = list(
    N_STEP = 10,
    START = 0.025,
    END = Inf,
    METHOD = "max",
    PRECISION = 1e-6,
    MINITER = 20,
    MAXITER = 1000
);

EigenResult = function(n_step, vector, value) {
    return (list(
        n_step = n_step,
        vector = vector,
        value = value
    ));
}

ExponentialResult = function(n_step, vector, value, I_seq) {
    return (list(
        n_step = n_step,
        vector = vector,
        value = value,
        I_seq = I_seq
    ));
}

# Classical power method for estimation of dominant eigenvalue and eigenvector.
# u: Estimated dominant eigenvalue to be updated each round.
# b: Estimated dominant eigenvector to be updated each round.
# The algorithm stops until the required precision is obtained or until the number of iterations exceeds the maximum.
# Setting precision = None makes the algorithm proceed until the maximum number of iterations.
power_method = function(
    A, miniter = Const$MINITER, maxiter = Const$MAXITER, precision = Const$PRECISION
) {

    n_node = nrow(A);
    # degrees = A %*% rep(1, n_node);
    # vector = degrees / sum(degrees);
    vector = rep(1, n_node);

    value = (t(vector) %*% A %*% vector / (t(vector) %*% vector));
    for (i in 1:maxiter) {
        new_vector = A %*% vector;
        new_value = (t(vector) %*% A %*% vector) / (t(vector) %*% vector);
        if (i > miniter & abs((new_value - value) / new_value) < precision) {
            print("== POWER METHOD ==");
            print(sprintf("# Algorithm successfully converged over %s iterations.", i));
            print(sprintf("# Dominant eigenvalue = %s.", new_value));
            return (EigenResult(n_step = i, vector = new_vector, value = new_value));
        } else {
            vector = new_vector / (sum(new_vector ** 2) ** 0.5);
            value = new_value;
        }
    }
    print("== POWER METHOD ==");
    print(sprintf("# Algorithm failed to converge after %s iterations.", maxiter));
    print(sprintf("# Last estimated dominant eigenvalue = {value}.", value));
    return (EigenResult(n_step = maxiter, vector = vector, value = value));

}

# === Pseudocode ===
# Loop until convergence or other stop criteria:
## x[n + 1] = (A + identity) @ x[n] = x[n] + A @ x[n]
## I[n + 1] = sum(x[n + 1])
## slope[n + 1] = ln(I[n + 1]) - ln(I[n])
## x[n + 1] = x[n + 1] / norm(x[n + 1])
## I[n + 1] = I[n + 1] / norm(I[n + 1])
# End loop.
# corrected_final_slope = exp(final_slope) - 1
# Return corrected_final_slope, final_x and step no. to convergence.
# ======
# One major problem with regular exponential epidemic models is that for large graphs, the total infection severity
# grows quickly and easily exceeds the maximum number range for most computers. To remedy this, a normalization step
# is included so that x[n + 1] will always be scaled down by a factor of its norm at each step, so that explosive
# growth does not occur. Note that this will not affect the estimation of slopes, because within each step, x[n + 1]
# and x[n] will be scaled down by the same factor, as well as I[n + 1] and I[n], which will not affect computation
# in the context of logarithmic values.
# Nonetheless, normalization will destroy the evolution of severity as computed from the original exponential model. Therefore,
# an option is included to disable the normalization step so that other details of the evolution of severity can be extracted
# and studied. However, one must note that under this mode, large graphs are liable to produce overflow errors.
exponential_method = function(
    A, miniter = Const$MINITER, maxiter = Const$MAXITER, precision = Const$PRECISION,
    normalize = TRUE, verbose = TRUE
) {

    n_node = nrow(A);
    # degrees = A %*% rep(1, n_node);
    #densest_nodes = [i for i in range(n_node) if degrees[i] == np.max(degrees[i])]

    # x = degrees / sum(degrees);
    x = rep(1, n_node)
    # x = np.repeat(1 / n_node, n_node)
    # Partitions initial severity of 1 uniformly over nodes with the densest connections.
    # x = np.zeros(shape=(n_node,))
    # x[np.ix_(densest_nodes)] = 1 / len(densest_nodes)
    I = sum(x);
    if (!normalize) {
        I_seq = c(I);
    }
    # Initializing with slope = 0 is feasible, since the Perron-Frobenius theorem guarantees a positive eigenvalue
    # for the adjacency matrix, with the fact that total severity strictly increases over time.
    # slope = 0
    estimate = 0;
    transition_matrix = A + diag(n_node);

    for (i in 1:maxiter) {
        # x_new = x + A @ x * alpha
        x_new = transition_matrix %*% x;
        I_new = sum(x_new);
        estimate_new = I_new / I - 1;
        # slope_new = (math.log(I_new) - math.log(I)) / alpha
        # if precision is not None and i > miniter and abs((slope_new - slope) / slope_new) < precision:
        if (i > miniter & abs((estimate_new - estimate) / estimate_new) < precision) {
            # corrected_slope_new = apply_correction(slope_new)
            if (verbose) {
                print("== EPIDEMIC-BASED EXPONENTIAL METHOD ==");
                print(sprintf("# Algorithm successfully converged over %s iterations.", i));
                # print(f"# Dominant eigenvalue = {corrected_slope_new}.")
                print(sprintf("# Dominant eigenvalue = %s.", estimate_new));
            }
            if (normalize) {
                # return EigenResult(n_step = i, vector = x_new, value = corrected_slope_new)
                return (EigenResult(n_step = i, vector = x_new, value = estimate_new));
            } else {
                I_seq = c(I_seq, I_new);
                #I_seq.append(I_new);
                # return ExponentialResult(n_step = i, vector = x_new, value = corrected_slope_new,
                #     I_seq = np.array(I_seq))
                return (ExponentialResult(n_step = i, vector = x_new, value = estimate_new,
                    I_seq = I_seq));
            }
        } else if (normalize) {
            x = x_new / (sum(x_new ** 2) ** 0.5);
            I = sum(x);
            # slope = slope_new
            estimate = estimate_new;
        } else {
            x = x_new;
            I = I_new;
            I_seq = c(I_seq, I_new);
            #I_seq.append(I_new)
            # slope = slope_new
            estimate = estimate_new;
        }
    }
    # corrected_slope = apply_correction(slope)
    if (verbose) {
        print("== EPIDEMIC-BASED EXPONENTIAL METHOD ==")
        print(sprintf("# Algorithm failed to converge after %s iterations.", maxiter));
        # print(f"# Last estimated dominant eigenvalue = {math.sqrt(corrected_slope}.")
        print(sprintf("# Last estimated dominant eigenvalue = %s.", estimate));
    }
    if (normalize) {
        # return EigenResult(n_step = maxiter, vector = x, value = corrected_slope)
        return (EigenResult(n_step = maxiter, vector = x, value = estimate));
    } else {
        # return ExponentialResult(n_step = maxiter, vector = x, value = corrected_slope,
        #     I_seq = np.array(I_seq))
        return (ExponentialResult(n_step = maxiter, vector = x, value = estimate,
            I_seq = I_seq));
    }

}

print("Module successfully run.")

matrix_0 = rbind(
    c(0, 1, 1, 1, 1, 0),
    c(1, 0, 0, 0, 0, 1),
    c(1, 0, 0, 0, 0, 1),
    c(1, 0, 0, 0, 0, 1),
    c(1, 0, 0, 0, 0, 1),
    c(0, 1, 1, 1, 1, 0)
)
print("MATRIX 0")
power_method(matrix_0)
exponential_method(matrix_0)

matrix_1 = rbind(
    c(0.  , 0.  , 0.  , 1.  , 0.5 , 0.  , 0.  ),
    c(1.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.25),
    c(0.  , 1.  , 0.  , 0.  , 0.  , 0.  , 0.25),
    c(0.  , 0.  , 1.  , 0.  , 0.  , 0.  , 0.  ),
    c(0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.25),
    c(0.  , 0.  , 0.  , 0.  , 0.5 , 0.  , 0.25),
    c(0.  , 0.  , 0.  , 0.  , 0.  , 1.  , 0.  )
)
print("MATRIX 1")
power_method(matrix_1)
exponential_method(matrix_1)

matrix_2 = rbind(
    c(0, 1, 0, 0, 0, 0, 0, 0, 0),
    c(1, 0, 1, 1, 0, 0, 0, 0, 0),
    c(0, 1, 0, 0, 0, 0, 0, 0, 0),
    c(0, 1, 0, 0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 1, 0, 1, 1, 1, 1),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 0, 0, 0, 0)
)
print("MATRIX 2")
power_method(matrix_2)
exponential_method(matrix_2)

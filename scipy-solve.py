from scipy.optimize import linprog

obj = [-20, -12, -40, -25]

lhs_ineq = [[1, 1, 1, 1],
            [3, 2, 1, 0],
            [0, 1 , 2, 3]]  # Yellow constraint left side

rhs_ineq = [50,  # Red constraint right side
            100,  # Blue constraint right side
            90]  # Yellow constraint right side

bnd = [(0, float("inf")),
       (0, float("inf")),
       (0, float("inf")),
       (0, float("inf"))]

opt = linprog(c=obj, A_ub=lhs_ineq, b_ub=rhs_ineq, bounds=bnd, method="revised simplex")

print(opt)


from pulp import LpMaximize, LpProblem, LpStatus, lpSum, LpVariable

model = LpProblem(name="test", sense=LpMaximize)
x = LpVariable(name="x", lowBound=0)
y = LpVariable(name="y", lowBound=0)

model += (2 * x + y <= 20, "red_constraint")
model += (4 * x - 5 * y >= -10, "blue_constraint")
model += (-x + 2 * y >= -2, "yellow_constraint")
model += (-x + 5 * y == 15, "green_constraint")

model += lpSum([x, 2 * y])
print(model)
status = model.solve()
print(model.variables()[0].value())


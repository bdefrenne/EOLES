from pyomo.environ import *

model = ConcreteModel()

# declare decision variables
model.x = Var(domain=NonNegativeReals)

# declare objective
model.profit = Objective(
    expr=40 * model.x,
    sense=maximize)

# declare constraints
model.demand = Constraint(expr=model.x <= 40)
model.laborA = Constraint(expr=model.x <= 80)
model.laborB = Constraint(expr=2 * model.x <= 100)

# solve
SolverFactory('cbc').solve(model).write()
print("Profit = ", model.profit(), " per week")
print("X = ", model.x(), " units per week")

variable_dims = {
    'capacity': ['tec_total'],
    'charge_capacity': ['storing'],
    'volume': ['storing'],
    'generation_nonvre': ['nonvre', 'time'],
    'storage': ['storing', 'time'],
    'stored': ['storing', 'time'],
    'reserve': ['frr', 'time']}

#: Constraint definitions as dictionnary from constraint name
#: to a tuple with a list of sets and a rule function.
constraint_definitions = {
    'c_adequacy': (['time'], self.adequacy_rule),
    'c_max_generation_nonvre':
        (['nonvre', 'time'], self.max_generation_nonvre_rule),
    'c_storing':
        (['storing', 'time_without_last'], self.storing_rule),
    'c_storage_refilled': (['storing'], self.storage_refilled_rule),
    'c_max_stored': (['storing', 'time'], self.max_stored_rule),
    'c_max_charge_capacity':
        (['storing'], self.max_charge_capacity_rule),
    'c_max_storage': (['storing', 'time'], self.max_storage_rule),
    'c_max_yearly_generation':
        (['biogas'], self.max_yearly_generation_rule),
    'c_max_monthly_generation_lake':
        (['months'], self.max_monthly_generation_lake_rule),
    'c_max_generation_frr':
        (['frr', 'time'], self.max_generation_frr_rule),
    'c_reserve': (['time'], self.reserve_rule)}

"""Solve optimization problem."""
# Initialize model
model = ConcreteModel()

# Cost
model.cost = Objective(rule=cost_rule(model))

set_names = ['technology', 'storing', 'vre', 'nonvre', 'frr',
             'time', 'time_without_last', 'months', 'biogas']

def cost_rule(model):
    """Get cost (G€).

    :param m: Model.
    :type m: :py:class:`pyomo.Model`

    :returns: Expression.
    :rtype: :py:class:`pyo.Expression`
    """
    # CAPEX
    expr = sum(get_cost_capex_expr(model, tec)
               for tec in model.technology)

    # fixed OM costs
    expr += sum(get_cost_om_fixed_expr(model, tec)
                for tec in model.technology)

    # Storage volume costs
    expr += sum(get_cost_volume_capex_expr(model, tec)
                for tec in model.storing)

    # Charging capacity costs
    expr += sum(get_cost_charge_capacity_expr(model, tec)
                for tec in model.storing)

    # Generation costs
    expr += sum(get_cost_om_variable_full_expr(model, tec)
                for tec in model.nonvre)

    # Convert
    expr /= 1.e3

    return expr


def update_sets_coords(model):
    """Update sets, coordinates with loaded input."""
    # Update sets and coordinates with loaded input (for time index)
    for set_name in set_names:
        initialize = getattr(self.input, set_name)
        self._add_set_coord(set_name, ordered=True, initialize=initialize)

def add_set_coord(model, set_name, ordered=True, initialize=None):
    """Add set and coordinates to model.

    :param set_name: Set name.
    :param ordered: Whether set is ordered.
    :param initialize: Initialization array.
    :type set_name: str
    :type ordered: bool
    :type initialize: sequence
    """
    # Add set
    s = Set(ordered=ordered, initialize=initialize)
    setattr(model, set_name, s)

# Add constraints
for cstr_name, cstr_def in constraint_definitions.items():
    indices = [getattr(model, dim) for dim in cstr_def[0]]
    constraint = pyo.Constraint(*indices, rule=cstr_def[1])
    setattr(self.model, cstr_name, constraint)

# Create `dual` and `rc` suffixes component on the instance
# so the solver plugin will collect duals and reduced costs
self.model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)
self.model.rc = pyo.Suffix(direction=pyo.Suffix.IMPORT)

# Create solver
options = None
cfg_solver = self.cfg.get('solver')
if cfg_solver is not None:
    solver_name = cfg_solver.get('name') or 'cbc'
options = cfg_solver.get('options')
self.solver = pyo.SolverFactory(solver_name)

# Set solver options
if options is not None:
    for
opt, val in options.items():
self.solver.options[opt] = val
tee = False if self.cfg.get('no_verbose') else True

# Solve
self.info('Solving optimization problem')
self.results = self.solver.solve(self.model, tee=tee)

# Check
if ((self.results.solver.status == pyo.SolverStatus.ok) and
        (self.results.solver.termination_condition ==
         pyo.TerminationCondition.optimal)):
    self.info("Solution feasible and optimal.")
elif (self.results.solver.termination_condition ==
      pyo.TerminationCondition.infeasible):
    self.warning("Solution infeasible!")
# raise RuntimeError(str(self.results.solver))
else:
# something else is wrong
    raise RuntimeError(str(self.results.solver))

# Parse solution
ds = support.parse_solution(
    self.model, self.variable_dims, self.cfg, **kwargs)

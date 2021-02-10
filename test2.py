from pyomo.environ import *
import pandas as pd
from pyomo.repn.plugins.baron_writer import NonNegativeReals

""" DATA """
production_profile = pd.read_csv("inputs/vre_profiles2006.csv", index_col=[0, 1], squeeze=True, header=None)
existing_capacities = pd.read_csv("inputs/existing_capas.csv", index_col=[0, 1], squeeze=True, header=None)
maximum_capacities = pd.read_csv("inputs/max_capas.csv", index_col=[0, 1], squeeze=True, header=None)
capex = pd.read_csv("inputs/annuities.csv", index_col=[0, 1], squeeze=True, header=None)
capex_storage = pd.read_csv("inputs/str_annuities.csv", index_col=[0, 1], squeeze=True, header=None)
demand = pd.read_csv("inputs/demand2006.csv", index_col=[0, 1], squeeze=True, header=None)
lake_inflows = pd.read_csv("inputs/lake_inflows.csv", index_col=[0, 1], squeeze=True, header=None)
gene_river = pd.read_csv("inputs/run_of_river.csv", index_col=[0, 1], squeeze=True, header=None)
reserve_requirements = pd.read_csv("inputs/reserve_requirements.csv", index_col=[0, 1], squeeze=True, header=None)
fOM = pd.read_csv("inputs/fO&M.csv", index_col=[0, 1], squeeze=True, header=None)
vOM = pd.read_csv("inputs/vO&M.csv", index_col=[0, 1], squeeze=True, header=None)

model = ConcreteModel()

set_names = ['technology', 'storing', 'vre', 'nonvre', 'frr',
                   'time', 'time_without_last', 'months', 'biogas']

# Define the sets
first_hour = 0
last_hour = 8759
model.time = RangeSet(first_hour, last_hour)
model.technologies = Set(initialize=["pv", "onshore", "offshore", "river", "lake", "biogas", "gas", "phs", "battery", "methanation"])
model.generations = Set(initialize=["pv", "onshore", "offshore", "river", "lake", "biogas"])
model.vre = Set(initialize=["pv", "onshore", "offshore"])
model.ncomb = Set(initialize=["pv", "onshore", "offshore", "river", "lake", "phs", "battery"])
model.comb = Set(initialize=["biogas", "methanation"])
model.storages = Set(initialize=["phs", "battery", "methanation"])
model.frr = Set(initialize=["lake", "phs", "battery", "gas"])

# Define variables
model.capacities = Var(model.technologies, within=NonNegativeReals, initialize=0)
model.generation = Var(((gen, hour) for gen in model.generations for hour in model.time), within=NonNegativeReals, initialize=0)
model.storage_power_capacity = Var(model.storages, within=NonNegativeReals, initialize=0)
model.storage_energy_capacity = Var(model.storages, within=NonNegativeReals, initialize=0)
model.storage = Var(((storage, hour) for storage in model.storages for hour in model.time), within=NonNegativeReals, initialize=0)
model.reserve = Var(((reserve, hour) for reserve in model.frr for hour in model.time), within=NonNegativeReals, initialize=0)


def generation_vre_constraint_rule(model, time, tec):
    return model.generation[tec, t] = model.capacities[tec] * production_profile[tec][t]

model.generation_vre_constraint = Constraint(model.time, model.vre, rule=generation_vre_constraint_rule)

pyo.Constraint(*indices, rule=cstr_def[1])

gene_vre(vre,h)..                GENE(vre,h)             =e=     CAPA(vre)*load_factor(vre,h);
gene_capa(tec,h)..               CAPA(tec)               =g=     GENE(tec,h);
combustion(h)..                  GENE('gas',h)           =e=     sum(comb,GENE(comb,h));
capa_frr(frr,h)..                CAPA(frr)               =g=     GENE(frr,h) + RSV(frr,h);
storing(h,h+1,str)..             STORED(str,h+1)         =e=     STORED(str,h) + STORAGE(str,h)*eta_in(str) - GENE(str,h)/eta_out(str);
storage_const(str,first,last)..  STORED(str,first)       =e=     STORED(str,last)+STORAGE(str,last)*eta_in(str)-GENE(str,last)/eta_out(str);
lake_res(m)..                    lake_inflows(m)         =g=     sum(h$(month(h) = ord(m)),GENE('lake',h))/1000;
stored_cap(str,h)..              STORED(str,h)           =l=     CAPACITY(str);
storage_capa1(str,h)..           S(str)                  =g=     STORAGE(str,h);
storage_capa2(str)..             S(str)                  =l=     CAPA(str);
biogas_const..                   sum(h,GENE('biogas',h)) =l=     max_biogas*1000;
reserves(h)..                    sum(frr, RSV(frr,h))    =e=     sum(vre,epsilon(vre)*CAPA(vre))+ demand(h)*load_uncertainty*(1+delta);
adequacy(h)..                    sum(ncomb,GENE(ncomb,h))+GENE('gas',h)    =g=     demand(h)+sum(str,STORAGE(str,h)) ;
obj..                            COST                    =e=     (sum(tec,(CAPA(tec)-capa_ex(tec))*capex(tec))+ sum(str,CAPACITY(str)*capex_en(str))+sum(tec,(CAPA(tec)*fOM(tec)))+ sum(str,S(str)*(s_capex(str)+s_opex(str))) + sum((tec,h),GENE(tec,h)*vOM(tec)))/1000;





variable_names = [
    'capacityfactor', 'demand', 'inflows_lake',
    'reserve_requirement', 'capacity_existing', 'capacity_max',
    'capex', 'volume_capex', 'om_fixed', 'om_variable',
    'efficiency_charging', 'efficiency_discharging',
    'max_generation_year', 'uncertainty', 'variation_factor',
    'charge_capex', 'charge_om_fixed'
]

variable_dimensions = {
    'capacity': ['tec_total'],
    'charge_capacity': ['storing'],
    'volume': ['storing'],
    'generation_nonvre': ['nonvre', 'time'],
    'storage': ['storing', 'time'],
    'stored': ['storing', 'time'],
    'reserve': ['frr', 'time']}

constraint_definitions = {
    'c_adequacy': (['time'], adequacy_rule),
    'c_max_generation_nonvre':
        (['nonvre', 'time'],  max_generation_nonvre_rule),
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


update_sets_coords()

# Add variables to model
add_state_variables()
add_constraints()


def update_sets_coords():
    """Update sets, coordinates with loaded input."""
    # Update sets and coordinates with loaded input (for time index)
    for set_name in self._set_names:
        initialize = getattr(self.input, set_name)
        add_set_coord(set_name, initialize=initialize)


def add_set_coord(set_name, initialize=None):
    # Add set
    s = Set(ordered=True, initialize=initialize)
    setattr(model, set_name, s)




s = Set(ordered=True, initialize=initialize)



model.time = Set(ordered=True, initialize=range(0, 8759))
model.adequacy_rule = Constraint(*[getattr(model, dim) for dim in cstr_def[0]], rule=cstr_def[1])

def add_state_variables():
    """Add state variable to model."""
    for variable_name, dimensions in variable_dimensions.items():
        indices = [getattr(model, dimension) for dimension in dimensions]
        variable_value = Var(*indices, domain=NonNegativeReals)
        setattr(model, variable_name, variable_value)

for cstr_name, cstr_def in self.constraint_definitions.items():
    indices = [getattr(self.model, dim) for dim in cstr_def[0]]
    constraint = Constraint(*indices, rule=cstr_def[1])
    setattr(self.model, cstr_name, constraint)





def cost_rule(model):
    """ Get cost (G€) """
    expr = sum(get_new_capacity_capex_expression(model, tec)
               for tec in model.technology)

    expr += sum(get_fixed_om_cost_expression(model, tec)
                for tec in model.technology)

    expr += sum(get_storage_volume_cost_expression(model, tec)
                for tec in model.storing)

    expr += sum(get_storage_capacity_cost_expression(model, tec)
                for tec in model.storing)

    expr += sum(get_total_generation_cost_expression(model, tec)
                for tec in model.nonvre)

    return expr / 1.e3


def adequacy_rule(model, t):
    """Get constraint for adequacy."""
    # Get total generation
    generation = get_generation_expression(model, t)
    stored = + sum(model.storage[storing, t] for storing in model.storing)

    return generation >= demand[t] + stored

def get_generation_expression(model, t):
    return sum(production_profile[tec][t] * model.capacity[tec] for tec in model.vre)\
    + sum(model.generation_nonvre[tec, t] for tec in model.nonvre)

def get_new_capacity_capex_expression(model, tec):
    """Get CAPEX cost of additional capacity for technology.

    :param m: Model.
    :param tec: Technology.
    :type m: :py:class:`pyomo.Model`
    :type tec: str

    :returns: Expression.
    :rtype: :py:class:`Expression`
    """
    # Get new capacity
    total_capacity = model.capacity[tec]
    new_capacity = total_capacity - existing_capacities[tec]

    # Return product
    return capex[tec] * new_capacity


def get_fixed_om_cost_expression(model, tec):
    """Get fixed O&M cost for technology.

    :param m: Model.
    :param tec: Technology.
    :type m: :py:class:`pyomo.Model`
    :type tec: str

    :returns: Expression.
    :rtype: :py:class:`Expression`
    """
    # Get total capacity for technology
    total_capacity = model.capacity[tec]

    # Return product
    return fOM[tec] * total_capacity


def get_total_generation_cost_expression(model, tec):
    """Get generation cost over full time series for technology.

    :param m: Model.
    :param tec: Technology.
    :type m: :py:class:`pyomo.Model`
    :type tec: str

    :returns: Expression.
    :rtype: :py:class:`Expression`
    """
    return sum(get_generation_cost_expression(model, tec, t) for t in model.time)


def get_generation_cost_expression(model, tec, t):
    """Get generation cost at timestamp for technology.

    :param m: Model.
    :param tec: Technology.
    :param t: Timestamp.
    :type m: :py:class:`pyomo.Model`
    :type tec: str
    :type t: :py:class:`pyomo.Set`

    :returns: Expression.
    :rtype: :py:class:`Expression`
    """
    return model.generation_nonvre[tec, t] * vOM['om_variable'][tec]


def get_storage_volume_cost_expression(model, tec):
    """Get storage-volume cost for technology.

    :param m: Model.
    :param tec: Technology.
    :type m: :py:class:`pyomo.Model`
    :type tec: str

    :returns: Expression.
    :rtype: :py:class:`Expression`
    """
    return capex_storage[tec] * model.volume[tec]


def get_storage_capacity_cost_expression(model, tec):
    """Get charging-capacity cost for technology.

    :param m: Model.
    :param tec: Technology.
    :type m: :py:class:`pyomo.Model`
    :type tec: str

    :returns: Expression.
    :rtype: :py:class:`Expression`
    """
    charge_capacity_cost_tec = (capex[tec] +
                                fOM[tec])

    return model.charge_capacity[tec] * charge_capacity_cost_tec






def _update_state_variables(self, **kwargs):
    """Update variables with loaded input."""
    for variable_name, dims in self.variable_dims.items():
        self._add_state_variable(variable_name, dims, **kwargs)

def _add_state_variable(self, variable_name, dims, **kwargs):
    """Add state variable to model.

    :param variable_name: Variable name.
    :param dims: Variable dimensions.
    :type variable_name: str
    :type variable_dimensions: sequence of :py:class:`str`
    """
    indices = [getattr(self.model, dim) for dim in dims]
    var = Var(*indices, domain=NonNegativeReals)
    setattr(self.model, variable_name, var)


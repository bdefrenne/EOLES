from pathlib import Path
from collections import OrderedDict
import numpy as np
import pandas as pd
import xarray as xr
import pyomo.environ as pyo


class Optimizer(OptimizerBase):
    """EOLES e4clim implementation."""

    def __init__(self, app, cfg=None, **kwargs):
        """Build optimizer context.

        :param app: Context application.
        :param cfg: Optimizer configuration.
        """
        # Initialize as OptimizerBase
        name = 'eoles'
        super(Optimizer, self).__init__(app, name, cfg=cfg, **kwargs)

        #: Demand component manager.
        self.context_component_demand = None
        #: Generation component managers.
        self.context_component_generation = None
        #: Generation component names.
        self.component_generation_names = None
        # Add available components
        self._add_available_components(**kwargs)

        #: Whether problem is regional
        self.regional = bool(self.input.regional)

        #: Model.
        self.model = None
        #: Solver.
        self.solver = None
        #: Results.
        self.results = None
        #: Coordinates.
        self.coords = None

        #: Sets.
        self._set_names = ['technology', 'storing', 'vre', 'nonvre', 'frr',
                           'time', 'time_without_last', 'months', 'biogas']

        #: Coordinates.
        self.coords = OrderedDict({})

        #: Variable dimensions.
        self.variable_dims = {
            'capacity': ['tec_total'],
            'charge_capacity': ['storing'],
            'volume': ['storing'],
            'generation_nonvre': ['nonvre', 'time'],
            'storage': ['storing', 'time'],
            'stored': ['storing', 'time'],
            'reserve': ['frr', 'time']}

        #: Constraint definitions as dictionnary from constraint name
        #: to a tuple with a list of sets and a rule function.
        self.constraint_definitions = {
            'c_adequacy': (['time'], self.adequacy_rule),
            'c_max_generation_nonvre':
            (['nonvre', 'time'],  self.max_generation_nonvre_rule),
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

        if self.regional:
            # Add  regional sets
            self._set_names += ['region', 'vre_total',
                                'vre_regional', 'tec_total']

            # Add/update variable dimensions
            self.variable_dims.update({
                'capacity_regional': ['vre_regional', 'region']})

            # Add/update constraints
            self.constraint_definitions.update({
                'c_max_capacity': (['vre_regional'], self.max_capacity_rule)})

        # Update variable names (as groups)
        support.update_solution_variables(
            self.solution, self.variable_dims, self.cfg, **kwargs)

    def get_new_input(self, **kwargs) -> None:
        """Get new input.

        :param returns: New input.
        """
        return Input(self, **kwargs)

    def get_new_solution(self, **kwargs) -> None:
        """Get new solution.

        :param returns: New solution.
        """
        return Solution(self, **kwargs)

    def solve(self, *args, **kwargs):
        """Solve optimization problem."""
        # Initialize model
        self.model = pyo.ConcreteModel()

        # Update sets, coordinates and variables with loaded input
        self._update_sets_coords(**kwargs)
        # Add variables to model
        self._update_state_variables(**kwargs)

        # Fix variables
        self._fix_variables()

        # Bound variables
        self._bound_variables()

        # Cost
        self.info('Defining cost')
        self.model.cost = pyo.Objective(rule=self.cost_rule)

        # Add constraints
        self.info('Defining constraints')
        for cstr_name, cstr_def in self.constraint_definitions.items():
            indices = [getattr(self.model, dim) for dim in cstr_def[0]]
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
            for opt, val in options.items():
                self.solver.options[opt] = val
        tee = False if self.cfg.get('no_verbose') else True
                                                                                                        """ INFO ?? """
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

        return ds

    def max_generation_nonvre_rule(self, m, tec, t):
        """Get constraint on maximum power for non-VRE technologies.

        :param m: Model.
        :param tec: Technology.
        :param t: Timestamp.
        :type m: :py:class:`pyomo.Model`
        :type tec: str
        :type t: :py:class:`pandas.Timestamp`

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return m.generation_nonvre[tec, t] <= m.capacity[tec]

    def storing_rule(self, m, tec, t):
        """Get constraint on storing.

        :param m: Model.
        :param tec: Technology.
        :param t: Timestamp.
        :type m: :py:class:`pyomo.Model`
        :type tec: str
        :type t: :py:class:`pandas.Timestamp`

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        charge = m.storage[tec, t] * self.input['efficiency_charging'][tec]
        discharge = m.generation_nonvre[tec, t] / self.input[
            'efficiency_discharging'][tec]
        flux = charge - discharge
        tp1 = self.input.next_time[t]

        return m.stored[tec, tp1] == m.stored[tec, t] + flux

    def storage_refilled_rule(self, m, tec):
        """Get constraint on stored energy to be higher at the end
        than at the start.

        :param m: Model.
        :param tec: Technology.
        :type m: :py:class:`pyomo.Model`
        :type tec: str

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return m.stored[tec, m.time[-1]] >= m.stored[tec, m.time[1]]

    def max_stored_rule(self, m, tec, t):
        """Get constraint on maximum stored volume.

        :param m: Model.
        :param tec: Technology.
        :param t: Timestamp.
        :type m: :py:class:`pyomo.Model`
        :type tec: str
        :type t: :py:class:`pandas.Timestamp`

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return m.stored[tec, t] <= m.volume[tec]

    def max_charge_capacity_rule(self, m, tec):
        """Get constraint on charging capacity.

        :param m: Model.
        :param tec: Technology.
        :type m: :py:class:`pyomo.Model`
        :type tec: str

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return m.charge_capacity[tec] <= m.capacity[tec]

    def max_storage_rule(self, m, tec, t):
        """Get constraint on storage power.

        :param m: Model.
        :param tec: Technology.
        :param t: Timestamp.
        :type m: :py:class:`pyomo.Model`
        :type tec: str
        :type t: :py:class:`pandas.Timestamp`

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return m.storage[tec, t] <= m.charge_capacity[tec]

    def max_yearly_generation_rule(self, m, tec):
        """Get constraint on maximum yearly biogas generation.

        :param m: Model.
        :param tec: Technology.
        :type m: :py:class:`pyomo.Model`
        :type tec: str

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        fact = 24 * 365 / len(m.time)
        gen_year_biogas = sum(m.generation_nonvre[tec, t]
                              for t in m.time) * fact

        return gen_year_biogas <= self.input['max_generation_year'][tec]

    def max_monthly_generation_lake_rule(self, m, month):
        """Get constraint on maximum monthly lake generation.

        :param m: Model.
        :param month: Month
        :type m: :py:class:`pyomo.Model`
        :type month: int

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        tec = 'lake'
        inflows = self.input['inflows_{}'.format(tec)][month]
        gen_lake_month = sum(m.generation_nonvre[tec, t]
                             for t in self.input.month_timestamps[month])
        return gen_lake_month <= inflows * 1000

    def max_generation_frr_rule(self, m, tec, t):
        """Get constraint on maximum generation including reserves
        for FRR technologies.

        :param m: Model.
        :param tec: Technology.
        :param t: Timestamp.
        :type m: :py:class:`pyomo.Model`
        :type tec: str
        :type t: :py:class:`pandas.Timestamp`

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return (m.generation_nonvre[tec, t] + m.reserve[tec, t] <=
                m.capacity[tec])

    def reserve_rule(self, m, t):
        """Get constraint for reserves.

        :param m: Model.
        :param t: Timestamp.
        :type m: :py:class:`pyomo.Model`
        :type t: :py:class:`pandas.Timestamp`

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        res_req = sum(self.input['reserve_requirement'][tec] *
                      _get_capacity_total(m, tec) for tec in m.vre)
        load_req = (self.input['demand'][t] * self.input['uncertainty'][
            'demand'] * (1 + self.input['variation_factor']['demand']))

        return sum(m.reserve[tec, t] for tec in m.frr) == res_req + load_req

    def max_capacity_rule(self, m, tec):
        """Get constraint on capacity.

        :param m: Model.
        :param tec: Technology.
        :type m: :py:class:`pyomo.Model`
        :type tec: str

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        # Get total capacity for technology
        cap = _get_capacity_total(m, tec)

        return cap <= self.input['capacity_max'][tec]

    def _get_cost_capex_expr(self, m, tec):
        """Get CAPEX cost of additional capacity for technology.

        :param m: Model.
        :param tec: Technology.
        :type m: :py:class:`pyomo.Model`
        :type tec: str

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        # Get new capacity
        cap = _get_capacity_total(m, tec)
        capacity_new = cap - self.input['capacity_existing'][tec]

        # Return product
        return self.input['capex'][tec] * capacity_new

    def _get_cost_om_fixed_expr(self, m, tec):
        """Get fixed O&M cost for technology.

        :param m: Model.
        :param tec: Technology.
        :type m: :py:class:`pyomo.Model`
        :type tec: str

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        # Get total capacity for technology
        cap = _get_capacity_total(m, tec)

        # Return product
        return self.input['om_fixed'][tec] * cap

    def _get_cost_volume_capex_expr(self, m, tec):
        """Get storage-volume cost for technology.

        :param m: Model.
        :param tec: Technology.
        :type m: :py:class:`pyomo.Model`
        :type tec: str

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return self.input['volume_capex'][tec] * m.volume[tec]

    def _get_cost_charge_capacity_expr(self, m, tec):
        """Get charging-capacity cost for technology.

        :param m: Model.
        :param tec: Technology.
        :type m: :py:class:`pyomo.Model`
        :type tec: str

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        charge_capacity_cost_tec = (self.input['charge_capex'][tec] +
                                    self.input['charge_om_fixed'][tec])

        return m.charge_capacity[tec] * charge_capacity_cost_tec

    def _get_cost_om_variable_full_expr(self, m, tec):
        """Get generation cost over full time series for technology.

        :param m: Model.
        :param tec: Technology.
        :type m: :py:class:`pyomo.Model`
        :type tec: str

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return sum(self._get_cost_om_variable_expr(m, tec, t) for t in m.time)

    def _get_cost_om_variable_expr(self, m, tec, t):
        """Get generation cost at timestamp for technology.

        :param m: Model.
        :param tec: Technology.
        :param t: Timestamp.
        :type m: :py:class:`pyomo.Model`
        :type tec: str
        :type t: :py:class:`pyomo.Set`

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return m.generation_nonvre[tec, t] * self.input['om_variable'][tec]

    def _get_generation_total_expr(self, m, t):
        """Get total generation.

        :param m: Model.
        :param t: Timestamp.
        :type m: :py:class:`pyomo.Model`
        :type t: :py:class:`pyomo.Set`

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        # Get total non-VRE generation
        generation_nonvre_total = sum(m.generation_nonvre[tec, t]
                                      for tec in m.nonvre)

        # Get total VRE generation
        generation_vre_total = sum(
            m.capacity[tec] * self.input['capacityfactor'][tec][t] for tec in m.vre_total)

        # Get total VRE generation from regional
        if self.regional:
            generation_vre_total += sum(
                self._get_generation_vre_regional_expr(m, tec, reg, t)
                for tec in m.vre_regional for reg in m.region)

        # Return sum
        return generation_nonvre_total + generation_vre_total

    def _get_generation_vre_expr(self, m, tec, t):
        """Get constraint on VRE generation.

        :param m: Model.
        :param tec: Technology.
        :param t: Timestamp.
        :type m: :py:class:`pyomo.Model`
        :type tec: str
        :type t: :py:class:`pandas.Timestamp`

        :returns: Expression.
        :rtype: :py:class:`pyo.Expression`
        """
        return m.capacity[tec] * self.input['capacityfactor'][tec][t]


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
        var = pyo.Var(*indices, domain=pyo.NonNegativeReals)
        setattr(self.model, variable_name, var)

    def _add_set_coord(self, set_name, ordered=True, initialize=None,
                       **kwargs):
        """Add set and coordinates to model.

        :param set_name: Set name.
        :param ordered: Whether set is ordered.
        :param initialize: Initialization array.
        :type set_name: str
        :type ordered: bool
        :type initialize: sequence
        """
        # Add set
        s = pyo.Set(ordered=ordered, initialize=initialize)
        setattr(self.model, set_name, s)

        # Add coordinate
        self.coords[set_name] = set_name, s.value_list

    def _fix_variables(self):
        """Fix variables from configuration."""
        cfg_fix = self.cfg.get('fix')
        if cfg_fix is not None:
            for variable_name, component_values in cfg_fix.items():
                for component_name, value in component_values.items():
                    if isinstance(value, str):
                        # Get value from input data for variable and component
                        value = self.input[value][component_name]
                    getattr(self.model, variable_name)[
                        component_name].fix(value)

    def _bound_variables(self):
        """Bound variables from configuration."""
        cfg_bound = self.cfg.get('bound')
        if cfg_bound is not None:
            self._bound_variables_side(cfg_bound, 'upper')
            self._bound_variables_side(cfg_bound, 'lower')

    def _bound_variables_side(self, cfg_bound, side):
        """Bound variables up or low from configuration.

        :param cfg_bound: Bounds configuration.
        :param side: Bound side (`'upper'` or `'lower'`).
        :type cfg_bound: mapping
        :type side: str
        """
        cfg_side = cfg_bound.get(side) or {}
        for variable_name, component_values in cfg_side.items():
            for component_name, value in (component_values or {}).items():
                if isinstance(value, str):
                    # Get value from input data for variable and component
                    value = self.input[value][component_name]
                if side == 'upper':
                    getattr(self.model, variable_name)[
                        component_name].setub(value)
                elif side == 'lower':
                    getattr(self.model, variable_name)[
                        component_name].setlb(value)

    def _add_available_components(self, **kwargs):
        """Add available components."""
        cfg_comp = self.cfg.get('component') or {}
        context_component_name = cfg_comp.get('demand')
        if context_component_name is not None:
            # Add demand component manager.
            self.context_component_demand = self.med.context_components[
                context_component_name]

        cfg_comp_cf = cfg_comp.get('capacityfactor')
        if cfg_comp_cf is not None:
            # Initialize generation component managers
            self.context_component_generation = OrderedDict()

            # Initialize generation component names
            self.component_generation_names = set()

            # Add generation components
            for context_component_name in cfg_comp_cf:
                context_component = self.med.context_components[
                    context_component_name]
                self.context_component_generation[
                    context_component_name] = context_component
                self.component_generation_names.add(
                    context_component.component_name)


def _get_capacity_total(m, tec):
    """Get total capacity for technology.

    :param m: Model.
    :param tec: Technology.
    :type m: :py:class:`pyomo.Model`
    :type tec: str

    :returns: Expression.
    :rtype: :py:class:`pyo.Expression`
    """
    if tec in m.vre_regional:
        cap = sum(m.capacity_regional[tec, reg] for reg in m.region)
    else:
        cap = m.capacity[tec]

    return cap


class Input(InputBase):
    """Optimization problem input data source."""

    def __init__(self, optimizer, **kwargs):
        """Initialize input data source.

        :param optimizer: Optimizer.
        :type optimizer: : py:class:`.optimization.OptimizerBase`
        """
        super(Input, self).__init__(optimizer, **kwargs)

        #: Whether problem is regional
        self.regional = bool(self.cfg.get('regional_components'))

        # Add variables
        variable_names = [
            'capacityfactor', 'demand', 'inflows_lake',
            'reserve_requirement', 'capacity_existing', 'capacity_max',
            'capex', 'volume_capex', 'om_fixed', 'om_variable',
            'efficiency_charging', 'efficiency_discharging',
            'max_generation_year', 'uncertainty', 'variation_factor',
            'charge_capex', 'charge_om_fixed'
        ]
        if self.regional:
            variable_names.append('capacity_factor_regional')
        self.update_variables(variable_names)

        # Time
        raw_start_date = pd.Timestamp('2000-01-01')
        raw_end_date = pd.Timestamp('2018-01-01')

        #: Raw time used to read.
        time_full = pd.date_range(raw_start_date, raw_end_date,
                                  freq='H', closed='left', name='time')
        # Remove bissextile extra day
        self.time_full = time_full[~((time_full.day == 29) &
                                     (time_full.month == 2))]
        self.months_full = pd.Index(np.unique(_get_month(self.time_full)),
                                    name='months')

        # Select period
        self._start_date = pd.Timestamp(
            self.cfg.get('start_date') or raw_start_date)
        self._end_date = (pd.Timestamp(
            self.cfg.get('end_date') or raw_end_date)
            - pd.Timedelta('1 second'))
        self.time = self.time_full[(self.time_full >= self._start_date) &
                                   (self.time_full <= self._end_date)]
        month_of_hour = _get_month(self.time)
        self.months = pd.Index(np.unique(month_of_hour))
        month_timestamps = pd.Series(self.time, index=month_of_hour)
        self.month_timestamps = pd.Series(
            dict([(group[0], group[1].tolist())
                  for group in month_timestamps.groupby(level=0)]))
        self.time_without_last = self.time[:-1]
        self.next_time = dict(zip(self.time[:-1], self.time[1:]))

        # Set
        regional_components = self.cfg.get('regional_components')
        self.vre_regional = (regional_components if self.regional else
                             pd.Index([]))
        for index_name, values in self.cfg['indices'].items():
            setattr(self, index_name, pd.Index(values, name=index_name))

        # Get non-VRE index
        self.nonvre = self.technology.difference(self.vre)
        self.nonvre.name = 'nonvre'

        # Remove regional VREs
        self.tec_total = self.technology.difference(self.vre_regional)
        self.tec_total.name = 'tec_total'
        self.vre_total = self.vre.difference(self.vre_regional)
        self.vre_total.name = 'vre_total'

        # Biogas singleton
        self.biogas = pd.Index(['biogas'], name='biogas')

        context_component_names = list(
            self.optimizer.context_component_generation)
        #: Region index.
        self.region = (pd.Index(self.optimizer.context_component_generation[
            context_component_names[0]].place_names)
            if self.regional else [])

        # Map variables to their index
        self.variable_indices = OrderedDict({
            'capacityfactor': [self.vre, self.time_full],
            'demand': [self.time_full],
            'inflows_lake': [self.months_full],
            'reserve_requirement': [self.vre],
            'capacity_existing': [self.technology],
            'capacity_max': [self.technology],
            'capex': [self.technology],
            'volume_capex': [self.storing],
            'om_fixed': [self.technology],
            'om_variable': [self.technology],
            'charge_capex': [self.storing],
            'charge_om_fixed': [self.storing],
            'efficiency_charging': [self.storing],
            'efficiency_discharging': [self.storing]})
        if self.regional:
            self.variable_indices['capacity_factor_regional'] = [
                self.vre_regional, self.region, self.time_full]

        #: VRE capacities.
        self.capacity_vre = None

    def parse(self, **kwargs):
        """Parse optimization problem input data.

        :returns: Dataset.
        :rtype: dict
        """
        # Add constants
        ds = self.cfg['constants'].copy()

        # Read arrays
        for variable_name in self.variable_component_names:
            if variable_name not in self.cfg['constants']:
                filepath_list = ensure_collection(
                    self.cfg['input_path'].get(variable_name), list)
                if filepath_list is not None:
                    filepath = Path(*filepath_list)
                    read_csv_kwargs = self.cfg.get(
                        'load_read_csv_kwargs') or {}
                    indices = self.variable_indices[variable_name]
                    ds[variable_name] = _read_csv_variable(
                        filepath, indices, self.time_full,
                        transform=_manage_months, months_full=self.months_full,
                        months=self.months, **read_csv_kwargs)

        # Add available-component input data
        ds_new = self.parse_available_components(
            regional=self.regional, dem_conv=1.e-3, **kwargs)
        for variable_name, s in ds_new.items():
            if self.regional and (variable_name == 'capacityfactor'):
                ds[variable_name + '_regional'] = s
            else:
                ds[variable_name].loc[s.index] = s

        return ds

    def get_costs_per_technology(self, ds, **kwargs):
        """Get splitted costs per technology.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: Dataset of costs per technology (M€).
        :rtype: dict
        """
        ds_cost = {}
        ds_cost['capex'] = self.get_cost_capex(ds, **kwargs)
        ds_cost['om_fixed'] = self.get_cost_om_fixed(ds, **kwargs)
        ds_cost['volume'] = self.get_cost_volume_capex(ds, **kwargs)
        ds_cost['charge_capacity'] = self.get_cost_charge_capacity(
            ds, **kwargs)
        ds_cost['om_variable'] = self.get_cost_om_variable_full(ds, **kwargs)

        return ds_cost

    @support.get_data
    def get_cost_capex(self, ds, **kwargs):
        """Get CAPEX cost of additional capacity for technology.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: CAPEX cost (M€).
        :rtype: :py:class:`xarray.DataArray`
        """
        da_cap = self.get_capacity(ds, **kwargs)
        capacity_new = da_cap - xr.DataArray(self['capacity_existing'])

        return xr.DataArray(self['capex']) * capacity_new

    @support.get_data
    def get_cost_om_fixed(self, ds, **kwargs):
        """Get O&M cost per technology.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: O&M cost (M€).
        :rtype: :py:class:`xarray.DataArray`
        """
        da_cap = self.get_capacity(ds, **kwargs)

        return xr.DataArray(self['om_fixed']) * da_cap

    @support.get_data
    def get_cost_volume_capex(self, ds, **kwargs):
        """Get storage-volume cost per technology.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: Storage-volume cost (M€).
        :rtype: :py:class:`xarray.DataArray`
        """
        return xr.DataArray(self['volume_capex']) * ds['volume']

    @support.get_data
    def get_cost_charge_capacity(self, ds, **kwargs):
        """Get charging-capacity cost per technology.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: Charging-capacity cost (M€).
        :rtype: :py:class:`xarray.DataArray`
        """
        charge_capacity_cost = self['charge_capex'] + self['charge_om_fixed']

        return ds['charge_capacity'] * xr.DataArray(charge_capacity_cost)

    @support.get_data
    def get_cost_om_variable_full(self, ds, **kwargs):
        """Get generation cost over full time series per technology.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: Generation cost over full time series (M€).
        :rtype: :py:class:`xarray.DataArray`
        """
        return self.get_cost_om_variable(ds, **kwargs).sum('time')

    @support.get_data
    def get_cost_om_variable(self, ds, **kwargs):
        """Get generation cost over full per technology and timestamp.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: Generation cost (M€).
        :rtype: :py:class:`xarray.DataArray`
        """
        # Get variable O&M costs
        da_om_variable = xr.DataArray(self['om_variable']).sel(
            technology=self.nonvre).rename(technology='nonvre')

        # Return product
        return ds['generation_nonvre'] * da_om_variable

    def get_generation(self, ds, **kwargs):
        """Get generation (VRE and non-VRE).

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: Generation (GWh).
        :rtype: :py:class:`xarray.DataArray`
        """
        # Get VRE generation
        da_gen_vre = self.get_generation_vre(ds, **kwargs)

        # Get non-VRE generation
        dim = 'technology'
        da_gen_nonvre = ds['generation_nonvre'].rename(nonvre=dim)

        # Return concatenation
        return xr.concat([da_gen_nonvre, da_gen_vre], dim=dim)

    def get_generation_vre(self, ds, **kwargs):
        """Get VRE generation.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: VRE generation (GWh).
        :rtype: :py:class:`xarray.DataArray`
        """
        # Get non-regional generation
        da_cf_total = self.get_capacity_factor(regional=False, **kwargs)
        if self.regional:
            da_cap_total = ds['capacity'].rename(tec_total='technology')
            da_cap_total = da_cap_total.sel(technology=self.vre_total)
        else:
            da_cap_total = ds['capacity'].sel(technology=self.vre)
        da_gen = da_cap_total * da_cf_total

        if self.regional:
            # Get regional generation
            da_cf_regional = self.get_capacity_factor(regional=True, **kwargs)
            da_cap_regional = ds['capacity_regional'].rename(
                vre_regional='technology')
            da_gen_regional = (da_cap_regional * da_cf_regional).sum('region')

            # Concatenate
            da_gen = xr.concat([da_gen, da_gen_regional], dim='technology')

        return da_gen

    @support.get_data
    def get_capacity_factor(self, regional=False, **kwargs):
        """Get (regional) capacity factors.

        :param regional: Whether VREs are regional.
        :type regional: bool

        :returns: Capacity factors.
        :rtype: :py:class:`xarray.DataArray`
        """
        variable_name = 'capacityfactor'
        if regional:
            variable_name += '_regional'

        return self[variable_name].to_xarray().rename(vre='technology')

    @support.get_data
    def get_capacity(self, ds, **kwargs):
        """Get capacity.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: Capacity factors.
        :rtype: :py:class:`xarray.DataArray`
        """
        da_cap = ds['capacity']
        if self.regional:
            # Get non-regional capacity
            da_cap = da_cap.rename(tec_total='technology')
            da_cap = da_cap.sel(technology=self.tec_total)

            # Get regional capacity
            da_cap_regional = ds['capacity_regional'].rename(
                vre_regional='technology').sum('region')
            da_cap = xr.concat([da_cap, da_cap_regional], dim='technology')

        return da_cap

    @support.get_data
    def get_demand(self, **kwargs):
        """Get demand.

        :returns: Demand (GWh).
        :rtype: :py:class:`xarray.DataArray`
        """
        return xr.DataArray(self['demand'])

    def get_load(self, ds, **kwargs):
        """Get load as demand plus storage.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: Load (GWh).
        :rtype: :py:class:`xarray.DataArray`
        """
        return self.get_demand(**kwargs) + ds['storage'].sum('storing')

    def get_residual_load(self, ds, **kwargs):
        """Get load as demand plus storage.

        :param ds: Solution dataset.
        :type ds: mapping

        :returns: Load (GWh).
        :rtype: :py:class:`xarray.DataArray`
        """
        # Get load
        da_load = self.get_load(ds, **kwargs)

        # Get VRE generation
        da_gen_vre = self.get_generation_vre(ds, **kwargs)

        # Return the difference
        return da_load - da_gen_vre.sum('technology')

    def set_capacity_vre(self, capacity_vre, **kwargs):
        """Set VRE capacities.

        :param capacity_vre: VRE capacities.
        :type capacity_vre: :py:class:`xarray.DataArray`
        """
        # Update VRE capacity in coarse-grain optimizer
        self.capacity_vre = capacity_vre

        # Enable load task
        self.task_mng.set_all(True)
        self.optimizer.solution.task_mng.set_all(True)

    def get_data_postfix(self, **kwargs):
        """Get data postfix.

        :returns: Postfix.
        :rtype: str
        """
        # Get user-defined postfix
        postfix = self.cfg.get('postfix') or ''

        return postfix

def _get_month(time, **kwargs):
    """Convert time stamp to month index.

    :param time: Time index.
    :type time: :py:class:`pandas.DatetimeIndex`

    :returns: Month index.
    :rtype: :py:class:`pandas.Index`
    """
    return time.month + time.year * 12


def _get_index_names(indices):
    index_names = [idx.name for idx in indices]

    return index_names


def _read_csv_variable(filepath, indices, time_index, transform=None, **kwargs):

    # Read CSV
    s = pd.read_csv(
        filepath, index_col=list(range(len(indices))),
        **kwargs).squeeze()

    _manage_time_index(time_index)

    if transform is not None:
        transform(s, indices, **kwargs)

    # Set named (multi-)index
    index_names = _get_index_names(indices)
    if len(indices) > 1:
        s.index.names = index_names
    else:
        s.index.name = index_names[0]

    # Sort index
    s = s.sort_index()

    return s


def _manage_months(s, indices, months_full=None, months=None):
    index_names = _get_index_names(indices)
    if 'months' in index_names:
        # Make month index
        s.index = months_full[s.index - 1]

        # Select months
        s = s.loc[months]


def _manage_time_index(s, indices, time_index, start_date=None, end_date=None):
    index_names = [idx.name for idx in indices]
    if 'time' in index_names:
        if len(indices) > 1:
            # Make time index
            time_index = time_index[s.index.levels[1]]
            s.index = pd.MultiIndex.from_product([
                s.index.levels[0], time_index])

            # Select period from MultiIndex
            s = s.loc[(slice(None), slice(start_date, end_date))]
        else:
            # Make time index
            s.index = time_index[s.index]
            s.name = 'time'

            # Select period from single index
            s = s.loc[slice(start_date, end_date)]

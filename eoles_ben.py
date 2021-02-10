from pulp import LpMaximize, LpMinimize, LpProblem, LpStatus, lpSum, LpVariable

import pyomo.environ as pyo
import pandas as pd

vre_profiles_2006_data = pd.read_csv("inputs/vre_profiles2006.csv", index_col=[0,1], squeeze=True)
vre_profiles_2006_data.head()

vre_profiles_2006_data = pd.read_csv("inputs/vre_profiles2006.csv", names=["Type", "Hour", "Production Profile"])
vre_profiles_2006_data = vre_profiles_2006_data.pivot(index='Hour', columns='Type')
# Preview the first 5 lines of the loaded data
vre_profiles_2006_data.describe()

prob = LpProblem("EOLES Optimization", LpMinimize)
Technologies = ["pv", "onshore", "offshore", "river", "lake", "biogas", "phs", "battery", "methanation"]

self.model = pyo.ConcreteModel()
self.model.cost = pyo.Objective(rule=self.cost_rule)





def adequacy_rule(m, t):
    """Get constraint for adequacy.

    :param m: Model.
    :param t: Timestamp.
    :type m: :py:class:`pyomo.Model`
    :type t: :py:class:`pyomo.Set`

    :returns: Expression.
    :rtype: :py:class:`pyo.Expression`
    """
    # Get total generation
    generation_total = self._get_generation_total_expr(m, t)

    # Get total storage
    storage_total = sum(m.storage[tec, t] for tec in m.storing)

    return generation_total >= self.input['demand'][t] + storage_total
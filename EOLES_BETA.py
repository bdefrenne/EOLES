from pyomo.environ import *
import pandas as pd
# from pyomo.repn.plugins.baron_writer import NonNegativeReals

# TODO: Faire fonctionner toutes les contraintes
# TODO: Donner un nom clair aux variables
# TODO: Optimiser sur 18 ans et pas 1 an


def main():
    # LOAD THE DATA
    production_profile = pd.read_csv(
        "inputs/vre_profiles2006.csv", index_col=[0, 1], squeeze=True, header=None)
    existing_capacities = pd.read_csv(
        "inputs/existing_capas.csv", index_col=0, squeeze=True, header=None)
    maximum_capacities = pd.read_csv(
        "inputs/max_capas.csv", index_col=0, squeeze=True, header=None)
    capex = pd.read_csv("inputs/annuities.csv", index_col=0,
                        squeeze=True, header=None)
    capex_storage = pd.read_csv(
        "inputs/str_annuities.csv", index_col=0, squeeze=True, header=None)
    demand = pd.read_csv("inputs/demand2006.csv", index_col=0,
                         squeeze=True, header=None)
    lake_inflows = pd.read_csv("inputs/lake_inflows.csv",
                               index_col=0, squeeze=True, header=None)
    gene_river = pd.read_csv("inputs/run_of_river.csv",
                             index_col=0, squeeze=True, header=None)
    reserve_requirements = pd.read_csv(
        "inputs/reserve_requirements.csv", index_col=0, squeeze=True, header=None)
    fOM = pd.read_csv("inputs/fO&M.csv", index_col=0,
                      squeeze=True, header=None)
    vOM = pd.read_csv("inputs/vO&M.csv", index_col=0,
                      squeeze=True, header=None)

    # PARAMETERS
    months_hours = {'jan': range(0, 744), 'feb': range(744, 1440), 'mar': range(1440, 2184), 'apr': range(2184, 2904),
                    'may': range(2904, 3648), 'jun': range(3648, 4368), 'jul': range(4368, 5112),
                    'aug': range(5112, 5856), 'sep': range(5856, 6576), 'oct': range(6576, 7320), 'nov': range(7320, 8040),
                    'dec': range(8039, 8760)}
    s_capex = {'phs': 0, 'battery': 0, 'methanation': 84.16086}
    s_opex = {'phs': 0, 'battery': 0, 'methanation': 59.25}
    eta_in = {'phs': 0.95, 'battery': 0.9, 'methanation': 0.59}
    eta_out = {'phs': 0.9, 'battery': 0.95, 'methanation': 0.45}
    pump_capa = 9.3
    max_phs = 0.18
    max_biogas = 15 * 1000
    load_uncertainty = 0.01
    delta = 0.1

    # INIT MODEL
    model = ConcreteModel()

    # DEFINE THE SETS
    first_hour = 0
    last_hour = 8759
    model.time = RangeSet(first_hour, last_hour)
    model.time_without_last = RangeSet(first_hour, last_hour - 1)
    model.technologies = Set(
        initialize=["pv", "onshore", "offshore", "river", "lake", "biogas", "gas", "phs", "battery", "methanation"])
    model.months = Set(initialize=["jan", "feb", "mar", "apr",
                                   "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"])
    model.generations = Set(
        initialize=["pv", "onshore", "offshore", "river", "lake", "biogas", "gas"])
    model.vre = Set(initialize=["pv", "onshore", "offshore"])
    model.ncomb = Set(initialize=["pv", "onshore",
                                  "offshore", "river", "lake", "phs", "battery"])
    model.comb = Set(initialize=["biogas", "methanation"])
    model.storages = Set(initialize=["phs", "battery", "methanation"])
    model.frr = Set(initialize=["lake", "phs", "battery", "gas"])

    # Define variables
    model.capacities = Var(
        model.technologies, within=NonNegativeReals, initialize=0)
    model.generation = Var(((gen, hour) for gen in model.technologies for hour in model.time), within=NonNegativeReals,
                           initialize=0)
    model.storage_power_capacity = Var(
        model.storages, within=NonNegativeReals, initialize=0)
    model.storage_energy_capacity = Var(
        model.storages, within=NonNegativeReals, initialize=0)  # STORED
    model.storage = Var(((storage, hour) for storage in model.storages for hour in model.time), within=NonNegativeReals,
                        initialize=0)  # 'hourly electricity input of battery storage GW'
    model.stored = Var(((storage, hour) for storage in model.storages for hour in model.time), within=NonNegativeReals,
                       initialize=0)  # 'energy stored in each storage technology in GWh' = Stage of charge
    model.reserve = Var(((reserve, hour) for reserve in model.frr for hour in model.time), within=NonNegativeReals,
                        initialize=0)

    # gene_vre(vre,h)..                GENE(vre,h)             =e=     CAPA(vre)*load_factor(vre,h);

    def generation_vre_constraint_rule(model, t, vre):
        # RENAME apacity_factor
        return model.generation[vre, t] == model.capacities[vre] * production_profile[vre][t]

    # gene_capa(tec,h)..               CAPA(tec)               =g=     GENE(tec,h);

    def generation_capacity_constraint_rule(model, t, vre_tecs):
        return model.capacities[vre_tecs] >= production_profile[vre_tecs][t]

    # combustion(h)..                  GENE('gas',h)           =e=     sum(comb,GENE(comb,h));

    def combustion_constraint_rule(model, t):
        return model.generation['gas', t] == sum(model.generation[tec, t] for tec in model.comb)

    # capa_frr(frr,h)..                CAPA(frr)               =g=     GENE(frr,h) + RSV(frr,h);

    def frr_capacity_constraint_rule(model, t, frr_tecs):
        return model.capacities[frr_tecs] >= model.generation[frr_tecs, t] + model.reserve[frr_tecs, t]

    # storing(h,h+1,str)..             STORED(str,h+1)         =e=     STORED(str,h) + STORAGE(str,h)*eta_in(str) - GENE(str,h)/eta_out(str);

    def storing_constraint_rule(model, t, tec):
        return (model.stored[tec, t + 1] == model.stored[tec, t]
                + model.storage[tec, t] * eta_in[tec]
                - model.generation[tec, t] / eta_out[tec])

    # storage_const(str,first,last)..  STORED(str,first)       =e=     STORED(str,last)+STORAGE(str,last)*eta_in(str)-GENE(str,last)/eta_out(str);

    def storage_constraint_rule(model, tec):
        return (model.stored[tec, first_hour] == model.stored[tec, last_hour]
                + model.storage[tec, last_hour] * eta_in[tec]
                - model.generation[tec, last_hour] / eta_out[tec])

    # lake_res(m)..                    lake_inflows(m)         =g=     sum(h$(month(h) = ord(m)),GENE('lake',h))/1000;

    def lake_reserve_constraint_rule(model, months):
        gen_lake_month = sum(model.generation['lake', t]
                             for t in months_hours[months])
        return gen_lake_month <= lake_inflows[months] * 1000

    # stored_cap(str,h)..              STORED(str,h)           =l=     CAPACITY(str);

    def stored_capacity_constraint(model, t, storage_tecs):
        return model.stored[storage_tecs, t] <= model.capacities[storage_tecs]

    # storage_capa1(str,h)..           S(str)                  =g=     STORAGE(str,h);

    def storage_capacity_1_constraint_rule(model, t, storage_tecs):
        return model.storage_power_capacity[storage_tecs] >= model.storage[storage_tecs, t]

    # storage_capa2(str)..             S(str)                  =l=     CAPA(str);

    def storage_capacity_2_constraint_rule(model, storage_tecs):
        return model.storage_power_capacity[storage_tecs] <= model.capacities[storage_tecs]

    # biogas_const..                   sum(h,GENE('biogas',h)) =l=     max_biogas*1000;

    def biogas_constraint_rule(model):
        fact = 24 * 365 / len(model.time)
        gen_year_biogas = sum(model.generation['biogas', hour]
                              for hour in model.time) * fact
        return gen_year_biogas <= max_biogas

    # reserves(h)..                    sum(frr, RSV(frr,h))    =e=     sum(vre,epsilon(vre)*CAPA(vre))+ demand(h)*load_uncertainty*(1+delta);

    def reserves_constraint_rule(model, t):
        return sum(model.reserve[frr, t] for frr in model.frr) == sum(
            reserve_requirements[vre] * model.capacities[vre] for vre in model.vre) + demand[t] * load_uncertainty * (1 + delta)

    # adequacy(h)..                    sum(ncomb,GENE(ncomb,h))+GENE('gas',h)    =g=     demand(h)+sum(str,STORAGE(str,h)) ;

    def adequacy_constraint_rule(model, t):
        return sum(model.generation[ncomb, t] for ncomb in model.ncomb) + model.generation["gas", t] >= demand[t]\
            + sum(model.storage[storage, t] for storage in model.storages)

    # obj..   COST               =e=   (sum(tec,(CAPA(tec)-capa_ex(tec))*capex(tec))+ sum(str,CAPACITY(str)*capex_en(str))+sum(tec,(CAPA(tec)*fOM(tec)))+ sum(str,S(str)*(s_capex(str)+s_opex(str))) + sum((tec,h),GENE(tec,h)*vOM(tec)))/1000;

    def objective_rule(model):
        return sum((model.capacities[tec] - existing_capacities[tec]) * capex[tec] for tec in model.technologies) \
            + sum(model.capacities[storage] * capex_storage[storage] for storage in model.storages)\
            + sum(model.capacities[tec] * fOM[tec] for tec in model.technologies)\
            + sum(model.storage_power_capacity[storage] * (s_opex[storage] + s_capex[storage]) for storage in model.storages)\
            + sum(sum(model.generation[tec, hour] * vOM[tec]
                      for hour in model.time) for tec in model.technologies)

    print('Defining constraints')
    model.generation_vre_constraint = Constraint(
        model.time, model.vre, rule=generation_vre_constraint_rule)
    model.generation_capacity_constraint = Constraint(
        model.time, model.vre, rule=generation_capacity_constraint_rule)
    model.combustion = Constraint(model.time, rule=combustion_constraint_rule)
    model.frr_capacity_constraint = Constraint(
        model.time, model.frr, rule=frr_capacity_constraint_rule)

    print('Defining storing')
    model.storing_constraint = Constraint(
        model.time_without_last, model.storages, rule=storing_constraint_rule)

    print('Defining storage')
    model.storage_constraint = Constraint(
        model.storages, rule=storage_constraint_rule)

    print('Defining lake')
    model.lake_reserve_constraint = Constraint(
        model.months, rule=lake_reserve_constraint_rule)

    print('Defining stored capacity')
    model.stored_capacity_constraint = Constraint(
        model.time, model.storages, rule=stored_capacity_constraint)
    print('Defining storage capacity 1')
    model.storage_capacity_1_constraint = Constraint(
        model.time, model.storages, rule=storage_capacity_1_constraint_rule)
    print('Defining storage capacity 2')
    model.storage_capacity_2_constraint = Constraint(
        model.storages, rule=storage_capacity_2_constraint_rule)
    print('Defining biogas')
    model.biogas_constraint = Constraint(
        rule=biogas_constraint_rule)
    print('Defining reserves')
    model.reserves_constraint = Constraint(
        model.time, rule=reserves_constraint_rule)
    print('Defining adequacy')
    model.adequacy_constraint = Constraint(
        model.time, rule=adequacy_constraint_rule)

    print('Defining objective')
    model.objective = Objective(rule=objective_rule)

    # Define solver
    opt = SolverFactory('ipopt')

    print('Solving')
    opt.solve(model, tee=True)


if __name__ == "__main__":
    main()

from petabtests import *
from petab.C import *
import petab

import pandas as pd


test_id = 5

# problem --------------------------------------------------------------------

model = DEFAULT_MODEL_FILE

condition_df = pd.DataFrame(data={
    CONDITION_ID: ['c0', 'c1'],
    'offset_A': ['offset_A_c0', 'offset_A_c1'],
}).set_index([CONDITION_ID])

measurement_df = pd.DataFrame(data={
    OBSERVABLE_ID: ['obs_a', 'obs_a'],
    SIMULATION_CONDITION_ID: ['c0', 'c1'],
    TIME: [10, 10],
    MEASUREMENT: [2.1, 3.2]
})

observable_df = pd.DataFrame(data={
    OBSERVABLE_ID: ['obs_a'],
    OBSERVABLE_FORMULA: ['A + offset_A'],
    NOISE_FORMULA: [1]
}).set_index([OBSERVABLE_ID])

parameter_df = pd.DataFrame(data={
    PARAMETER_ID: ['a0', 'b0', 'k1', 'k2', 'offset_A_c0', 'offset_A_c1'],
    PARAMETER_SCALE: [LIN] * 6,
    LOWER_BOUND: [0] * 6,
    UPPER_BOUND: [10] * 6,
    NOMINAL_VALUE: [1, 0, 0.8, 0.6, 2, 3],
    ESTIMATE: [1] * 6,
}).set_index(PARAMETER_ID)


# write files

write_problem(test_id=test_id,
              parameter_df=parameter_df,
              condition_dfs=[condition_df],
              observable_dfs=[observable_df],
              measurement_dfs=[measurement_df],
              sbml_files=['conversion_modified.xml'])

# solutions ------------------------------------------------------------------

simulation_df = measurement_df.copy(deep=True).rename(
    columns={MEASUREMENT: SIMULATION})
simulation_df[SIMULATION] = [analytical_a(10, 1, 0, 0.8, 0.6) + offset
                             for offset in [2, 3]]

chi2 = petab.calculate_chi2(
    measurement_df, simulation_df, observable_df, parameter_df)

llh = petab.calculate_llh(
    measurement_df, simulation_df, observable_df, parameter_df)
print(llh)
# write files

write_solution(test_id=test_id,
               chi2=chi2,
               llh=llh,
               simulation_dfs=[simulation_df])

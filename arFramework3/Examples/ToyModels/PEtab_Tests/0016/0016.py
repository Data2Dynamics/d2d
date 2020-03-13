from petabtests import *
from petab.C import *
import petab

import pandas as pd


test_id = 16

# problem --------------------------------------------------------------------

model = DEFAULT_MODEL_FILE

condition_df = pd.DataFrame(data={
    CONDITION_ID: ['c0'],
}).set_index([CONDITION_ID])

measurement_df = pd.DataFrame(data={
    OBSERVABLE_ID: ['obs_a', 'obs_b'],
    SIMULATION_CONDITION_ID: ['c0', 'c0'],
    TIME: [10, 10],
    MEASUREMENT: [0.2, 0.8]
})

observable_df = pd.DataFrame(data={
    OBSERVABLE_ID: ['obs_a', 'obs_b'],
    OBSERVABLE_FORMULA: ['A', 'B'],
    OBSERVABLE_TRANSFORMATION: [LIN, LOG],
    NOISE_FORMULA: [0.5, 0.7]
}).set_index([OBSERVABLE_ID])

parameter_df = pd.DataFrame(data={
    PARAMETER_ID: ['a0', 'b0', 'k1', 'k2'],
    PARAMETER_SCALE: [LIN] * 4,
    LOWER_BOUND: [0] * 4,
    UPPER_BOUND: [10] * 4,
    NOMINAL_VALUE: [1, 0, 0.8, 0.6],
    ESTIMATE: [1] * 4,
}).set_index(PARAMETER_ID)


# write files

write_problem(test_id=test_id,
              parameter_df=parameter_df,
              condition_dfs=[condition_df],
              observable_dfs=[observable_df],
              measurement_dfs=[measurement_df])

# solutions ------------------------------------------------------------------

simulation_df = measurement_df.copy(deep=True).rename(
    columns={MEASUREMENT: SIMULATION})
simulation_df[SIMULATION] = [
    analytical_a(10, 1, 0, 0.8, 0.6),
    analytical_b(10, 1, 0, 0.8, 0.6),
]

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

import numpy as np

from multiple_pulses_x import WrightFisherPopulationX

np.random.seed(0)

N = 30
recombination_distance = 10
ms = np.array([0.5, 0.1, 1.0])
Ts = np.array([5, 10, 11])
sources = ['New York', 'California', 'Montana']

wf = WrightFisherPopulationX(population_size=N, migration_times=Ts,
                             migration_sources=sources,
                             migration_probabilities=ms)
sources, end_points = wf.simulate_tracts(recombination_distance)

tract_lengths = np.diff(end_points)

assert len(sources) + 1 == len(end_points), (len(sources), len(end_points))

assert sources == ['New York', 'Montana', 'New York', 'Montana', 'New York',
                   'Montana', 'New York', 'Montana', 'California', 'Montana',
                   'California', 'Montana', 'New York', 'Montana',
                   'California', 'Montana', 'New York'], sources

assert list(tract_lengths) == [0.6937052355503548, 0.5072101120308574,
                               1.198369873588488, 0.1264673205330764,
                               3.549498613756529, 0.24512368850725075,
                               0.16582053991618917, 0.16873083609897055,
                               0.08690603778626116, 0.5879640412822909,
                               0.001071203147483324, 0.041461254783071944,
                               0.9555001254702935, 0.3133848760821536,
                               0.24114168151579563, 0.5564163572109422,
                               0.5612282027399917], list(tract_lengths)

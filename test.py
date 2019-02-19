import numpy as np

from multiple_pulses import WrightFisherPopulation

np.random.seed(0)

N = 30
recombination_distance = 10
ms = np.array([0.5, 0.1, 1.0])
Ts = np.array([5, 10, 11])
sources = ['New York', 'California', 'Montana']

wf = WrightFisherPopulation(population_size=N, migration_times=Ts,
                            migration_sources=sources,
                            migration_probabilities=ms)
sources, end_points = wf.simulate_tracts(recombination_distance)

tract_lengths = np.diff(end_points)

assert len(sources) + 1 == len(end_points), (len(sources), len(end_points))

assert sources == ['New York', 'Montana', 'New York', 'Montana', 'New York',
                   'Montana', 'New York', 'Montana', 'California', 'Montana',
                   'New York', 'Montana', 'New York', 'Montana', 'California',
                   'New York', 'Montana', 'New York', 'Montana', 'New York',
                   'Montana', 'New York', 'Montana'], sources

assert list(tract_lengths) == [0.01167693415614067, 0.004085117336844989,
                               0.23152448655873453, 0.027605737526490037,
                               0.51768293690225110, 0.076540191559524300,
                               0.37154370977900375, 0.307656327259098000,
                               0.28413347663820265, 0.686775190110418400,
                               1.27106413677881000, 0.327957906775034200,
                               0.14595079786078635, 0.043219244172644444,
                               0.02309194406985071, 1.046841705061321000,
                               0.33333299721227494, 0.642401802807860600,
                               1.49783666832221480, 0.174885326138086940,
                               0.32459367147154694, 0.072292927473903030,
                               1.57730676402895750], list(tract_lengths)

import sys
import argparse as ap
import numpy as np


class Population(object):
    def __init__(self, population_size, migration_times,
                 migration_probabilities, migration_sources):
        self.population_size = population_size

        assert migration_probabilities[-1] == 1.0, (
            "Final migration "
            "probability must be 1.0. It is %f" % migration_probabilities[-1])
        self.migration_probabilities = dict(zip(migration_times,
                                                migration_probabilities))
        self.migration_sources = dict(zip(migration_times, migration_sources))
        self.population = {}

        sys.setrecursionlimit(3 * migration_times[-1] + 100)

    def simulate_mosaic(self, recombination_distance):
        local_ancestors = [self.get_person(generation=1)]
        while local_ancestors[-1].source_population is None:
            local_ancestors.append(local_ancestors[-1].get_parent(0, self))

        tile_sources = []
        tile_boundaries = [0]

        position = local_ancestors[-1].get_tile_length()

        while position < recombination_distance:
            recombiner = np.random.choice(local_ancestors[:-1])
            while local_ancestors[-1] is not recombiner:
                local_ancestors.pop().last_seen_at = position

            recombiner.recombine(position)

            while local_ancestors[-1].source_population is None:
                ancestor = local_ancestors[-1].get_parent(position, self)
                local_ancestors.append(ancestor)

            tile_sources.append(ancestor.source_population)
            tile_boundaries.append(position)
            position += local_ancestors[-1].get_tile_length()

        tile_sources.append(local_ancestors[-1].source_population)
        tile_boundaries.append(recombination_distance)
        return tile_sources, tile_boundaries

    def simulate_tracts(self, r):
        tile_sources, tile_boundaries = self.simulate_mosaic(r)
        assert len(tile_sources) + 1 == len(tile_boundaries)

        current_source = tile_sources[0]

        tract_sources = [current_source]
        tract_boundaries = [0]
        for source, boundary in zip(tile_sources, tile_boundaries[:-1]):
            if source == current_source:
                continue
            tract_sources.append(source)
            tract_boundaries.append(boundary)
            current_source = source

        tract_boundaries.append(tile_boundaries[-1])

        return tract_sources, tract_boundaries

    def get_person(self, generation):
        person_id = (generation, np.random.randint(2 * self.population_size))

        if person_id in self.population:
            return self.population[person_id]

        source_population = None
        if generation in self.migration_probabilities:
            if np.random.uniform() < self.migration_probabilities[generation]:
                source_population = self.migration_sources[generation]

        self.population[person_id] = Person(generation, source_population)
        return self.population[person_id]


class Person(object):
    def __init__(self, generation, source_population):
        self.generation = generation
        self.source_population = source_population

        self.last_seen_at = -np.inf
        self.copying = None
        self.not_copying = None

    def swap_copying(self):
        self.copying, self.not_copying = self.not_copying, self.copying

    def recombine(self, position):
        self.last_seen_at = position
        self.swap_copying()

    def get_parent(self, position, population):
        distance_from_last_seen = position - self.last_seen_at
        self.last_seen_at = position

        swap_probability = .5 - .5 * np.exp(-2 * distance_from_last_seen)
        if np.random.uniform() < swap_probability:
            self.swap_copying()
        if self.copying is None:
            self.copying = population.get_person(self.generation + 1)
        return self.copying

    def get_tile_length(self):
        assert self.source_population is not None
        return np.random.exponential(1.0 / self.generation)


if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Simulate migrant tracts.')
    parser.add_argument(
        '-N', help='effective number of diploid people in population')
    parser.add_argument('-r', help='recombination distance, in Morgans')
    parser.add_argument('-m', nargs='+', help='migration probabilities')
    parser.add_argument('-T', nargs='+', help='migration times')
    parser.add_argument('-s', nargs='+', help='source population labels')
    parser.add_argument('-q', help='quantidade de cromossomos')
    parser.add_argument('-c', help='numero do chr')
    args = parser.parse_args()
    quantidade = int(args.q)
    cromossomo = str(args.c)

    for j in range(quantidade):
        N = int(args.N)
        r = float(args.r)
        sources = args.s
        Ts = np.array(args.T, dtype='int')
        ms = np.array(args.m, dtype='float')

        wf = WrightFisherPopulation(N=N)
        sim_count = 1
        sources, end_points = wf.simulate_tracts(
            r=r, ms=ms, Ts=Ts, sources=sources)
        tract_lengths = np.diff([0] + end_points)

        for i in range(len(sources)):
            print ('%s\t%f\t%s\t%s' %
                   (sources[i], tract_lengths[i], j, cromossomo))

import sys
import argparse as ap
import numpy as np


class WrightFisherPopulation:
    def __init__(self, N):
        self.N = N
        self.population = {}

    def simulate_mosaic(self, r, ms, Ts, sources):
        sys.setrecursionlimit(3 * Ts[-1] + 100)
        self.ms = dict(zip(Ts, ms))
        self.sources = dict(zip(Ts, sources))
        self.population = {}

        proband = Person(wf=self, child=None)
        contributor, next_position = proband.get_ancestor(position=0)
        tile_sources = [contributor.source]
        tile_boundaries = []
        while next_position < r:
            tile_boundaries.append(next_position)
            contributor, next_recombination = (
                contributor.child.recombine(next_position))
            next_position += next_recombination
            tile_sources.append(contributor.source)
        tile_boundaries.append(r)
        return tile_sources, tile_boundaries

    def simulate_tracts(self, r, ms, Ts, sources):
        tile_sources, tile_boundaries = self.simulate_mosaic(r, ms, Ts,
                                                             sources)
        switch_points = [s != s_next
                         for s, s_next in
                         zip(tile_sources[:-1], tile_sources[1:])]
        switch_points.append(True)

        tract_sources = [tile_sources[i] for i in range(len(tile_sources))
                         if switch_points[i]]
        tract_boundaries = [tile_boundaries[i] for i in
                            range(len(tile_sources)) if switch_points[i]]
        return tract_sources, tract_boundaries


class Person:
    def __init__(self, wf=None, child=None):
        self.wf = wf or child.wf
        self.generation = child.generation + 1 if child else 0
        self.child = child
        self.copying = None
        self.not_copying = None
        self.last_seen_at = -np.inf
        try:
            if np.random.uniform() < self.wf.ms[self.generation]:
                self.from_source_population = True
                self.source = self.wf.sources[self.generation]
            elif self.generation == max(self.wf.ms):
                self.from_source_population = True
                self.source = ''
            else:
                self.from_source_population = False
        except KeyError:
            self.from_source_population = False

    def swap_copying(self):
        self.copying, self.not_copying = self.not_copying, self.copying

    def get_copying(self):
        if self.copying is None:
            parent_ID = np.random.randint(2 * self.wf.N)
            parent_key = (self.generation + 1, parent_ID)

            if parent_key not in self.wf.population:
                self.wf.population[parent_key] = Person(child=self)
            self.copying = self.wf.population[(self.generation + 1, parent_ID)]
        return self.copying

    def recombine(self, position):
        if np.random.uniform() < 1.0 / (self.generation + 1):
            self.swap_copying()
            return self.get_copying().get_ancestor(position)
        else:
            self.last_seen_at = position
            return self.child.recombine(position)

    def get_ancestor(self, position):
        if self.from_source_population:
            return self, np.random.exponential(1.0 / self.generation)
        else:
            distance_from_last_seen = position - self.last_seen_at
            swap_probability = .5 - .5 * np.exp(-2 * distance_from_last_seen)
            if np.random.uniform() < swap_probability:
                self.swap_copying()
                self.get_copying().child = self
            return self.get_copying().get_ancestor(position)


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

        # print i;
        wf = WrightFisherPopulation(N=N)
        sim_count = 1
        sources, end_points = wf.simulate_tracts(
            r=r, ms=ms, Ts=Ts, sources=sources)
        tract_lengths = np.diff([0] + end_points)

        for i in range(len(sources)):
            print ('%s\t%f\t%s\t%s' %
                   (sources[i], tract_lengths[i], j, cromossomo))

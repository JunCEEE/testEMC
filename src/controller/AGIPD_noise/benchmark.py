"""diffr.h5 reading benchmark"""

import pytest
import h5py
from SimEx.Analysis.DiffractionAnalysis import DiffractionAnalysis

diffr_data = '/gpfs/exfel/data/user/juncheng/EMCProject/src/controller/s10/diffr.h5'


def SimExRead(idxs):
    diffr_analyzer = DiffractionAnalysis(diffr_data, idxs, poissonize=False)
    diffr_pattern = diffr_analyzer.numpyPattern()


def h5pyRead(idxs):
    with h5py.File(diffr_data, 'r') as f:
        for idx in idxs:
            grp_name = 'data/{:07}/diffr'.format(idx)
            pattern = f[grp_name][...]


def test_h5pyRead_10(benchmark):
    result = benchmark(h5pyRead, range(1, 10))


def test_h5pyRead_100(benchmark):
    result = benchmark(h5pyRead, range(1, 50))


def test_h5pyRead_200(benchmark):
    result = benchmark(h5pyRead, range(1, 100))


def test_SimExRead_10(benchmark):
    result = benchmark(SimExRead, range(1, 10))


def test_SimExRead_100(benchmark):
    result = benchmark(SimExRead, range(1, 50))


def test_SimExRead_200(benchmark):
    result = benchmark(SimExRead, range(1, 100))

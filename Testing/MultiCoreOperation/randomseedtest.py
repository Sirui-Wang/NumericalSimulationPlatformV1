import multiprocessing as mp
import time

import numpy as np


def init_worker(shared_data):
    global a, b, c
    a, b, c = shared_data


def worker(Simulation):
    global a, b, c
    np.random.seed(Simulation)
    result = np.random.normal(a, b, c)
    # print(os.getpid())
    # print(a,b,c)
    time.sleep(0.1)
    return result


def main():
    Simulations = np.arange(0, 100, 1)
    overallResult = np.zeros(100)
    with mp.Pool(processes=4, initializer=init_worker, initargs=((0, 1, 100),)) as pool:
        for result in pool.imap_unordered(worker, Simulations, chunksize=10):
            overallResult = np.add(overallResult, result)
    return overallResult


if __name__ == '__main__':
    a = main()
    b = main()
    c = main()
    print(a.sort() == b.sort())
    print(a.sort() == c.sort())
    print(a.sum())
    print(b.sum())
    print(c.sum())

import multiprocessing as mp
import time

import numpy as np
from matplotlib import pyplot as plt

import LongOperation


def worker(Inputs):
    IterRange = Inputs
    Array1 = np.random.normal(0, 1 / 10, len(IterRange))
    Array2 = np.zeros(len(IterRange))
    Array2[2] = 1
    ResultArray = LongOperation.SomeFunc((Array1, Array2))
    return ResultArray


def multi(IterRange):
    pool = mp.pool(mp.cpu_count() - 4)
    Superpositioned = pool.map(worker, IterRange)
    pool.close()
    pool.join()
    return Superpositioned


def main():
    start_time = time.time()
    IterRange = np.arange(0, 10)
    Superpositioned = np.zeros(len(IterRange) * 2 - 1)
    np.random.seed(0)
    Superpositioned = multi(IterRange, Superpositioned)
    print(Superpositioned)
    plt.figure("Multitest")
    plt.plot(Superpositioned)
    print("--- %s seconds ---" % (time.time() - start_time))
    plt.show()


if __name__ == "__main__":
    main()

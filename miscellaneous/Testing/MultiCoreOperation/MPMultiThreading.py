import multiprocessing as mp
import time
# from multiprocessing.pool import ThreadPool
from multiprocessing import Pool

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

import LongOperation


def newFunc(IterRange, Superpostioned):
    for i in tqdm(IterRange):
        Array1 = np.random.normal(0, 1 / 10, 10)
        Array2 = np.zeros(10)
        Array2[2] = 1
        ResultArray = LongOperation.SomeFunc((Array1, Array2))
        Superpostioned = np.add(Superpostioned, ResultArray)
    return Superpostioned

def main():
    start_time = time.time()
    IterRange = np.arange(0, 10)
    Superpositioned = np.zeros(len(IterRange) * 2 - 1)
    np.random.seed(0)
    mp.set_start_method("fork")
    pool = Pool(processes=2)
    threads = []
    for i in [0,5]:
        p = pool.apply_async(newFunc, (IterRange[i:i+5], Superpositioned))
        threads.append(p)
    for thread in threads:
        Superpositioned = np.add(Superpositioned, thread.get())
    print(Superpositioned)
    plt.figure("test")
    plt.plot(Superpositioned)
    print("--- %s seconds ---" % (time.time() - start_time))
    plt.show()


if __name__ == "__main__":
    main()

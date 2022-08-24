import time

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

import LongOperation


def main():
    start_time = time.time()
    IterRange = np.arange(0, 10)
    Superpositioned = np.zeros(len(IterRange) * 2 - 1)
    np.random.seed(0)
    for i in tqdm(IterRange):
        Array1 = np.random.normal(0, 1 / 10, len(IterRange))
        Array2 = np.zeros(len(IterRange))
        Array2[2] = 1
        ResultArray = LongOperation.SomeFunc((Array1, Array2))
        Superpositioned = np.add(Superpositioned, ResultArray)
    print(Superpositioned)
    plt.figure("test")
    plt.plot(Superpositioned)
    print("--- %s seconds ---" % (time.time() - start_time))
    plt.show()


if __name__ == "__main__":
    main()

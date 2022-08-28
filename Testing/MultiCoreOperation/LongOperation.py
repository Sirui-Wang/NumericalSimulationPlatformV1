import time

import numpy as np


def SomeFunc(Tuple1):
    """Takes tuple as input, tuple contains multiple numpy array"""
    Array1, Array2 = Tuple1
    ConvolutedArray = np.convolve(Array1, Array2)
    time.sleep(1)
    return ConvolutedArray


if __name__ == "__main__":
    SomeFunc()

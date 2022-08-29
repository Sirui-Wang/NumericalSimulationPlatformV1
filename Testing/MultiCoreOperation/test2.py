# SuperFastPython.com
# example of parallel map() with the process pool
from multiprocessing.pool import Pool
from random import random
from time import sleep


# task executed in a worker process
def task(identifier):
    # generate a value
    value = random()
    # report a message
    print(f'Task {identifier} executing with {value}', flush=True)
    # block for a moment
    sleep(value)
    # return the generated value
    return value


# protect the entry point
if __name__ == '__main__':
    # create and configure the process pool
    with Pool(processes=2) as pool:
        # execute tasks in order
        for result in pool.map(task, range(100)):
            print(f'Got result: {result}', flush=True)
    # process pool is closed automatically

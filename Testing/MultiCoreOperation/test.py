# SuperFastPython.com
# example of parallel map_async() with the process pool
from multiprocessing.pool import Pool
from random import random
from time import sleep


# task executed in a worker process
def task(identifier):
    # generate a value
    value = random() * 10
    # report a message
    print(f'Task {identifier} executing with {value}', flush=True)
    # block for a moment
    sleep(value)
    # return the generated value
    return value


# protect the entry point
if __name__ == '__main__':
    # create and configure the process pool
    with Pool() as pool:
        # issues tasks to process pool
        result = pool.map_async(task, range(10))
        # iterate results
        print("continued main1")
        for result in result.get():
            print(f'Got result: {result}', flush=True)
    # process pool is closed automatically
    print("Continued main")

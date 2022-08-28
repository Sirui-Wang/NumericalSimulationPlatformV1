import multiprocessing
import os
import time


def task_sleep(sleep_duration, task_number):
    time.sleep(sleep_duration)
    print(f"Task {task_number} done (slept for {sleep_duration}s)! "f"Process ID: {os.getpid()} "f"time: {time.ctime(time.time())}\n")

    return task_number

def multi_processing_startmap():
    time_start = time.time()

    # Create pool of workers
    pool = multiprocessing.Pool(3)

    tasks=[(1,i) for i in [1,2,3,4,5]]
    # Map pool of workers to process
    a=pool.starmap(func=task_sleep, iterable=tasks)
    print(a)
    # Wait until workers complete execution
    # process pool CAN closed automatically
    pool.close()

    time_end = time.time()
    print(f"Time elapsed: {round(time_end - time_start, 2)}s")

if __name__ == "__main__":
    multi_processing_startmap()

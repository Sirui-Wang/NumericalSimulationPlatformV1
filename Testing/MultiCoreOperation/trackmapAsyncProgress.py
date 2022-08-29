import multiprocessing
import time


def track_job(job, update_interval=3):
    while job._number_left > 0:
        print("Tasks remaining = {0}".format(
            job._number_left * job._chunksize))
        time.sleep(update_interval)


def hi(x):  # This must be defined before `p` if we are to use in the interpreter
    time.sleep(x // 2)
    return x


if __name__ == '__main__':
    a = [x for x in range(50)]
    p = multiprocessing.Pool()
    res = p.map_async(hi, a)
    track_job(res)

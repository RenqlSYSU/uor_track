"""
Purpose: Caculate the JJA East Asian summer monsoon index and output netcdf data in CMIP datasets

Created on August 25 2021
@author: Shan He
"""

import numpy
import xarray
from datetime import datetime
from scipy import stats
from multiprocessing import Process, Queue
import sys

def anom_dtrend(x, task_q, input_q, done_q):
    '''
    Calculate monthly anomalies and remove linear trend
    '''
    clim = x.groupby('time.month').mean('time')
    anom = x.groupby('time.month') - clim
    xvals = anom['time'].copy(data=numpy.linspace(0, anom.time.size, anom.time.size))
    while task_q.get() is not None:
        j = input_q.get()
        for i in range(anom.lat.size):
            slope, intercept, tmp, tmp, tmp = stats.linregress(xvals, anom[:,i,j])
            anom[:,i,j].values -= slope * xvals + intercept
        done_q.put((j,anom[:,:,j]))


worker_n = 50 # !!!

if __name__ == '__main__':
    
    ua = xarray.open_dataset('### NC filename ###')['ua'][:,::-1] # !!!
    task_q = Queue()
    input_q = Queue()
    done_q = Queue()
    for j in range(ua.lon.size):
        task_q.put(1)
        input_q.put(j)
    for n in range(worker_n):
        task_q.put(None)
    workers = [Process(target = anom_dtrend,
               args = (ua.copy(), task_q, input_q, done_q)) for n in range(worker_n)]
    for w in workers:
        w.start()
    t_start = datetime.now()
    t_last = t_start
    anom = ua.copy(data=numpy.empty(ua.shape))
    n_done = 0
    while n_done < ua.lon.size:
        j, anomY = done_q.get()
        anom.values[:,:,j] = anomY
        n_done += 1
        t_now = datetime.now()
        if (t_now - t_last).total_seconds() > 30:
            t_last = t_now
            dt = (t_now - t_start) / n_done * (ua.lon.size - n_done)
            print("PROGRESS: %d/%d complete, predicted completion at %s"
                  % (n_done, ua.lon.size, t_now + dt))
    for w in workers:
        w.join()
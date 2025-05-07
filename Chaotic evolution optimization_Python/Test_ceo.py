import Benchmark
import numpy as np
import matplotlib.pyplot as plt
from CEO import CEO
import time

func = Benchmark.Sphere
Dim = 100 # Dimensions of the optimization problem
Varmin = -100*np.ones([1,Dim])
Varmax = 100*np.ones([1,Dim])
MaxFES = Dim*1e4 # Maximum number of evaluation functions
N = 10 # number of chaotic samples
Np = 30 # population size  (Np is set to an even number greater than 2)

time_start = time.time()
ceo = CEO(func, Np, Dim, Varmin, Varmax, N, MaxFES)
Best, fBest, history = ceo.optimize()
time_end = time.time()
print("running time:{:.5f}".format(time_end - time_start))

plt.semilogy(history,'b-')
plt.xlabel('Iterations')
plt.ylabel('fitness')
plt.show()

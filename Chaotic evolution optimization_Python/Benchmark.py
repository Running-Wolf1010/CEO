#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from functools import reduce
from operator import mul
import numpy as np
import math


def Sphere(s):
	square = map(lambda y: y * y ,s)
	return sum(square)

def Rosenbrock(s):
	square = map(lambda x, y: 100 * (y - x**2)**2 + (x - 1)**2 , s[:len(s)-1],s[1:len(s)])
	return  sum(square)

def Rastrigin(s):
	square = map(lambda x: x**2 - 10 * math.cos(2*math.pi*x) + 10, s)
	return  float(sum(square))

def Michalewicz(s):
	n = len(s)
	t1 = -sum(map(lambda x, y: math.sin(x) * (math.sin( y*x**2/math.pi ))**20, s, range(1,n+1)))
	return t1

def Griewank(s): 
	t1 = sum(map(lambda x: 1/4000 * x**2, s))
	n = len(s)
	t2 = map(lambda x, y: math.cos(x/np.sqrt(y)), s, range(1,n+1))
	t3 = reduce(mul, t2)
	return t1 - t3 + 1

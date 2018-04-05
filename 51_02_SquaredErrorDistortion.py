# python3
import sys
from math import *
import numpy as np

'''
Squared Error Distortion Problem: Compute the squared error distortion of a set of data points with respect to a set of 
centers. 
  Input: A set of points Data and a set of centers Centers.
  Output: The squared error distortion Distortion(Data, Centers). 

Code Challenge: Solve the Squared Error Distortion Problem.
     Input: Integers k and m, followed by a set of centers Centers and a set of points Data.
     Output: The squared error distortion Distortion(Data, Centers).
Sample Input:
2 2
2.31 4.55
5.96 9.08
--------
3.42 6.03
6.23 8.25
4.76 1.64
4.47 4.33
3.95 7.61
8.93 2.97
9.74 4.03
1.73 1.28
9.72 5.01
7.27 3.77
Sample Output:
18.246

To address limitations of MaxDistance, we will introduce a new scoring function. Given a set Data of n data points and a set 
Centers of k centers, the squared error distortion of Data and Centers, denoted Distortion(Data, Centers), is defined as the 
mean squared distance from each data point to its nearest center, 

Distortion(Data,Centers) = (1/n) ∑all points DataPoint in Datad(DataPoint, Centers)2 .

Note that whereas MaxDistance(Data, Centers) only accounts for the length of the single red segment in the figure below, the 
squared error distortion accounts for the length of all segments.
'''

def readFromFile():
    f = open('input.txt', 'r')
    raw = [d.split() for d in f.read().strip().split('--------')]
    k, m = int(raw[0][0]), int(raw[0][1])
    centers = np.zeros((k, m))
    for i in range(k):
        centers[i, ] = [float(d) for d in raw[0][2+i*m:2+(i+1)*m]]
    n = len(raw[1])//m
    data = np.zeros((n, m))
    for i in range(n):
        data[i, ] = [float(d) for d in raw[1][i*m:(i+1)*m]]
    return data, centers

def getSquaredDist(p1, p2):
    return sum((p1-p2)**2)

def calDistortion(data, centers):
    n, m = data.shape
    k = centers.shape[0]
    distortion = sum([min([getSquaredDist(data[i, ], centers[j, ]) for j in range(k)]) for i in range(n)])/n
    return distortion

def solve():
    data, centers = readFromFile()
    distortion = calDistortion(data, centers)
    print(distortion)

if __name__ == '__main__':
    solve()
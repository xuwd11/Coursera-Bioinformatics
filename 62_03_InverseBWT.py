# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Inverse Burrows-Wheeler Transform Problem: Reconstruct a string from its Burrows-Wheeler transform.
     Input: A string Transform (with a single "$$" symbol).
     Output: The string Text such that BWT(Text) = Transform.

Sample Input:
TTCCTAACG$A
Sample Output:
TACATCACGT$
'''

def InverseBWT(bwt):
    l = len(bwt)
    output = [''] * l
    count = dict()
    shortcut = [0] * l
    byteStart = dict()
    for i in range(l):
        lastChar = bwt[i]
        currCount = count.get(lastChar, 0)
        shortcut[i] = currCount
        count[lastChar] = currCount + 1

    currIndex = 0
    firstCol = []
    for char, currCount in sorted(count.items(), key = lambda x:x[0]):
        firstCol += [char] * currCount
        byteStart[char] = currIndex
        currIndex += currCount

    currIndex = 0
    for i in range(l):
        output[l-i-1] = firstCol[currIndex]
        currIndex = byteStart[bwt[currIndex]] + shortcut[currIndex]
    
    return ''.join(output)

if __name__ == '__main__':
    #text = sys.stdin.read().strip()
    f = open('input.txt', 'r')
    text = f.read().strip()
    f.close()
    ori = InverseBWT(text)
    print(ori)
    f = open('result.txt', 'w')
    f.write(ori)
    f.close()
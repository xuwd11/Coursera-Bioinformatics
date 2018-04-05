# python3

import sys

'''
Solve the Change Problem. The DPChange pseudocode is reproduced below for your convenience.
     Input: An integer money and an array Coins = (coin1, ..., coind).
     Output: The minimum number of coins with denominations Coins that changes money.
Sample Input:
40
50,25,20,10,5,1
Sample Output:
2

Pseudocode:
   DPChange(money, Coins)
      MinNumCoins(0) ← 0
      for m ← 1 to money
         MinNumCoins(m) ← ∞
         for i ← 1 to |Coins|
            if m ≥ coini
               if MinNumCoins(m - coini) + 1 < MinNumCoins(m)
                  MinNumCoins(m) ← MinNumCoins(m - coini) + 1
      output MinNumCoins(money)
'''

class MinNumCoins:
    def __init__(self):
        self._input()
        print(self.DPChange(self.money, self.coins))
    
    def _input(self):
        data = sys.stdin.read().strip().split()
        self.money = int(data[0])
        self.coins = [int(c) for c in data[1].split(',')]

    def DPChange(self, money, coins):
        minNumCoins = [0]
        for m in range(1, money + 1):
            globalMin = float('inf')
            for coin in coins:
                if m >= coin:
                    currMin = minNumCoins[m-coin] + 1
                    if currMin < globalMin:
                        globalMin = currMin
            minNumCoins.append(globalMin)
        return minNumCoins[money]                       

if __name__ == "__main__":
    MinNumCoins()
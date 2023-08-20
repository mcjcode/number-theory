class FenwickTree:
    def __init__(self, size):
        self.size = size+1
        self.tree = [0]*(size+1)
    def add(self, idx, val):
        while idx <= self.size:
            self.tree[idx] += val
            idx += idx&-idx
    def prefix_sum(self, idx):
        res = 0
        while idx>0:
            res += self.tree[idx]
            idx -= idx&-idx
        return res

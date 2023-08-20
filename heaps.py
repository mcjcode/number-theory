class heapset():
    """
    A heap that also supports quick membership testing,
    position lookup and deletion of arbitrary elements by position.
    """
    
    def __init__(self):
        self.h = []
        self.hmap = {}
        
    def __bool__(self):
        """
        is there something in the heap?
        """
        return bool(self.h)
    
    def __contains__(self, elt):
        """
        is the given element in the heap
        """
        return elt in self.hmap

    def __len__(self):
        return len(self.hmap)
    
    def __getitem__(self, elt):
        """
        Return the priority of the given element.
        """
        pos = self.hmap[elt]
        return self.h[pos][0]
            
    def heappush(self, priority, elt):
        """
        Push the given element into the heap with the given priority
        """
        if elt in self.hmap:
            raise ValueError(f'pushed element {elt} already in heapset')

        self.h.append((float('inf'), elt))
        self.hmap[elt] = len(self.h)-1
        self.raise_priority(priority, elt)
        
    def raise_priority(self, new_priority, elt):
        """
        take the element and give it more priority (i.e. lower
        numerical value)
        """
        n = self.hmap[elt]
        old_priority, _ = self.h[n]
        if new_priority > old_priority:
            raise ValueError('new priority {new_priority} is not <= old priority {priority}')
        self.h[n] = new_priority, elt
        n0 = (n-1)//2
        while n and self.h[n0] > self.h[n]:
            self.hmap[self.h[n][1]] = n0
            self.hmap[self.h[n0][1]] = n
            self.h[n], self.h[n0] = self.h[n0], self.h[n]
            n = n0
            n0 = (n-1)//2
        
    def heappop(self, n=0):
        """
        pop the element at a given location in the heap (the first element
        by default.)
        """
        priority, elt = self.h[n]
        del self.hmap[elt]

        if len(self.h)==n+1: # we popped the last element of the heap
            self.h.pop()
            return priority, elt
        
        self.h[n] = self.h.pop()
        self.hmap[self.h[n][1]] = n
        
        children = [(self.h[x],x) for x in [2*n+1, 2*n+2] if x < len(self.h)]
        while children:
            hmin, n1 = min(children)
            if self.h[n] <= hmin:
                break
            self.hmap[self.h[n][1]] = n1
            self.hmap[self.h[n1][1]] = n
            self.h[n],    self.h[n1]    = self.h[n1],    self.h[n]
            n = n1
            children = [(self.h[x],x) for x in [2*n+1, 2*n+2] if x < len(self.h)]
        return priority, elt

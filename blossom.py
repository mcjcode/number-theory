# Drew from many sources in compiling this code together including CMU & Stanford
# To all those with other open source code on the internet I thank you
# Find the lowest common ancestor in the blossom tree

__all__=['max_matching']

def lca(match, base, p, a, b):
    used = [False] * len(match)
    while True:
        a = base[a]
        used[a] = True
        if match[a] == -1:
            break
        a = p[match[a]]
    while True:
        b = base[b]
        if used[b]:
            return b
        b = p[match[b]]

# Mark the path from v to the base of the blossom
def mark_path(match, base, blossom, p, v, b, children):
    while base[v] != b:
        blossom[base[v]] = blossom[base[match[v]]] = True
        p[v] = children
        children = match[v]
        v = p[match[v]]

def find_path(graph, match, p, root):
    n = len(graph)
    used = [False] * n
    p[:] = [-1] * n
    base = list(range(n))
    used[root] = True
    q = [root]
    while q:
        v = q.pop(0)
        for to in graph[v]:
            if base[v] == base[to] or match[v] == to:
                continue
            if to == root or (match[to] != -1 and p[match[to]] != -1):
                curbase = lca(match, base, p, v, to)
                blossom = [False] * n
                mark_path(match, base, blossom, p, v, curbase, to)
                mark_path(match, base, blossom, p, to, curbase, v)
                for i in range(n):
                    if blossom[base[i]]:
                        base[i] = curbase
                        if not used[i]:
                            used[i] = True
                            q.append(i)
            elif p[to] == -1:
                p[to] = v
                if match[to] == -1:
                    return to
                to = match[to]
                used[to] = True
                q.append(to)
    return -1


# Implementation of Blossom Algorithm
def max_matching(graph):
    n = len(graph)
    match = [-1] * n
    p = [0] * n
    for i in range(n):
        if match[i] == -1:
            v = find_path(graph, match, p, i)
            while v != -1:
                pv = p[v]
                ppv = match[pv]
                match[v] = pv
                match[pv] = v
                v = ppv
    # Returns number of pairs in graph
    return sum(1 for x in match if x != -1) // 2

# Implementation of Donald Knuth's 'dancing links' algorithm
# for solving exact cover problems.
#
# https://arxiv.org/pdf/cs/0011047.pdf

from dataclasses import dataclass
import numpy as np

# forward declarations:
class Column:
    pass
class Node:
    pass

@dataclass(eq=False)
class Column:
    """
    Represents a primary (mandatory) or secondary (optional) constraint.
    """
    right: Column = None
    left: Column = None
    up: Node = None
    down: Node = None
    size: int = 0
    name: str = ''


@dataclass(eq=False)
class Node:
    """
    Represents the membership of an element (a row) in a constraint (a column)
    """
    right: Node = None
    left: Node = None
    up: Node = None
    down: Node = None
    column: Column = None
    name: str = ''

def wire_dlx_from_dlx_spec(primary_names, secondary_names, rows):
    """
    primary_names: the names of the columns that must be covered (exactly once)
    secondary_names: the names of the columns that *may* be covered (at most once)
    rows: a list of rows: each row is a list of covered columns.
    """
    m = len(rows)
    names = primary_names + secondary_names
    nPrimary = len(primary_names)
    nSecondary = len(secondary_names)
    n = nPrimary + nSecondary
    names = primary_names + secondary_names
    sizes = [len([1 for row in rows if name in row]) for name in names]

    h = Node(name='h')
    columns = [ Column(*size_name) for size_name in zip(sizes, names) ]
    for j in range(nPrimary-1):
        columns[j].right = columns[j+1]
        columns[j+1].left = columns[j]
    columns[nPrimary-1].right = h
    h.left = columns[nPrimary-1]
    columns[0].left = h
    h.right = columns[0]
    top = columns.copy()

    for row in rows:
        row.sort(key = lambda x: names.index(x))
        first_in_row = True
        for name in row:
            j = names.index(name)
            node = Node(up=top[j], column=columns[j], name="")
            top[j].down = node
            top[j] = node
            if first_in_row:
                first_in_row_node = node
                first_in_row = False
            else:
                prev_node.right = node
                node.left = prev_node
            prev_node = node
        prev_node.right = first_in_row_node
        first_in_row_node.left = prev_node
    for j, node in enumerate(top):
        node.down = columns[j]
        columns[j].up = node

    for col in columns:
        size = 0
        next = col.down
        while next is not col:
            size += 1
            next = next.down
        col.size = size

    return h, columns


def wire_dlx_from_matrix(arr, names=None, secondary_columns=[]):
    m, n = arr.shape
    #
    # make the columns.
    #
    h = Node(name='h')
    
    sizes = arr.sum(axis=0).tolist()
    if not names:
        names = [ str(i) for i in range(1,n+1) ]
    columns = [ Column(size=size, name=name) for (size, name) in zip(sizes, names) ]
    primary_columns = [j for j in range(n) if j not in secondary_columns]
    for j in primary_columns[:-1]:
        columns[j].right = columns[j+1]
        columns[j+1].left = columns[j]
    columns[primary_columns[-1]].right = h
    h.left = columns[primary_columns[-1]]
    columns[primary_columns[0]].left = h
    h.right = columns[primary_columns[0]]

    # top will point to the last node n each column so far.
    top = columns.copy()
    
    for i in range(m):
        first_in_row = True
        for j in range(n):
            if arr[i,j]:
                node = Node(column=columns[j], name=(i,j))
                # wire node to the node immediately above
                node.up = top[j]
                top[j].down = node
                # now update top[j]
                top[j] = node
                if first_in_row:
                    first_in_row_node = node
                    first_in_row = False
                else:
                    prev_node.right = node
                    node.left = prev_node
                prev_node = node
        prev_node.right = first_in_row_node
        first_in_row_node.left = prev_node

    for j, node in enumerate(top):
        node.down = columns[j]
        columns[j].up = node
        
    for col in columns:
        size = 0
        next = col.down
        while next is not col:
            size += 1
            next = next.down
        col.size = size
        
    return h, columns


def cover_column(c):
    #print(f'covering column {c.name}')
    if c.left is not None:
        # this is a primary column
        c.right.left = c.left
        c.left.right = c.right
    i = c.down
    while i is not c:
        j = i.right
        while j is not i:
            j.down.up = j.up
            j.up.down = j.down
            j.column.size -= 1
            j = j.right
        i = i.down
    return


def uncover_column(c):
    #print(f'uncovering column {c.name}')
    i = c.up
    while i is not c:
        j = i.left
        while j is not i:
            j.column.size += 1
            j.down.up = j
            j.up.down = j
            j = j.left
        i = i.up
    
    if c.left is not None:
        c.right.left = c
        c.left.right = c
    return


def print_solution(O):
    for o in O:
        print(o.column.name,end=' ')
        nxt = o.right
        while nxt is not o:
            print(nxt.column.name,end=' ')
            nxt = nxt.right
        print()
    print()


def subset_from_row(r):
    subset = [r.column.name]
    elt = r.right
    while elt is not r:
        subset.append(elt.column.name)
        elt = elt.right
    return subset


def algorithmx_dlx(h, columns, O=[], print_solutions=False):
    if h.right == h: # we found a solution
        if print_solutions:
            print_solution(O)
        return 1
    ans = 0
    # choose the column with the smallest 'size'
    # and hence the smallest branching factor.
    ci = h.right
    c = ci
    while ci is not h:
        if ci.size < c.size:
            c = ci
        ci = ci.right
            
    cover_column(c)
    r = c.down
    while r != c:
        #subset = subset_from_row(r)
        #print(f'Using subset {subset}')
        O.append(r)
        j = r.right
        while j is not r:
            cover_column(j.column)
            j = j.right
        ans += algorithmx_dlx(h, columns, O, print_solutions=print_solutions)
        r = O.pop()
        c = r.column
        j = r.left
        while j is not r:
            uncover_column(j.column)
            j = j.left
        r = r.down
    
    uncover_column(c)
    return ans

def test1():
    """
    Here we solve Knuth's toy example that he uses to
    first demonstrate the dancing links idea
    """
    import numpy as np
    matrix = """0 0 1 0 1 1 0
1 0 0 1 0 0 1
0 1 1 0 0 1 0
1 0 0 1 0 0 0
0 1 0 0 0 0 1
0 0 0 1 1 0 1"""
    a = np.array([[int(xx) for xx in row.split(' ')] for row in matrix.split('\n')])
    print('  A B C D E F G')
    print(a)
    h, cols = wire_dlx_from_matrix(a, names = list('ABCDEFG'))
    ans = algorithmx_dlx(h, cols, O=[], print_solutions=True)
    print(ans)
    assert ans == 1

    
def nqueens(n, print_solutions=False):
    """
    Return the number of ways to place n queens on an
    nxn chessboard such that no queen is attacking
    another queen.
    
    Do this by dancing links - represent the problem as
    an exact cover with 'primary' and 'secondary' columns
    where 'primary' columns represent items that must be covered
    exactly once while secondary columns only need to be
    covered *at most* once.
    """
    nr = n
    nc = n
    nd = 2*n-1
    ns = 2*n-1
    
    a = np.zeros((nr*nc,nr+nc+nd+ns),dtype=int)
    subset_index = 0
    for i in range(n):
        for j in range(n):
            #
            # when we place a queen at (i,j) it
            # occupies one row, one column, one
            # (/) diagonal and one (\) anti-diagonal
            #
            a[subset_index,i                 ] = 1
            a[subset_index,nr      +j        ] = 1
            a[subset_index,nr+nc   +i+j      ] = 1
            a[subset_index,nr+nc+nd+(i-j)+n-1] = 1
            subset_index += 1

    # every row and every column has to be occupied
    #
    # but we only require that there be one or fewer
    # 
    secondary_columns = list(range(nr+nc, nr+nc+nd+ns))
    
    names =  [f'R{i}' for i in range(n)]
    names += [f'C{i}' for i in range(n)] 
    names += [f'D{i}' for i in range(2*n-1)]
    names += [f'S{i}' for i in range(2*n-1)]
    
    h, columns = wire_dlx_from_matrix(a, names, secondary_columns)
    ans = algorithmx_dlx(h, columns, print_solutions=print_solutions)
    
    return ans


def nweakqueens_rectangle(m, n, w=0):
    """
    Return the number of ways to place m 'weak' queens on an
    m x n chessboard such that no queen is attacking
    another queen.
    
    Do this by dancing links - represent the problem as
    an exact cover with 'primary' and 'secondary' columns
    where 'primary' columns represent items that must be covered
    exactly once while secondary columns only need to be
    covered *at most* once.
    """
    if m > n:
        return 0
        
    # primary columns
    nr = m          # row (-) constraints
    
    # secondary columns
    nc = n * m      # column (|) constraints
    nd = n * m      # diagonal (/) indexed by upper right hand corner
    ns = n * m      # anti-diagonal (\) indexed by upper left hand corner
    
    a = np.zeros((n*m, nr+nc+nd+ns), dtype=int)
    subset_index = 0
    for i in range(m):
        for j in range(n):
            # set the (-) row constraint
            a[subset_index, i] = 1
            #
            # set the (|) column, (/) diagonal, and 
            # (\) anti-diagonal constraints
            #
            for k in range(n-w):
                if k<=i:              a[subset_index, nr +           n*(i-k)+j    ] = 1
                if k<=i and k<=j:     a[subset_index, nr + nc +      n*(i-k)+(j-k)] = 1
                if k<=i and k+j<=n-1: a[subset_index, nr + nc + nd + n*(i-k)+(j+k)] = 1
            subset_index += 1

    # every row has to be occupied
    # but we only require that there be one or fewer
    # in the other constraint sets.
    #
    secondary_columns = list(range(nr, nr+nc+nd+ns))
    
    names =  [f'R{i}'       for i in range(m)]
    names += [f'C[{i},{j}]' for i in range(m) for j in range(n)]
    names += [f'D[{i},{j}]' for i in range(m) for j in range(n)]
    names += [f'S[{i},{j}]' for i in range(m) for j in range(n)]
    
    h, columns = wire_dlx_from_matrix(a, names, secondary_columns)
    ans = algorithmx_dlx(h, columns, print_solutions=False)
    
    return ans

def nweakqueens(n, w=0):
    return nweakqueens_rectangle(n, n, w)

def test2():
    assert nqueens(8) == 92


def problem_des_menages(n, print_solutions=False):
    seats    = ['S'+str(i) for i in range(n)]
    husbands = ['H'+str(i) for i in range(n)]

    primary_names = seats + husbands

    rows = []
    for i, seat in enumerate(seats):
        for j, husband in enumerate(husbands):
            if (j!=i) and (j!=(i-1)%n) :
                rows.append([seat, husband])

    h, columns = wire_dlx_from_dlx_spec(primary_names, [], rows)
    ans = algorithmx_dlx(h, columns, print_solutions=print_solutions)
    return ans;

def test_problem_des_menages():
    assert problem_des_menages(5)==13


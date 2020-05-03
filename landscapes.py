import csv
import math
import itertools
import re

## Exact landscapes retain all information on piecewise cuts
## For lookup speed, the list of x,y values is kept in a binary search tree
class _LinearNode:
    """ Binary search tree node that stores linear parameters to successor node. """
    def __init__(self,x,y,m=None):
        self.left, self.right = None, None
        self.x, self.y = x, y
        if not m == None:
            self.m, self.b = m, y - m*x
    def get_prev(self,root): # returns in-order previous node
        if not self.left == None:
            return self.left.get_rightmost()
        prev = None
        while not root == None:
            if self.x > root.x:
                prev = root
                root = root.right
            elif self.x < root.x:
                root = root.left
            else:
                break
        return prev
    def get_next(self,root): # returns in-order successor node
        if not self.right == None:
            return self.right.get_leftmost()
        succ = None
        while not root == None:
            if self.x < root.x:
                succ = root
                root = root.left
            elif self.x > root.x:
                root = root.right
            else:
                break
        return succ
    def get_leftmost(self): # returns leftmost node of subtree with root self
        current = self
        while not current == None:
            if current.left == None:
                break
            current = current.left
        return current
    def get_rightmost(self): # returns rightmost node of subtree with root self
        current = self
        while not current == None:
            if current.right == None:
                break
            current = current.right
        return current
    def __iter__(self): # Ordered traversal
        if self.left:
            for node in self.left:
                yield node
        yield self
        if self.right:
            for node in self.right:
                yield node
    def _pairwise_iterator(self): # Ordered traversal with consecutive pairs
        import itertools
        i, j = itertools.tee(self)
        next(j,None)
        return zip(i,j)
    
 
class LandscapeFunction:
    def __call__(self, x):
        return self.evaluate(x)
    def evaluate(self,x):
        """Returns the value of the function at input value x."""
        if (x < self.xmin) or (x > self.xmax):
            return 0
        elif self._cache == None:
            return self._evaluate(x)
        elif x in self._cache:
            return self._cache[x]
        else:
            y = self._evaluate(x)
            self._cache[x] = y
            return y
    def convert_to_fixed_lookup(self,mesh=None):
        """
        Converts to a fixed lookup landscape function.
        
        mesh (list of x values): Lookup points in the domain.
        """
        if mesh == None:
            mesh = self.get_xvalues()
        if self.__class__.__name__ == 'LandscapeFunction_Fixed_Lookup' and mesh == self.get_xvalues():
            return self
        else:
            return LandscapeFunction_Fixed_Lookup(mesh,[self.evaluate(x) for x in mesh])
    def _pairwise_iterator(self):
        """Iterates as (i, i+1)."""
        import itertools
        i, j = itertools.tee(self)
        next(j,None)
        return zip(i,j)
    def integrate(self,x_values=None):
        """
        Compute the integral of the function. If x_values are given, these will be the mesh of the integral. If not, the mesh is all corners of the landscape function.
        """
        integral = 0
        y2 = None
        if X == None:
            pwiter = self._pairwise_iterator()
        else:
            pwiter = zip(x_values,x_values[1:])
        for x1, x2 in pwiter:
            if not x1 == x2:
                if y2 == None:
                    y1 = self.evaluate(x1)
                else:
                    y1 = y2 # this enforces only n evaluations
                y2 = self.evaluate(x2)
                if not y1 == y2:
                    integral += 0.5*(x2 - x1)*(y2 + y1)
                else:
                    integral += x2*y1 - y1*x2
        return integral
    def plot(self):
        import matplotlib.pyplot as plt
        X = self.get_xvalues()
        Y = [self(x) for x in X]
        plt.plot(X,Y)        
    def __add__(self,other):
        return LandscapeFunction_LinearCombination([1.0,1.0],[self,other])
    def __sub__(self,other):
        return LandscapeFunction_LinearCombination([1.0,-1.0],[self,other])
    def __mul__(self,other):
        return LandscapeFunction_LinearCombination([float(other)],[self])
    def __truediv__(self,other):
        return LandscapeFunction_LinearCombination([1.0/float(other)],[self])
    def __matmul__(self,other): #Integral of product
        return abs(self - other).integrate()

class LandscapeFunction_Fixed_Lookup(LandscapeFunction):
    # Cannot interpolate between points
    def __init__(self, X, Y): #X does not have to be sorted
        self.xmin, self.xmax = min(X), max(X)
        self._cache = dict(zip(X,Y))
        self.__sorted_xvalues = None
    def get_xvalues(self):
        if not self.__sorted_xvalues: 
            return self._cache.keys()
        else:
            return self.__sorted_xvalues
    def _evaluate(self,x):
        # Convert to an interpolating landscape function.
        # This is inefficient. For optimal performance, do not evaluate a fixed lookup landscape function outside of its domain.
        if not self.__sorted_xvalues:
            self.__sorted_xvalues = [x for x in sorted(self._cache)]
        self = LandscapeFunction_Interpolating(self.__sorted_xvalues,[self._cache[x] for x in self.__sorted_xvalues],0)
        return self.evaluate(x)
    def __iter__(self):
        if not self.__sorted_xvalues:
            self.__sorted_xvalues = [x for x in sorted(self._cache)]
        for x in self.__sorted_xvalues:
            yield x
    def __matmul__(self,other): #add outside of domain
        integral = 0
        # Adds integral of other that is below the domain of self
        X = [x for x in other.get_xvalues() if x <= self.xmin]
        if len(X) > 1:
            integral += other.integrate(X)
        # Adds integral of other that is above the domain of self
        X = [x for x in other.get_xvalues() if x >= self.xmax]
        if len(X) > 1:
            integral += other.integrate(X)
        # Adds integral of |self-other| on the domain of self
        integral += abs(self - other).integrate(self.get_xvalues())
        return integral

class LandscapeFunction_Zero(LandscapeFunction):
    def __init__(self):
        self.xmin, self.xmax = 0, 0
        self._cache = None
    def get_xvalues(self):
        return []
    def _evaluate(self,x):
        return 0
    def __iter__(self):
        yield 0

class LandscapeFunction_Interpolating(LandscapeFunction):
    def __init__(self, X, Y, already_sorted=False, cache_values='local'):
        """
        
        cache_values controls memoization.
            - 'none'    don't cache anything. Use for when you will only evaluate between x-values with no reason to expect repeated x-values.
            - 'local'   cache value of function only on its x values. Use when function will be evaluated on its set of x-values frequently, and rarely in between. Recommended for landscape functions in ensembles.
            - 'all'     cache all x values after every evaluation.
        """
        self.root = None
        self.cache_values = cache_values
        if cache_values == 'none':
            self._cache = None
        else:
            self._cache = dict(zip(X,Y))
        self.xmin, self.xmax = min(X), max(X)
        if already_sorted:
            self.SortedListToTree(X, Y)
        else:
            self.ListToTree(X,Y)
    def insert(self,x,y,m=None,delay_update=False):
        if self.xmax < x:
            self.xmax = x
        if self.xmin > x:
            self.xmin = x
        if self.root == None:
            self.root = _LinearNode(x,y,m)
        else:
            self._insert(self.root,x,y,m,delay_update)
    def _insert(self,node,x,y,m=None,delay_update=False): 
        if x < node.x:
            if node.left == None:
                node.left = _LinearNode(x,y,m)
                if not delay_update and m == None: # update linear parameters for new node
                    node.left.m = (node.y - y)/(node.x - x)
                    node.left.b = y - node.left.m * x
                if not delay_update: # update linear parameters for node previous to new node
                    prev = node.left.get_prev(self.root)
                    if not prev == None:
                        prev.m = (y - prev.y)/(x - prev.x)
                        prev.b = prev.y - prev.m * prev.x
            else:
                self._insert(node.left,x,y,m)
        elif x > node.x:
            if node.right == None:
                node.right = _LinearNode(x,y,m)
                if not delay_update and m == None: # update linear parameters for new node
                    succ = node.right.get_next(self.root)
                    if not succ == None:
                        node.right.m = (succ.y - y)/(succ.x - x)
                        node.right.b = y - node.right.m * x
                if not delay_update: # update linear parameters for node
                    node.m = (y - node.y)/(x - node.x)
                    node.b = node.y - node.m * node.x
            else:
                self._insert(node.right,x,y,m)
        else: # if node with same x value already exists, overwrites
            if not (node.y == y) or not delay_update:
                node.y = y
                if m == None: # update linear parameters for successor node
                    succ = node.get_next(self.root)
                    if not succ == None:
                        node.m = (succ.y - y)/(succ.x - x)
                        node.b = y - node.m * x
                else:
                    node.m, node.b = m, y - m * x
                # update linear parameters for previous node
                prev = node.get_prev(self.root)
                if not prev == None:
                    prev.m = (y - prev.y)/(x - prev.x)
                    prev.b = prev.y - prev.m * prev.x
            if not self.cache_values == 'none':
                self._cache[node.x] = node.y
    def update_all(self):
        for node1,node2 in self.root._pairwise_iterator():
            node1.m = (node2.y - node1.y)/(node2.x - node1.x)
            node1.b = node1.y - node1.m * node1.x
        self.root.get_rightmost().m, self.root.get_rightmost().b = 0.0, 0.0
    def SortedListToTree(self, X, Y):
        M = [(y2 - y1)/(x2 - x1) if x1 != x2 else 0.0 for x1, x2, y1, y2 in zip(X,X[1:],Y,Y[1:])]
        self._ListToTree(list(zip(X,Y,M+[0.0])),0,len(X)-1)
    def ListToTree(self, X, Y):
        self._ListToTree([(x,y,None) for x,y in zip(X,Y)],0,len(X)-1)
        self.update_all()
    def _ListToTree(self, nodes, a, b):
        if a > b:
            return
        mid = int(a + (b - a) / 2)
        self.insert(nodes[mid][0],nodes[mid][1],nodes[mid][2],delay_update=True)
        self._ListToTree(nodes,a,mid-1)
        self._ListToTree(nodes,mid+1,b)
    def evaluate(self, x):
        if (x < self.xmin) or (x > self.xmax):
            return 0
        if not self.cache_values == 'none' and x in self._cache:
            return self._cache[x]
        return self._evaluate(x,self.root)
    def _evaluate(self, x, node): # Find largest node.x below x and return linear interpolation
        if node == None:
            return None
        if x == node.x:
            if not self.cache_values == 'none':
                self._cache[x] = node.y
            return node.y
        if x > node.x:
            y = self._evaluate(x,node.right)
            if y == None:
                y = (node.m)*x + node.b
                if self.cache_values == 'all':
                    self._cache[x] = y
            return y
        if x < node.x:
            return self._evaluate(x,node.left)
    def get_xvalues(self): #result is sorted
        return tuple(a.x for a in self.root)
    def get_xyvalues(self):
        if self._cache == None:
            X, Y = zip(*((node.x, node.y) for node in self.root))
        else:
            X, Y = list(self._cache.keys()), list(self._cache.values())
        return X, Y
    def __abs__(self): # Changes self
        y2 = None
        ell = []
        for x1, x2 in self._pairwise_iterator():
            if y2 == None:
                y1 = self.evaluate(x1)
            else:
                y1 = y2 # this enforces only n evaluations
            y2 = self.evaluate(x2)
            if y2 < 0:
                if y1 > 0:
                    x0 = (y2*x1 - y1*x2)/(y2 - y1)
                    print(x1,x0,x2)
                    print(y1,0,y2)
                    print('')
                    ell += [(x0,0)]
                    #self.insert(x0,0,None,delay_update=True)#(y1 - y2)/(x2 - x1))
                    if not self.cache_values == 'none':
                        self._cache[x0] = 0
                ell += [(x1,-y1)]
                #self.insert(x1,-y1,None,delay_update=True)#(y2 - y1)/(x2 - x1))
                if not self.cache_values == 'none':
                    self._cache[x2] = -y2
            elif y1 < 0: #and y2 >= 0
                if y2 > 0:
                    x0 = (y2*x1 - y1*x2)/(y2 - y1)
                    print(x1,x0,x2)
                    print(y1,0,y2)
                    print('')
                    #self.insert(x0,0,None,delay_update=True)#(y2 - y1)/(x2 - x1))
                    ell += [(x0,0)]
                    if not self.cache_values == 'none':
                        self._cache[x0] = 0
                ell += [(x1,-y1)]
                #self.insert(x1,-y1,None,delay_update=True)#,(y1 - y2)/(x2 - x1))
                if not self.cache_values == 'none':
                    self._cache[x1] = -y1
        for l in ell:
            self.insert(l[0],l[1],None,delay_update=True)
        self.update_all()
        return self
    def __iter__(self):
        #Iterators of LandscapeFunctions yield sorted x values
        for node in self.root:
            yield node.x

class LandscapeFunction_LinearCombination(LandscapeFunction):
    def __init__(self,coefficients,landscape_functions,no_cacheQ=None):
        self.coefficients, self.landscape_functions = coefficients, landscape_functions
        self.xmin, self.xmax = min(landscape_function.xmin for landscape_function in self.landscape_functions), max(landscape_function.xmax for landscape_function in self.landscape_functions)
        if no_cacheQ == None: # Disable this for one-time use ensembles
            self._cache = {} 
        else:
            self._cache = None
    def _evaluate(self,x):
        return sum(coefficient*(landscape_function.evaluate(x)) for coefficient, landscape_function in zip(self.coefficients,self.landscape_functions))
    def get_xvalues(self):
        return sorted(list(set(x for x in self)))
    def collapse(self):
        return LandscapeFunction_Interpolating(*zip(*((x,self.evaluate(x)) for x in self.get_xvalues())))
    def __abs__(self):
        return abs(self.collapse())
    def __iadd__(self,other):
        self.coefficients += [1]
        self.landscape_functions += [other]
        if not self._cache == None:
            for x in self._cache:
                self._cache[x] += other.evaluate(x)
        return self
    def __isub__(self,other):
        self.coefficients += [-1]
        self.landscape_functions += [other]
        if not self._cache == None:
          for x in self._cache:
              self._cache[x] -= other.evaluate(x)
        return self
    def __imul__(self,other):
        self.coefficients = [float(other)*coefficient for coefficient in self.coefficients]
        if not self._cache == None:
          for x in self._cache:
              self._cache[x] *= float(other)
        return self
    def __itruediv__(self,other):
        self.coefficients = [float(1.0/other)*coefficient for coefficient in self.coefficients]
        if not self._cache == None:
          for x in self._cache:
              self._cache[x] *= (1.0/other)
          return self         
    def __iter__(self):
        #Iterates through sorted x-values of member LandscapeFunctions
        from heapq import merge
        for x in merge(*self.landscape_functions):
            yield x

class LandscapeFunction_Product():
    def __init__(self,L1,L2):
        self.L1, self.L2 = L1, L2
    def __iter__(self):
        from heapq import merge
        for x in merge(*[self.L1,self.L2]):
            yield x
    def _pairwise_iterator(self): # iterates as (i, i+1)
        import itertools
        i, j = itertools.tee(self)
        next(j,None)
        return zip(i,j)
    def integrate(self,X=None):
        integral = 0
        y2 = None
        if X == None:
            pwiter = self._pairwise_iterator()
        else:
            pwiter = zip(X,X[1:])
        for x1, x2 in pwiter:
            if not x1 == x2:
                if y2 == None:
                    y1 = self.L1.evaluate(x1)
                    z1 = self.L2.evaluate(x1)
                else:
                    y1 = y2 # this enforces only n evaluations
                    z1 = z2
                y2 = self.L1.evaluate(x2)
                z2 = self.L2.evaluate(x2)
                if not y1 == y2:
                    if not z1 == z2:
                        integral += (1.0/6.0)*(x2 - x1)*(y1*(2*z1+z2) + y2*(z1+2*z2))
                    else:
                        integral += (1.0/6.0)*(x2 - x1)*(y1 + y2)*z1 
                else:
                    if not z1 == z2:
                        integral += (1.0/6.0)*(x2 - x1)*y1*(z1 + z2)
                    else:
                        integral += (x2 - x1)*y1*z1
        return integral
  
 

def inner(landscape_function1, landscape_function2,X=None):
    return abs(landscape_function1 - landscape_function2).integrate(X)

class Landscape:
    def __init__(self,landscape_functions):
        self.landscape_functions = landscape_functions
    def __getitem__(self,index):
        if index < len(self.landscape_functions):
            return self.landscape_functions[index]
        else:
            return LandscapeFunction_Zero()
    def __len__(self):
        return len(self.landscape_functions)
    def evaluate(self,x,maxrank=None):
        if maxrank == None:
            return [landscape_function.evaluate(x) for landscape_function in self.landscape_functions]
        elif maxrank < len(self.landscape_functions):
            return [landscape_function.evaluate(x) for landscape_function in self.landscape_functions[:maxrank]]
        else:
            return [landscape_function.evaluate(x) for landscape_function in self.landscape_functions] + [0]*(len(self.landscape_functions) - maxrank)
    def convert_to_fixed_lookup(self,mesh=None):
        """
        - WITH a mesh, returns a Landscape_Fixed_Lookup object containing LandscapeFunction_Fixed_Lookup objects all generated using the same mesh. Using a mesh is strongly recommended.
        - WITHOUT a mesh, returns a Landscape object containing LandscapeFunction_Fixed_Lookup objects which are hashed over their own possibly different xvalue sets.
        - Note that landscapes will be 0 valued outside of the mesh, so the mesh should be between xmin and xmax.
        """
        if mesh == None:
            return Landscape([landscape_function.convert_to_fixed_lookup(mesh) for landscape_function in self.landscape_functions])
        else:
            return Landscape_Fixed_Lookup([landscape_function.convert_to_fixed_lookup(mesh) for landscape_function in self.landscape_functions],mesh)
    def get_range(self):
        return (min(landscape_function.xmin for landscape_function in self.landscape_functions), max(landscape_function.xmax for landscape_function in self.landscape_functions))
    def get_mesh(self,num_bins):
        xmin, xmax = self.get_range()
        delta = (xmax - xmin)/num_bins
        return [xmin + delta*k for k in range(num_bins + 1)]
    def write(self,filename):
        LandscapeWriter.write(self,filename)
    def __add__(self,other):
        return Landscape_LinearCombination([1.0,1.0],[self,other])
    def __sub__(self,other):
        return Landscape_LinearCombination([1.0,-1.0],[self,other])
    def __mul__(self,other):
        return Landscape_LinearCombination([float(other)],[self])
    def __truediv__(self,other):
        return Landscape_LinearCombination([1.0/float(other)],[self])
    def __abs__(self):
        return Landscape([abs(landscape_function) for landscape_function in self.landscape_functions])
    def __matmul__(self,other):
        return sum([self[k] @ other[k] for k in range(max(len(self),len(other)))])

class Landscape_Fixed_Lookup(Landscape):
    #Landscape_Fixed_Lookup objects are assumed to contain LandscapeFunction_Fixed_Lookup objects all generated over the same mesh.
    def __init__(self,landscape_functions,mesh):
        self.landscape_functions = landscape_functions
        self.mesh = mesh

class Landscape_LinearCombination(Landscape):
    def __init__(self,coefficients,landscapes):
        maxdepth = max([len(landscape) for landscape in landscapes])
        self.landscape_functions = [LandscapeFunction_LinearCombination(coefficients,[landscape[i] for landscape in landscapes]) for i in range(maxdepth)]
    def collapse(self):
        return Landscape([landscape_function.collapse() for landscape_function in self.landscape_functions]) 
    def __iadd__(self,other):
        for landscape_function in self.landscape_functions:
            landscape_function += other
        return self
    def __isub__(self,other):
        for landscape_function in self.landscape_functions:
            landscape_function -= other
        return self
    def __imul__(self,other):
        for landscape_function in self.landscape_functions:
            landscape_function *= other
        return self
    def __itruediv__(self,other):
        for landscape_function in self.landscape_functions:
            landscape_function /= other
        return self

def average(collection):
    if all([issubclass(entry.__class__,Landscape) for entry in collection]):
        return Landscape_LinearCombination([float(1)/float(len(collection))] * len(collection) , collection)
    if all([issubclass(entry.__class__,LandscapeFunction) for entry in collection]):
        return LandscapeFunction_LinearCombination([float(1)/float(len(collection))] * len(collection) , collection)


class Landscape_Reader:
    def __almostequal(x,y):
        EPSILON = 0.00000000001
        return abs(x-y) <= EPSILON
    def __from_Barcode(barcode): ## Generates exact landscape from barcodes, input as [birth time, death time] pairs
        def birth(a): return a[0]-a[1]
        def death(a): return a[0]+a[1]
        landscape_functions = []
        barcode = sorted(barcode,key= lambda x:(x[0],-x[1])) # sort primarily by 1st arg ascending, secondarily by 2nd arg descending
        barcode = [((p[0]+p[1])/2.0, (p[1]-p[0])/2.0) for p in barcode] # map to center, radius form
        while not len(barcode) == 0:
            L = [(birth(barcode[0]),0.0),barcode[0]]
            i = 1
            newbarcode = []
            while i < len(barcode):
                p = 1
                if (birth(barcode[i]) >= birth(L[-1])) and (death(barcode[i]) > death(L[-1])):
                    if birth(barcode[i]) < death(L[-1]):
                        pt = ((birth(barcode[i]) + death(L[-1]))/2.0, (death(L[-1]) - birth(barcode[i]))/2.0)
                        L.append(pt)
                        while (i+p < len(barcode)) and (Landscape_Reader.__almostequal(birth(pt),birth(barcode[i+p]))) and death(pt) <= death(barcode[i+p]):
                            newbarcode.append(barcode[i+p])
                            p = p + 1
                        newbarcode.append(pt)
                        while (i+p < len(barcode)) and (birth(pt) <= birth(barcode[i+p])) and (death(pt) >= death(barcode[i+p])):
                            newbarcode.append(barcode[i+p])
                            p = p + 1
                    else:
                        L.append((death(L[-1]),0.0))
                        L.append((birth(barcode[i]),0.0))
                    L.append(barcode[i])
                else:
                    newbarcode.append(barcode[i])
                i = i + p
            L.append((death(L[-1]),0.0))
            #remove duplicates from L
            seen = set()
            seen_add = seen.add
            L = [x for x in L if not (x in seen or seen_add(x))]
            landscape_functions.append(LandscapeFunction_Interpolating(*tuple(zip(*L))))
            barcode = newbarcode 
        return Landscape(landscape_functions)
    def __from_PointLists(landscape_pointlists):
        return Landscape([LandscapeFunction_Interpolating(*tuple(zip(*pointlist))) for pointlist in landscape_pointlists])
    def __read_bar_file(filename, ERRORMAX = 10):
        """
        Reads barcode data from .bar file, assumed created by perseus.
        Set ERRORMAX slightly over intended diameter of sample space to debug problems in .bar files.
        """
        data=[]
        ERRORMAX = 10
        with open(filename,'r') as barcodefile:
            barcodereader = csv.reader(barcodefile,delimiter=' ')
            for row in barcodereader:
                b,d = [float(x) for x in row]
                if b >= 0 and d >= 0 and b < ERRORMAX and d < ERRORMAX: #throw out infinite -1, 0 barcodes and anything that extends beyond sample diameter
                    data.append([b,d])
                else:
                    if b >= 0 and d >= 0:
                        print('ERRORMAX triggered')
        return Landscape_Reader.__from_Barcode(data)
    def __read_lan_file(filename):
        """
        Reads landscape data from .lan file.
        """
        data, current = [], None
        with open(filename,'r') as landscapefile:
            for line in landscapefile:
                newlandscape = re.compile("lambda_(\\d+)")
                newlandscapematches = newlandscape.findall(line)
                if len(newlandscapematches) > 0:
                    if not current == None:
                        data.append(current)
                    current = []
                    lasta = -1.0
                    continue
                newpoint = re.compile("([-\\d\\.e]+)\\s+([-\\d\\.e]+)")
                newpointmatches = newpoint.findall(line)
                if len(newpointmatches) > 0:
                    number = re.compile("(-*\\d\\.*\\d*)e*(-*[\\d]*)")
                    numbermatch = number.findall(newpointmatches[0][0])
                    if numbermatch[0][1] == '':
                        a = float(numbermatch[0][0])
                    else:
                        a = float(numbermatch[0][0])*10**float(numbermatch[0][1])
                    numbermatch = number.findall(newpointmatches[0][1])
                    if numbermatch[0][1] == '':
                        b = float(numbermatch[0][0])
                    else:
                        b = float(numbermatch[0][0])*10**float(numbermatch[0][1])
                    if not ( a == -1.0 or a == -0.5 or Landscape_Reader.__almostequal(a,lasta) ): # This disallows duplicate entries. Adjust EPSILON as necessary.
                        current.append([a,b])
                        lasta = a
            data.append(current)
        return Landscape_Reader.__from_PointLists(data)
    def read(filename):
        if filename[-3:] == 'bar':
            return Landscape_Reader.__read_bar_file(filename)
        elif filename[-3:] == 'lan':
            return Landscape_Reader.__read_lan_file(filename)
        else:
            return Landscape_Reader.__from_PointLists(filename)

# TODO
#class Landscape_Writer: 
#    def write(landscape,filename):
#       pass
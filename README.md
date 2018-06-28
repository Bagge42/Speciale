The purpose of this algorithm is to compute a proper equilibrium of bimatrix games. 

Usage (Requires a Python installation):
Run ipython in a cmd window within the folder containing the algorithm code.
Create two payoff matrices as NumPy arrays.
Execute the method findProperEquilibrium on the two payoff matrices.

Example:
$ ipython
In [1]: import numpy as np
In [2]: import specialePythonInt
In [3]: A=np.array(([1,2],[2,5]))
In [4]: B=np.array(([0,1],[0,2]))
In [5]: print(specialePythonInt.findProperEquilibrium(A,B))
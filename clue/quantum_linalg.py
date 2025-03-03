r"""
    Module for dedicated operations related with Linear Algebra in the Quantum setup

    In this module we include all the structures and code that is related with Linear Algebra that are 
    useful to study quantum circuits and related applications.

    In general this module will complement :mod:`linalg` by extending all the necessary classes
    of :class:`clue.linalg.Vector`, :class:`clue.linalg.Matrix` and :class:`clue.linalg.Subspace`
    to be able to compute quantum bisimulations as in the papers

    * Forward and Backward Constrained Bisimulations for Quantum Circuits (https://doi.org/10.1007/978-3-031-57249-4_17)
    * Forward and Backward Constrained Bisimulations for Quantum Circuits Using Decision Diagrams (https://doi.org/10.1145/3712711)
"""
from .linalg import Vector, Matrix, SparseRowMatrix, SparseVector

from sympy.polys.domains.domain import Domain
from .numerical_domains import CC

class DensityVector(Vector):
    def __init__(self, dim: int, field: Domain = CC):
        super().__init__(dim, field)
        self.__data: list[SparseVector] = [SparseVector(self.dim, self.field) for _ in range(self.dim)]

    def reduce(self, coef, vector):
        for i in range(self.dim):
            self[i].reduce(coef, vector[i])

    def scale(self, coef):
        for i in range(self.dim):
            self[i].scale(coef)

    def conjugate(self, *, _inplace=False):
        result = self if _inplace else [self[i].copy() for i in range(self.dim)]
        for i in range(self.dim):
            result[i] = result[i].conjugate()

        return result

    def inner_product(self, rhs, *, _conjugate = True):
        lhs = self.conjugate() if _conjugate else self # we conjugate the vector (in case the field is CC) if indicated by argument
        result = self.field.zero
        for i in range(self.dim):
            result += lhs.__data[i] * rhs.__data[i]

        return result
    
    def apply_matrix(self, matr):
        if isinstance(matr, DensityOperator):
            return NotImplemented
        elif isinstance(matr, SparseRowMatrix):
            return NotImplemented
        else:
            return NotImplemented
        #if we get densityoperator we do the inner loop (should go down to sparserowmatrix)
        #if the matric that we get is a sparserowmatrix (the same as in sparserowmatrix)

    def __add__(self, other):
        if self.dim != other.dim:
            return NotImplemented
        if self.field != other.field:
            return NotImplemented
        if not isinstance(other, DensityVector):
            return NotImplemented
        
        result = DensityVector(self.dim, self.field)
        for i in range(self.dim):
            result.__data[i] = self[i] + other[i]

        return result
     
    def __getitem__(self, i: int):
        if(i < 0 or i >= self.dim):
            raise IndexError(f"Element {i} out of dimension")
        return self.__data[i]
        


class DensityOperator(Matrix):
    r'''
        Class for representing superoperators in noisy quantum circuits.

        A superoperator `A` is the combination of several quantum noisy gates that can be described as follows:
        the noisy gate `((p_i, U_i))` applies to a quantum state the gate `U_i` with probability `p_i`.

        These quantum noisy gates can be represented with a matrix `A_i` that works over the space of density matrices
        (for `n` qbits, there are `N=2^n` quantum states and `2^{2n} = N^2` density matrices). Hence, these superoperator
        are matrices of dimension `N^2`.

        At the end of the day, when we combine several gates, we still get a set `((\pi_i, C_i))` where we get to apply
        the full circuit `C_i` with probability `\pi_i`.
    '''
    def __init__(self, circuits: tuple[SparseRowMatrix], probabilities: tuple, dim:int = None):
        # Same length of two arguments
        if len(circuits) != len(probabilities):
            raise ValueError(f"`circuits` and `probabilities` must be tuples of same length")
        if len(circuits) == 0: # no circuits - identity case - we use dimension
            if dim is None:
                raise ValueError(f"Identity matrix without dimension")
            super().__init__(dim, CC)
            self.__data = tuple()
        else:
            ## The circuits must have all the same dimension
            if not all(c.dim == circuits[0].dim for c in circuits[1:]):
                raise TypeError("We have different circuits in each probability")
            if any(not c.is_square() for c in circuits):
                raise TypeError(f"A circuit must always be a square matrix")
            
            N = circuits[0].nrows
            if dim is not None and dim != N**2:
                raise ValueError(f"Dimension provided with circuits is not compatible")
            
            super().__init__(N**2, CC)

            self.__data = tuple(zip(circuits,probabilities))
    
    def data(self):
        return self.__data

    @classmethod
    def eye(cls, dim: int):
        return cls(dim=dim)
    
        # * ``eye``: class method to create identity matrix
        # * ``transpose``: method to create the transpose of a matrix
        # * ``conjugate``: method to compute the conjugate (entry-wise) of a matrix
        # * ``_add_matrix_``: receives another matrix and computes the addition
        # * ``_add_matrix_inplace_``: same as before, but do computations inplace
        # * ``_matmul_``: performs matrix multiplication.
        # * ``scalar``: scales a matrix using a scalar number
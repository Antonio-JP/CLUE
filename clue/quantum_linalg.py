from __future__ import annotations
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

    @staticmethod
    def from_matrix(matrix: SparseRowMatrix) -> DensityVector:
        if not matrix.is_square():
            raise TypeError(f"DensityVectors are always square matrices")
        
        output = DensityVector(matrix.dim[0], matrix.field)
        output.__data = [matrix[i].copy() for i in range(matrix.nrows)] # we override the rows of the matrix

        return output

    @staticmethod
    def from_vector(vector: SparseVector) -> DensityVector:
        return DensityVector.from_matrix(vector.tensor(vector))
    
    @staticmethod
    def from_ensemble(vectors: tuple[SparseVector], probabilities: tuple[float]) -> DensityOperator:
        if len(vectors) <= 0 or len(vectors) != len(probabilities):
            raise TypeError(f"The input must be non-empty lists of same lengths")
        if sum(probabilities) != 1:
            raise ValueError(f"The probabilities must provide a valid finite distribution (i.e., add up to 1)")
        return sum(p*DensityVector.from_vector(v) for (p,v) in zip(vectors, probabilities))

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
            if matr.is_ensembled(): # base case -> sum of probabilities * apply circuits
                return NotImplemented # TODO: by Thomas
            else: # composed case -> we apply one by one
                v = self
                for operator in matr.operators():
                    v = v.apply_matrix(operator)
                return v
        elif isinstance(matr, SparseRowMatrix):
            return NotImplemented # TODO: by Thomas
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
        Class for representing super-operators in noisy quantum circuits.

        A super-operator `A` is the combination of several quantum noisy gates that can be described as follows:
        the noisy gate `((p_i, U_i))` applies to a quantum state the gate `U_i` with probability `p_i`.

        These quantum noisy gates can be represented with a matrix `A_i` that works over the space of density matrices
        (for `n` qbits, there are `N=2^n` quantum states and `2^{2n} = N^2` density matrices). Hence, these super-operator
        are matrices of dimension `N^2`.

        At the end of the day, when we combine several gates, we still get a set `((\pi_i, C_i))` where we get to apply
        the full circuit `C_i` with probability `\pi_i`.
    '''
    def __init__(self, *,
                circuits: tuple[SparseRowMatrix] = None, probabilities: tuple = None, 
                operators : tuple[DensityOperator] = None,
                dim:int = None):
        self.__data = None
        self.__operators = None
        # We have three options to create a density operator:
        ## it is a ensemble operator --> given by a tuple of circuits and probabilities
        if circuits is not None and probabilities is not None:
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
        elif circuits is not None or probabilities is not None:
            raise ValueError(f"Either both or none 'circuits' and 'probabilities' are provided.")
        elif operators != None:
            if any(not isinstance(op, DensityOperator) or not op.is_ensembled() for op in operators):
                raise ValueError(f"Composite operator: must have as pieces all ensembled density operators")
            elif any(op.dim != operators[0].dim for op in operators):
                raise TypeError(f"Composite operator: all operators must have the same dimension")
            
            super().__init__(operators[0].dim, CC)
            self.__operators = operators
        else:
            raise ValueError(f"Density Operator: incompatible input for class")
    
    def data(self):
        return self.__data

    def operators(self) -> tuple[DensityOperator]:
        r'''
            Return a tuple of ensembled density operators that represent self
        '''
        if self.is_ensembled():
            return (self,)
        return self.__operators
    
    ## Methods for Density Operators
    def is_ensembled(self) -> bool:
        return self.__operators is None
    
    def is_identity(self) -> bool:
        return self.__data is not None and len(self.__data) == 0

    ## Abstract methods from Matrix
    @classmethod
    def eye(cls, dim: int):
        return cls(dim=dim)
    
    def transpose(self) -> DensityOperator:
        if self.is_ensembled():
            circuits, probabilities = list(zip(*self.data()))
            return DensityOperator(circuits=tuple(M.transpose() for M in circuits), probabilities=probabilities)
        else:
            return DensityOperator(operators=tuple(op.transpose() for op in self.operators()[::-1]))
        
    def conjugate(self) -> DensityOperator:
        if self.is_ensembled():
            circuits, probabilities = list(zip(*self.data()))
            return DensityOperator(circuits=tuple(M.conjugate() for M in circuits), probabilities=probabilities)
        else:
            return DensityOperator(operators=tuple(op.conjugate() for op in self.operators()[::-1]))
        
    def _add_matrix_(self, other):
        if not self.is_ensembled():
            raise TypeError(f"Adding Density Operators not valid for not ensembled case")
        elif not isinstance(other, DensityOperator) or not other.is_ensembled():
            raise TypeError(f"Adding Density Operators not valid for not ensembled case")
        
        ## Both are ensembled
        self_circ, self_prob = list(zip(*self.data()))
        other_circ, other_prob = list(zip(*other.data()))

        return DensityOperator(circuits=self_circ + other_circ, probabilities=self_prob+other_prob)
    
    def _add_matrix_inplace_(self, other):
        if not self.is_ensembled():
            raise TypeError(f"Adding Density Operators not valid for not ensembled case")
        elif not isinstance(other, DensityOperator) or not other.is_ensembled():
            raise TypeError(f"Adding Density Operators not valid for not ensembled case")

        self.__data += other.data()

    def _matmul_(self, other: DensityOperator):
        # This is how actually we multiply two operators
        if not isinstance(other, DensityOperator):
            raise TypeError(f"The composition of Density operators are only valid for other density operators")
        if self.is_identity():
            return other
        elif other.is_identity():
            return self
        else:
            return DensityOperator(operators=self.operators()+other.operators())
        
    def scalar(self, other) -> DensityOperator:
        if self.is_ensembled():
            circuits, probabilities = list(zip(*self.data()))
            return DensityOperator(circuits=tuple(M.scale(other) for M in circuits), probabilities=probabilities)
        else:
            return DensityOperator(operators=tuple(op.scale(other) for op in self.operators()))

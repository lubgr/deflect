package bvp

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type cholesky struct{}

func (c *cholesky) SolveLinearSystem(a mat.Symmetric, b, x *mat.VecDense) (result error) {
	// Gonum panics when the input matrix is singular. Since we'd like to try and handle such a case
	// more gracefully, we turn this into an error instead.
	defer func() {
		if r := recover(); r != nil {
			result = fmt.Errorf("could not Cholesky-factorise coefficient matrix: %v", r)
		}
	}()

	var ch mat.Cholesky

	if ok := ch.Factorize(a); !ok {
		det := ch.Det()
		return fmt.Errorf("failed to compute Cholesky factorisation, deteterminant = %v", det)
	}

	if err := ch.SolveVecTo(x, b); err != nil {
		det := ch.Det()
		return fmt.Errorf("failed to solve Cholesky-factorised system, det = %v", det)
	}

	return nil
}

// NewCholeskySolver creates a Cholesky solver for symmetric positive definite coefficient matrices.
func NewCholeskySolver() EquationSolver {
	return &cholesky{}
}

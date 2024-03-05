package deflect

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

func applyNodalBC(bc NodalValue, indices map[Index]int, dst *mat.VecDense) error {
	index, found := indices[bc.Index]

	if !found {
		return fmt.Errorf("index %v/%v not found in lookup map", bc.NodalID, bc.Dof)
	}

	dst.SetVec(index, bc.Value)

	return nil
}

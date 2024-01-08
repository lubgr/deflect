package integr

import (
	"encoding/json"
	"errors"
	"path/filepath"
	"sync"
	"testing"
	"unsafe"

	"github.com/google/go-jsonnet"
	"github.com/lubgr/deflect/bvp"
)

func TestIntegration(t *testing.T) {
	mutex := sync.Mutex{}
	vm := jsonnet.MakeVM()

	for _, file := range allTestFiles("*.jsonnet", t) {
		t.Run(filepath.Base(file), func(t *testing.T) {
			t.Parallel()

			mutex.Lock()
			testJSON, err := vm.EvaluateFile(file)
			mutex.Unlock()

			if err != nil {
				t.Fatalf("Failed to evaluate '%v': %v", file, err)
			}

			for _, input := range jsonAsMultiPart(testJSON, t) {
				name := extractName(input, t)
				t.Run(name, func(t *testing.T) {
					runTestCase(input, t)
				})
			}
		})
	}
}

func allTestFiles(pattern string, t *testing.T) []string {
	t.Helper()

	files, err := filepath.Glob(pattern)

	if err != nil {
		t.Fatalf("Failed to glob '%v': %v", pattern, err)
	} else if len(files) == 0 {
		t.Fatalf("Expected to glob one or more test files with '%v', got none", pattern)
	}

	return files
}

// jsonAsMultiPart splits the given input JSON string into multiple JSON strings when it is an
// array (one string per array item). Otherwise, it returns input as a byte slice.
func jsonAsMultiPart(input string, t *testing.T) [][]byte {
	t.Helper()

	if input == "" {
		return nil
	}

	var multi []json.RawMessage
	if err := json.Unmarshal([]byte(input), &multi); err != nil {
		// Not an array, so return as is.
		return [][]byte{unsafe.Slice(unsafe.StringData(input), len(input))}
	}

	result := make([][]byte, len(multi))

	for i, data := range multi {
		// Calling MarshalJSON is probably unnecessary since json.RawMessage is []byte under the hood.
		// But it shouldn't hurt, so let's not do clever shortcuts.
		if bytes, err := data.MarshalJSON(); err != nil {
			t.Fatalf("Failed to marshal unparsed JSON: %v", err)
		} else {
			result[i] = bytes
		}
	}

	return result
}

func extractName(input []byte, t *testing.T) string {
	var tmp struct {
		Name *string
	}

	// Wasteful to parse this only to extract the description, but let's not worry unless it shows up
	// in a profiler.
	if err := json.Unmarshal(input, &tmp); err != nil {
		t.Fatalf("Couldn't extract description from JSON: %v", err)
	} else if tmp.Name == nil {
		return ""
	}

	return *tmp.Name
}

func runTestCase(input []byte, t *testing.T) {
	// Parsing is split into two phases. That's not the most efficient way, but it allows a cleaner
	// separation of concerns. Parsing of expectations happens only in tests anyhow.
	problem, errProblem := bvp.FromJSON(input)
	expect, errExpect := ExpectFromJSON(input)

	if err := errors.Join(errProblem, errExpect); err != nil {
		t.Fatalf("Failed to build BVP/expectations from JSON: %v", err)
	}

	layout, err := bvp.NewIndexLayout(&problem)

	if err != nil {
		for _, e := range expect {
			e.Failure(err, t)
		}
		return
	}

	solver := bvp.NewLinearProblemSolver()
	primary, secondary, err := solver.Solve(&problem, layout, bvp.NewCholeskySolver())

	for _, e := range expect {
		e.Failure(err, t)
		e.Primary(primary, t)
		e.Reaction(secondary, t)
		e.Interpolation(t)
	}
}

package deflect

// transform implements the classical map operation, similar to Python's map built-in.
// While 'Map' would be a better name, it is too close to Go's 'map' keyword, so we
// choose Transform as a compromise and in honour of C++'s tradition of poorly named
// standard library functions.
func transform[S ~[]E, E any, T any](op func(E) T, s S) []T {
	mapped := make([]T, 0, len(s))

	for _, item := range s {
		mapped = append(mapped, op(item))
	}

	return mapped
}

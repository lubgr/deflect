
.PHONY: default
default: build

.PHONY: build
build:
	go build ./...

.PHONY: fmt
fmt:
	goimports -d . | patch

.PHONY: lint
lint:
	pre-commit run --all-files

.PHONY: test
test:
	go test -race ./...

.PHONY: wasm
wasm:
	GOOS=wasip1 GOARCH=wasm go test -exec 'wazero run -mount=.:/:ro' ./...

.PHONY: coverage
coverage:
	go test -coverprofile=/tmp/coverage.out ./...
	go tool cover -html=/tmp/coverage.out

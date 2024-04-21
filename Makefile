
.PHONY: default
default: build

.PHONY: build
build:
	go build ./...

.PHONY: generate
generate:
	go generate ./...

.PHONY: lint
lint:
	pre-commit run --all-files

.PHONY: pre-commit
pre-commit:
	# If this succeeds, running git commit --no-verify should be okay
	./scripts/clean_staged_pre_commit.sh

.PHONY: test
test:
	go test -race -parallel 8 ./...

.PHONY: wasm
wasm:
	GOOS=wasip1 GOARCH=wasm go test -exec 'wazero run -mount=.:/:ro' ./...

.PHONY: coverage
coverage:
	go test -coverprofile=/tmp/coverage.out ./...
	go tool cover -html=/tmp/coverage.out

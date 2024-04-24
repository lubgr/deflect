#!/usr/bin/env bash
#
# This script encapsulates a pre-commit run with staged changes but without untracked files. This
# solves an issue where untracked files cause the pre-commit run to yield false positive or false
# negative results. If this script succeeds, running git commit with the --no-verify option should
# be fine.
#
# A different implementation could use git stash -u before running pre-commit hooks. However, this
# needs special care for the git stash pop operation at the end, since the initial stash could be
# empty. To not mess with a potentially existing stash, this script prefers to act on a clean
# checkout.

set -x
set -e

dir=`mktemp -d`

trap "rm -rf ${dir}" EXIT

git clone . "${dir}"
git diff --staged | git -C "${dir}" apply --index --allow-empty

pushd ${dir}

pre-commit run $@

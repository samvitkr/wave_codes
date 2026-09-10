#!/usr/bin/bash

set -e -o pipefail

export NINJA_STATUS="[ninja][%f/%t] "
echo '```' >> "$GITHUB_STEP_SUMMARY"
trap 'printf "%s\n" "\`\`\`" >> "$GITHUB_STEP_SUMMARY"' EXIT
"${@}" | tee >(grep -v '^\[ninja\]\[' >> "$GITHUB_STEP_SUMMARY")

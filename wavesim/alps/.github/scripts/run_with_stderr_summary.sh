#!/usr/bin/bash

set -e -o pipefail

echo '```' >> "$GITHUB_STEP_SUMMARY"
trap 'printf "%s\n" "\`\`\`" >> "$GITHUB_STEP_SUMMARY"' EXIT
"${@}" 2> >(tee -a "$GITHUB_STEP_SUMMARY" >&2)
#!/usr/bin/env bash
set -euo pipefail

script_dir=$(dirname "${BASH_SOURCE[0]}")
tmp_dir=$(mktemp -d "$script_dir/.StopTime.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

for exponent in 1 4 7; do
    for algorithm in 1 2 3; do
        output="$tmp_dir/Alg${algorithm}-StopTime-1e${exponent}.dat"
        : > "$output"

        for m in 50 100 200 400 800 1600; do
            input="$script_dir/../$m/Alg${algorithm}-TimeError-1e${exponent}.dat"
            if [[ ! -f "$input" ]]; then
                printf 'Missing input file: %s\n' "$input" >&2
                exit 1
            fi

            if ! elapsed=$(awk 'NF { value = $1; found = 1 }
                END { if (!found) exit 1; print value }' "$input"); then
                printf 'No data rows found in: %s\n' "$input" >&2
                exit 1
            fi

            printf '%s\t%s\n' "$m" "$elapsed" >> "$output"
        done
    done
done

for output in "$tmp_dir"/*.dat; do
    mv "$output" "$script_dir/"
done

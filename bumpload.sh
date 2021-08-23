#!/usr/bin/env bash
set -uo pipefail

show_help() {
    cat <<EOF
Usage: {$0##*/} <version_part>
EOF
}

die() {
    printf '%s\n' "$1" >&2
    exit 1
}

if [[ -n "$1" ]]; then
    bumpversion --list --allow-dirty "$1" || die "bumpversion failed"
    poetry build
    poetry publish -r local
else
    show_help
fi

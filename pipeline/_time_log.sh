#!/usr/bin/env bash
set -euo pipefail

mode=$1

case "$mode" in
    start)
        time_file=$2
        rule_name=$3
        threads=${4:-}

        {
            printf 'rule=%s\n' "$rule_name"
            printf 'host=%s\n' "$(hostname)"
            printf 'cwd=%s\n' "$(pwd)"
            if [[ -n "$threads" ]]; then
                printf 'threads=%s\n' "$threads"
            fi
            printf 'start_utc=%s\n' "$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
        } > "$time_file"
        ;;
    end)
        time_tmp=$2
        time_final=$3

        printf 'end_utc=%s\n' "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" >> "$time_tmp"
        mv "$time_tmp" "$time_final"
        ;;
    *)
        printf 'usage: %s <start|end> ...\n' "$0" >&2
        exit 1
        ;;
esac

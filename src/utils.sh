#!/usr/bin/env bash
# Bash utilities

set -eu

function log(){
    local program=${0##*/}
    local now=$(date +'%Y-%m-%dT%H:%M:%S:%3N')
    echo "[${now}] (${program}) $@" >&2
}
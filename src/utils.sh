#!/usr/bin/env bash
# Bash utilities

set -eu

function log(){
    local program=${0##*/}
    local now=$(date +'%Y-%m-%dT%H:%M:%S:%3N')
    echo -e "[${now}] (${program}) $@" >&2
}
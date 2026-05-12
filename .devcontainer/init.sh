#!/usr/bin/env bash
# init.sh — runs on the HOST before the container starts.
# Ensures bind-mount sources exist so Docker doesn't create them as root-owned directories.

set -e

mkdir -p "$HOME/.claude"
touch "$HOME/.zsh_history"
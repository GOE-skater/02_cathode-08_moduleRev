#!/bin/bash

cd "$(dirname "$0")" || exit 1

if [ -z "$1" ]; then
    echo "Usage: ./ginit.command <github-repository-name>"
    echo "Example: ./ginit.command 02_cathode-07_modularization"
    exit 1
fi

REPO_NAME="$1"
GITHUB_USER="GOE-skater"

REMOTE_URL="https://github.com/${GITHUB_USER}/${REPO_NAME}.git"

echo "Local directory:"
pwd
echo

echo "GitHub repository:"
echo "$REMOTE_URL"
echo

if [ ! -d ".git" ]; then
    echo "Initializing git repository..."
    git init
fi

if git remote get-url origin >/dev/null 2>&1; then
    echo "origin already exists."
    echo "Current origin:"
    git remote get-url origin
    echo
    echo "Updating origin to:"
    echo "$REMOTE_URL"
    git remote set-url origin "$REMOTE_URL"
else
    echo "Adding origin:"
    git remote add origin "$REMOTE_URL"
fi

git branch -M main

git add .

if git diff --cached --quiet; then
    echo "No changes to commit."
else
    git commit -m "Initial commit"
fi

git push -u origin main

echo
echo "Done."
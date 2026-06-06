#!/bin/bash

cd "$(dirname "$0")" || exit 1

COMMENT="$1"

if [ -z "$COMMENT" ]; then
    COMMENT="Update files"
fi

echo "Local directory:"
pwd
echo

if [ ! -d ".git" ]; then
    echo "Error: This directory is not a git repository."
    exit 1
fi

REMOTE_URL=$(git remote get-url origin 2>/dev/null)

if [ -z "$REMOTE_URL" ]; then
    echo "Error: origin is not set."
    echo "Run ginit.command first."
    exit 1
fi

echo "GitHub remote:"
echo "$REMOTE_URL"
echo

echo "Commit message:"
echo "$COMMENT"
echo

echo "Git status before commit:"
git status
echo

git add .

if git diff --cached --quiet; then
    echo "No changes to commit."
else
    git commit -m "$COMMENT"
fi

git push

echo
echo "Done."
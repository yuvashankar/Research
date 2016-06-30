#!/bin/bash
#this is a small script that should update git automatically without me having to pull add commit push"
# how you run it is: sh gitupdate.sh "branch-name" "commit-message"

git pull origin $1

git add *

git commit -m "$2"

git push origin $1



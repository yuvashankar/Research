#!/bin/bash
#this is a small script that should update git automatically without me having to pull add commit push"
# how you run it is: sh gitupdate.sh "branch-name"
#It will then prompt you for a commit message.. type that in and then it will do the rest

git pull origin $1

git add *

#ask user for inptu
echo "What's your commit message, no quotes needed"
read input_variable

git commit -m "$input_variable"

git push origin $1



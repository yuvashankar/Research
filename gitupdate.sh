#!/bin/bash
#this is a small script that should update git automatically without me having to pull add commit push"

git pull origin master

git add *

git commit -m "$1"

git push origin master



#!/bin/sh

LID=`date +%y%m%d`

#scp lzhenn@hqlx27.ust.hk:/home/lzhenn/cooperate/script/2107-case-slp-uv-spatial-pillow.py /home/users/qd201969/uor_track/script
scp /home/users/qd201969/ERA5-1HR-lev/* qd201969@login2.jasmin.ac.uk:/home/users/qd201969/ERA5-1HR-lev

echo "Mission1.1: add project file..."
cd ~/uor_track

find . -name "*.ncl" | xargs git add
find . -name "*.sh" | xargs git add
find . -name "*.f90" | xargs git add
find . -name "*.F90" | xargs git add
find . -name "*.F" | xargs git add
#find . -name "*.vbs" | xargs git add
find . -name "*.py" | xargs git add
find . -name "*.txt" | xargs git add
find . -name "*.csh" | xargs git add
find . -name "*.md" | xargs git add
find . -name "*.m" | xargs git add
find . -name "*.in" | xargs git add

#git add */script/*
#git add */SourceMods*

echo "Mission1.2: remove project's deleted file..."
git status | grep 'deleted' | cut -d ':' -f 2 | xargs git rm --cached

git commit -m "${LID}"
git push --force origin master


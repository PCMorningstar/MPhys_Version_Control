#!/usr/bin/bash

if [ $# -eq 0 ]
  then
    echo "No arguments provided"
    exit
fi

echo Renaming "MyCustomFrame" to "${1}"

## replace words in the files

sed -i "s/MyCustomFrame/$1/g" Root/LinkDef.h
sed -i "s/MyCustomFrame/$1/g" Root/MyCustomFrame.cc
sed -i "s/MyCustomFrame/$1/g" CMakeLists.txt
sed -i "s/MyCustomFrame/$1/g" MyCustomFrame/MyCustomFrame.h

## rename the files and folders
mv Root/MyCustomFrame.cc Root/${1}.cc
mv MyCustomFrame/MyCustomFrame.h MyCustomFrame/${1}.h
mv MyCustomFrame ${1}

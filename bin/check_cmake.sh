#!/bin/bash

bad=0

function check_directory() {
  dir=$1

  if [ -d $dir/src ]; then
    pushd $dir >/dev/null 2>&1

    #echo "Checking directory $dir"
    files=($(ls src/*.cc src/*.cpp 2>/dev/null))
    cmakes=($(grep -oE "src/.*.c(c|pp)" CMakeLists.txt))

    goodfiles=(`echo ${files[@]} ${cmakes[@]} | tr ' ' '\n' | sort | uniq -d`)
    badfiles=(`echo ${files[@]} ${goodfiles[@]} | tr ' ' '\n' | sort | uniq -u`)
    badcmakes=(`echo ${goodfiles[@]} ${cmakes[@]} | tr ' ' '\n' | sort | uniq -u`)

    for item in ${badfiles[@]}; do
      echo "File $dir/$item is not defined in CMakeLists.txt!"
      bad=1
    done

    for item in ${badcmakes[@]};do
      echo "File $dir/$item defined in CMakeLists.txt but not found in directory!"
      bad=1
    done

    popd >/dev/null 2>&1
  fi
}

for dir in $PWD/*;do
  check_directory ${dir%/}
done

exit $bad

for dir in */;do
    cd $dir
    dirname=`echo $dir|sed 's|/||g'`

    if [ -f CMakeLists.txt ];then
      mv CMakeLists.txt{,.bak}
    fi

    if [ -d src ]; then

    make_lib=0
    for file in src/*.cc;do
      if [[ $file =~ "_module.cc" ]]; then
        continue
      fi
      make_lib=1
      break
    done

    if [ $make_lib -eq 1 ];then
      echo 'cet_make_library(
      SOURCE' >CMakeLists.txt

      for file in src/*.cc;do
        if [[ $file =~ "_module.cc" ]]; then
          continue
        fi
        echo $file >>CMakeLists.txt
      done

      echo "LIBRARIES PUBLIC" >>CMakeLists.txt
      echo "cetlib::cetlib" >>CMakeLists.txt # For example
      echo ")" >>CMakeLists.txt
      echo >>CMakeLists.txt
    fi

    modcount=`ls src/*_module.cc 2>/dev/null |wc -l`
    if [ $modcount -gt 0 ]; then
      for file in src/*_module.cc;do
        modname=`echo $file|sed 's|src/||g;s/_module.cc//g'`
        echo "cet_build_plugin($modname art::module REG_SOURCE $file)" >>CMakeLists.txt
        echo >> CMakeLists.txt
      done
    fi

    echo >>CMakeLists.txt
    echo "install_source(SUBDIRS src)" >> CMakeLists.txt

    fi
    if [ -d inc ]; then
      echo "install_headers(SUBDIRS inc)" >> CMakeLists.txt
    fi
    if [ -d fcl ]; then
      echo "install_fhicl(SUBDIRS fcl SUBDIRNAME $dirname)" >> CMakeLists.txt
    fi

    cd ..
done

for dir in */;do
    cd $dir
    dirname=`echo $dir|sed 's|/||g'`

    if [ -f CMakeLists.txt ];then
      mv CMakeLists.txt{,.bak}
    fi
    
    lib_list=''
    for file in src/*.cc inc/*.hh;do
      if ! [ -f $file ]; then continue; fi
      if [[ $file =~ _.*\.cc ]]; then
        continue
      fi
      if [[ $file =~ _.*\.hh ]]; then
        continue
      fi
      lib_list="$lib_list`grep Offline/ $file|grep -v $dir|sed 's|.*Offline/|offline::|g;s|/.*||g'|uniq|tr '\n' ','`"
    done
    lib_list=`echo $lib_list|tr ',' '\n'|sort|uniq|awk '{print "      "$1}'`

    make_lib=0
    make_ilib=0
    for file in src/*.cc;do
      if ! [ -f $file ]; then continue; fi
      if [[ $file =~ _.*\.cc ]]; then
        continue
      fi
      make_lib=1
      break
    done
    for file in inc/*.hh;do
      if ! [ -f $file ]; then continue; fi
      make_ilib=1
      break
    done

    if [ $make_lib -eq 1 ];then
      echo 'cet_make_library(
    SOURCE' >CMakeLists.txt

      for file in src/*.cc;do
        if [[ $file =~ _.*\.cc ]]; then
          continue
        fi
        echo "      $file" >>CMakeLists.txt
      done

      echo "    LIBRARIES PUBLIC" >>CMakeLists.txt
      echo "$lib_list" >>CMakeLists.txt
      echo ")" >>CMakeLists.txt
      echo >>CMakeLists.txt

    elif [ $make_ilib -eq 1 ];then
      echo 'cet_make_library(INTERFACE
    SOURCE' >>CMakeLists.txt
      
      for file in inc/*.hh;do
        echo "      $file" >>CMakeLists.txt
      done
      
      echo "    LIBRARIES INTERFACE" >>CMakeLists.txt
      echo "$lib_list" >>CMakeLists.txt
      echo ")" >>CMakeLists.txt
      echo >>CMakeLists.txt

    fi
    
    execcount=`ls src/*_main.cc 2>/dev/null|wc -l`
    if [ $execcount -gt 0 ]; then
      for file in src/*_main.cc;do
        execname=`echo $file|sed 's|src/||g;s/_main.cc//g'`
        echo "cet_make_exec(NAME $execname
    SOURCE $file
    LIBRARIES" >>CMakeLists.txt
        if [ $make_lib -eq 1 ] || [ $make_ilib -eq 1 ]; then
          echo "      offline::$dirname" >>CMakeLists.txt
        fi
        file_lib_list="`grep Offline/ $file|grep -v $dir|sed 's|.*Offline/|offline::|g;s|/.*||g'|uniq|tr '\n' ','`"
        file_lib_list=`echo $file_lib_list|tr ',' '\n'|sort|uniq|awk '{print "      "$1}'`
        echo "$file_lib_list" >>CMakeLists.txt
        echo ")" >>CMakeLists.txt
        echo >> CMakeLists.txt
      done
    fi

    svccount=`ls src/*_service.cc 2>/dev/null|wc -l`
    if [ $svccount -gt 0 ]; then
      for file in src/*_service.cc;do
        svcname=`echo $file|sed 's|src/||g;s/_service.cc//g'`
        echo "cet_build_plugin($svcname art::service
    REG_SOURCE $file
    LIBRARIES REG" >>CMakeLists.txt
        if [ $make_lib -eq 1 ] || [ $make_ilib -eq 1 ]; then
          echo "      offline::$dirname" >>CMakeLists.txt
        fi
        file_lib_list="`grep Offline/ $file|grep -v $dir|sed 's|.*Offline/|offline::|g;s|/.*||g'|uniq|tr '\n' ','`"
        file_lib_list=`echo $file_lib_list|tr ',' '\n'|sort|uniq|awk '{print "      "$1}'`
        echo "$file_lib_list" >>CMakeLists.txt
        echo ")" >>CMakeLists.txt
        echo >> CMakeLists.txt
      done
    fi
    
    srccount=`ls src/*_source.cc 2>/dev/null |wc -l`
    if [ $srccount -gt 0 ]; then
      for file in src/*_module.cc;do
        srcname=`echo $file|sed 's|src/||g;s/_source.cc//g'`
        echo "cet_build_plugin($modname art::source
    REG_SOURCE $file
    LIBRARIES REG" >>CMakeLists.txt
        if [ $make_lib -eq 1 ] || [ $make_ilib -eq 1 ]; then
          echo "      offline::$dirname" >>CMakeLists.txt
        fi
        file_lib_list="`grep Offline/ $file|grep -v $dir|sed 's|.*Offline/|offline::|g;s|/.*||g'|uniq|tr '\n' ','`"
        file_lib_list=`echo $file_lib_list|tr ',' '\n'|sort|uniq|awk '{print "      "$1}'`
        echo "$file_lib_list" >>CMakeLists.txt
        echo ")" >>CMakeLists.txt
        echo >> CMakeLists.txt
      done
    fi

    modcount=`ls src/*_module.cc 2>/dev/null |wc -l`
    if [ $modcount -gt 0 ]; then
      for file in src/*_module.cc;do
        modname=`echo $file|sed 's|src/||g;s/_module.cc//g'`
        echo "cet_build_plugin($modname art::module
    REG_SOURCE $file
    LIBRARIES REG" >>CMakeLists.txt
        if [ $make_lib -eq 1 ] || [ $make_ilib -eq 1 ]; then
          echo "      offline::$dirname" >>CMakeLists.txt
        fi
        file_lib_list="`grep Offline/ $file|grep -v $dir|sed 's|.*Offline/|offline::|g;s|/.*||g'|uniq|tr '\n' ','`"
        file_lib_list=`echo $file_lib_list|tr ',' '\n'|sort|uniq|awk '{print "      "$1}'`
        echo "$file_lib_list" >>CMakeLists.txt
        echo ")" >>CMakeLists.txt
        echo >> CMakeLists.txt
      done
    fi
    
    toolcount=`ls src/*_tool.cc 2>/dev/null |wc -l`
    if [ $toolcount -gt 0 ]; then
      for file in src/*_tool.cc;do
        toolname=`echo $file|sed 's|src/||g;s/_tool.cc//g'`
        echo "cet_build_plugin($toolname art::tool
    REG_SOURCE $file
    LIBRARIES REG" >>CMakeLists.txt
        if [ $make_lib -eq 1 ] || [ $make_ilib -eq 1 ]; then
          echo "      offline::$dirname" >>CMakeLists.txt
        fi
        file_lib_list="`grep Offline/ $file|grep -v $dir|sed 's|.*Offline/|offline::|g;s|/.*||g'|uniq|tr '\n' ','`"
        file_lib_list=`echo $file_lib_list|tr ',' '\n'|sort|uniq|awk '{print "      "$1}'`
        echo "$file_lib_list" >>CMakeLists.txt
        echo ")" >>CMakeLists.txt
        echo >> CMakeLists.txt
      done
    fi

    if [ -f src/classes.h ]; then
      echo 'art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults'  >>CMakeLists.txt
      echo '    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml' >>CMakeLists.txt
      echo '    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h' >>CMakeLists.txt
      echo '     DICTIONARY_LIBRARIES' >>CMakeLists.txt
        if [ $make_lib -eq 1 ] || [ $make_ilib -eq 1 ]; then
          echo "      offline::$dirname" >>CMakeLists.txt
        fi
      echo ')' >>CMakeLists.txt
    fi
      
    if [ -d src ]; then
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

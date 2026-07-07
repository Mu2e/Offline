//
// swig interface file to wrap c++ code for python
//

%module DataProducts

%include "std_string.i"
%include "math.i"
%include "stdint.i"

%ignore *::operator<<;

%{
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
%}

%include "Offline/DataProducts/inc/StrawEnd.hh"
%include "Offline/DataProducts/inc/StrawId.hh"


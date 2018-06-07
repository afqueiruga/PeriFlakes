
%module bonds_lib
%{
#define SWIG_FILE_WITH_INIT
#include "line_intersection.h"
#include "bond_pressure.h"
%}

%include "kernel.h"
//%include "numpy.i"
//%init %{
//  import_array();
//%}

%include "line_intersection.h"
%include "bond_pressure.h"
/* Wrappers to call them */
%inline %{

%}

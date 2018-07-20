
%module ficticious_lib
%{
#define SWIG_FILE_WITH_INIT
#include "bobaru_n.h"
#include "bobaru_y.h"
#include "bobaru_n3.h"
%}

%include "kernel.h"
//%include "numpy.i"
//%init %{
//  import_array();
//%}

%include "bobaru_n.h"
%include "bobaru_y.h"
%include "bobaru_n3.h"
/* Wrappers to call them */
%inline %{

%}

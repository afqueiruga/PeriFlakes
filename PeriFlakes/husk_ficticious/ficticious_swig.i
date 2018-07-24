
%module ficticious_lib
%{
#define SWIG_FILE_WITH_INIT
#include "bobaru_n.h"
#include "bobaru_F3.h"
#include "bobaru_F.h"
#include "bobaru_n3.h"
#include "bobaru_y.h"
%}

%include "kernel.h"
//%include "numpy.i"
//%init %{
//  import_array();
//%}

%include "bobaru_n.h"
%include "bobaru_F3.h"
%include "bobaru_F.h"
%include "bobaru_n3.h"
%include "bobaru_y.h"
/* Wrappers to call them */
%inline %{

%}

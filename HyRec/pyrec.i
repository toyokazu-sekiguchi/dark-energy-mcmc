%module pyrec
%{
#include "pyrec.h"
%}

%typemap(in) double[3](double temp[3]) {
  if (PyTuple_Check($input)) {
    if (!PyArg_ParseTuple($input,"ddd",temp,temp+1,temp+2)) {
      PyErr_SetString(PyExc_TypeError,"tuple must have 3 elements");
      SWIG_fail;
    }
    $1 = &temp[0];
  } else {
    PyErr_SetString(PyExc_TypeError,"expected a tuple.");
    SWIG_fail;
  }
 }


extern void rec_build_history_wrap(double tcmb, double obh2, double odmh2, double okh2, double odeh2, 
				   double w0, double wa, double yp, double nnu, double mnu[3]);
extern double hyrec_xe(double a);
extern double hyrec_tm(double a);

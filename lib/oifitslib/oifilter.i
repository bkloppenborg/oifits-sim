// $Id: oifilter.i,v 1.2 2007-11-01 17:00:32 jsy1001 Exp $

// oifilter python module - SWIG interface definition file
%define DOCSTRING
"This module provides an object-oriented interface to the OIFITSlib
dataset filtering functionality. To filter a oifits.OiFits instance,
create an instance of the Filter class, make the appropriate
assignments to its data attributes, than pass the OiFits instance to
its apply() method (which returns a filtered version of the dataset).

>>> import oifits, oifilter
>>> o = oifits.OiFits('testdata.fits')
>>> f = oifilter.Filter()
>>> f.accept_t3amp = 0
>>> f.accept_t3phi = 0
>>> filtered = f.apply(o)
>>> print len(filtered.t3List)
0

Examples/tests of other Filter attributes:

>>> f = oifilter.Filter()
>>> f.arrname = 'COAST'
>>> print f.apply(o).arrayList[0].arrname
COAST
>>> f = oifilter.Filter()
>>> f.insname = 'COAST_NICMOS'
>>> print f.apply(o).wavelengthList[0].insname
COAST_NICMOS
>>> f = oifilter.Filter()
>>> f.target_id = 2
>>> print f.apply(o).numTargets
0
>>> f = oifilter.Filter()
>>> f.mjd_range = (51836.0, 51836.9590)
>>> print f.apply(o).vis2List[0].numrec
1
>>> f = oifilter.Filter()
>>> f.wave_range = (2.0e-6, 2.4e-6)
>>> filtered = f.apply(o)
>>> print len(filtered.visList) + len(filtered.vis2List) + len(filtered.t3List)
0
>>> f = oifilter.Filter()
>>> f.accept_vis = 0
>>> print len(f.apply(o).visList)
0
>>> f = oifilter.Filter()
>>> f.accept_vis2 = 0
>>> print len(f.apply(o).vis2List)
0
"
%enddef


%include "oifits_typemaps.i"

%module(docstring=DOCSTRING) oifilter
%pythoncode
%{
def _test():
  import doctest
  doctest.testmod()

if __name__ == '__main__':
  _test()
%}


%{
#include "oifilter.h"
%}

// Apply typemaps from oifits_typemaps.i
%apply char NULL_TERMINATED [ANY] {char [ANY]};
%apply double TUPLE_INPUT [ANY]  {double mjd_range [2]};
%apply double TUPLE_OUTPUT [ANY] {double mjd_range [2]};
%apply float TUPLE_INPUT [ANY]  {float wave_range [2]};
%apply float TUPLE_OUTPUT [ANY] {float wave_range [2]};


%rename(Filter) oi_filter_spec;
typedef struct {
  char arrname[FLEN_VALUE];
  char insname[FLEN_VALUE];
  int target_id;
  double mjd_range[2];
  float wave_range[2];
  int accept_vis;
  int accept_vis2;
  int accept_t3amp;
  int accept_t3phi;
} oi_filter_spec;

// Object-oriented interface to filter
%extend oi_filter_spec
{
  oi_filter_spec()
  {
    oi_filter_spec *self;
    self = (oi_filter_spec *) malloc(sizeof(oi_filter_spec));
    init_oi_filter(self);
    return self;
  }

  ~oi_filter_spec()
  {
    free(self);
  }
  
  const char *__str__()
  {
    return format_oi_filter(self);
  }
  
  oi_fits *apply(const oi_fits *inOi)
  {
    oi_fits *outOi = (oi_fits *) malloc(sizeof(oi_fits));
    apply_oi_filter(inOi, self, outOi);
    return outOi;
  }
  
}

// Local Variables:
// mode: C
// End:

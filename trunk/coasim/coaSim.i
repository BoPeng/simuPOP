%module coaSim

//
//interface file for the python/swig binding of coasim
//
//
/////////////// VECTORS AND MATRICES FOR PARAMETER INPUT/////////////

%include "std_vector.i"
%include "std_string.i"
%include "stl.i"

namespace std
{
  %template()     vector< unsigned int >;
}



///////////////////////// EXCEPTION HANDLING ////////////////////////

%include exception.i

// exceptions will be passed to and correctly handled by python
%exception
{
  try
  {
    $function
  }
  catch(core::illegal_position e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(core::empty_interval e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(core::interval_out_of_range e)
  {
    SWIG_exception(SWIG_IndexError, e.what());
  }
  catch(core::non_pos_pop_size e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(core::out_of_sequence e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(...)
  {
    SWIG_exception(SWIG_RuntimeError, "Unknown runtime error happened.");
  }
}


////////////////////////// INCLUDE FILES ////////////////////////////

%{
#include "all_markers.hh"
#include "builder.hh"
#include "builder_events.hh"
#include "compile_options.hh"
#include "configuration.hh"
#include "descender.hh"
#include "dist_funcs.hh"
#include "interval.hh"
#include "marker.hh"
#include "epochs.hh"
#include "micro_satellite_marker.hh"
#include "monitor.hh"
#include "node.hh"
#include "retired_interval.hh"
#include "simulator.hh"
#include "snp_marker.hh"
#include "trait_marker.hh" 
%}

namespace std
{
  %template(MarkerVec) vector< core::Marker* >;
  %template(EpochVec)  vector< core::Event* >;
}

////////////////////////// SWIG_INIT FUNCTION ///////////////////////
%inline
%{
  // required, but we have nothing to initialize yet.
  bool initialize()
  {
    return true;
  }  
%}

%init
%{
  initialize();
%}

////////////////////////// WRAP THESE CLASSES ///////////////////////

%include "all_markers.hh"
%include "builder.hh"
%include "marker.hh"
%include "configuration.hh"
%include "builder_events.hh"
%include "compile_options.hh"
%include "descender.hh"
%include "dist_funcs.hh"
%include "interval.hh"
%include "epochs.hh"
%include "micro_satellite_marker.hh"
%include "monitor.hh"
%include "node.hh"
%include "retired_interval.hh"
%include "simulator.hh"
%include "snp_marker.hh"
%include "trait_marker.hh"

//////////////////////// PYTHON UTILITY FUNCTION /////////////////////

%pythoncode %{

def blah():
  pass
  
%}

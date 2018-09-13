#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#define R_NO_REMAP

#include <Rinternals.h>

#include <covafill/Core>
#include <covafill/Tree>
#include "utils/convert.hpp"

#include "headers/fill.hpp"
#include "headers/tree.hpp"



extern "C" {
#include <R_ext/Rdynload.h>
#define CALLDEF(name,n) {#name, (DL_FUNC) &name, n}
  static const
  R_CallMethodDef callMethods[] = {
    CALLDEF(MakeFill,4),
    CALLDEF(getFillDim,1),
    CALLDEF(getFillDegree,1),
    CALLDEF(getFillBandwith,1),
    CALLDEF(setFillBandwith,2),
    CALLDEF(getFillObservations,1),
    CALLDEF(getFillCoordinates,1),
    CALLDEF(predictFill,2),
    CALLDEF(predictFillSE,2),
    CALLDEF(lnoResiduals,2),
    CALLDEF(MakeTree,5),
    CALLDEF(getTreeDim,1),
    CALLDEF(predictTree,2),
    {NULL,NULL,0}
  };

  void R_init_covafillr(DllInfo *info)
  {
    /* Register the .C and .Call routines.
       No .Fortran() or .External() routines,
       so pass those arrays as NULL.
    */
    R_registerRoutines(info,
  		       NULL, callMethods,
  		       NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info,TRUE);
  }

  void R_unload_covafillr(DllInfo *info)
  {

  }
}


#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ReadBED(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"ReadBED",            (DL_FUNC) &ReadBED,             6},
  {NULL, NULL, 0}
};

void R_init_SCAT(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

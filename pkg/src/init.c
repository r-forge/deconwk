/* Copyright (C) 2010 */
/* Berwin A Turlach <Berwin.Turlach@gmail.com> */

/* This program is free software; you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation; either version 2 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along */
/* with this program; if not, write to the Free Software Foundation, Inc., */
/* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */
  
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "DeconWK.h"

static
R_NativePrimitiveArgType getmin_t[] =
  {REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP};

static const
R_FortranMethodDef FortEntries[] = {
  {"getmin", (DL_FUNC) &F77_SUB(getmin), 7, getmin_t},
  {NULL, NULL, 0}
};

void attribute_visible
R_init_DeconWK(DllInfo *dll){
R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
R_useDynamicSymbols(dll, FALSE);
}

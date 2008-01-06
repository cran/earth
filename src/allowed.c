// allowed.c: routines for the "allowed" parameter of the R function earth().

#include "R.h"
#include "Rinternals.h"
#define printf Rprintf

#define INLINE inline

#ifndef bool
    typedef int bool;
    #define false 0
    #define true  1
#endif

#define Dirs_(iTerm,iPred)      Dirs[(iTerm) + (iPred)*(nMaxTerms)]

static SEXP AllowedFunc;
static SEXP AllowedEnv;

// Initialize the R function AllowedFunc from the Allowed functon
// argument which was passed into ForwardPassR.
// For efficiency, we initialize once here rather than in IsAllowed.
//
// The caller of ForwardPassR has already checked that Allowed is
// a function and has three arguments: degree, pred, parents.
//
// The "allowed" function has the following prototype:
//
//     allowed <- function(degree, pred, parents)
//     {
//         ...
//         TRUE   # return TRUE if not constrained
//     }
//
// where "degree" is the MARS term degree, with pred in the term.
//       "pred" is column index in the input matrix x
//       "parents" is an integer vector of parent predictors
//                 (it's a copy of Dirs[iParent,]

void InitAllowedFunc(SEXP Allowed, SEXP Env, int nPreds)
{
    if (isNull(Allowed))
        AllowedFunc = R_NilValue;
    else {
        AllowedEnv = Env;

        // the UNPROTECT for the PROTECT below is in FreeAllowedFunc()
        PROTECT(AllowedFunc = allocList(4));

        SEXP s = AllowedFunc;   // 1st element is the function
        SETCAR(s, Allowed);
        SET_TYPEOF(s, LANGSXP);

        s = CDR(s);             // 2nd element is "degree"
        SETCAR(s, allocVector(INTSXP, 1));

        s = CDR(s);             // 3rd element is "pred"
        SETCAR(s, allocVector(INTSXP, 1));

        s = CDR(s);             // 4th element is "parents"
        SETCAR(s, allocVector(INTSXP, nPreds));
    }
}

void FreeAllowedFunc(void)
{
     if (!isNull(AllowedFunc))
         UNPROTECT(1);          // matches PROTECT in InitAllowedFunc
}

// this uses the globals AllowedFunc and AllowedEnv

static INLINE bool EvalAllowedFunc(void)
{
    SEXP s = eval(AllowedFunc, AllowedEnv);

    int allowed;
    switch(TYPEOF(s)) {         // be fairly permissive with return type
        case LGLSXP:
            allowed = (int)(LOGICAL(s)[0]);
            break;
        case INTSXP:
            allowed = INTEGER(s)[0];
            break;
        case REALSXP:
            allowed = (int)(REAL(s)[0]);
            break;
        default:
            error("the \"allowed\" function returned a %s instead of a logical",
                  Rf_type2char(TYPEOF(s)));
            allowed = 0; // -Wall
            break;
    }
    if (LENGTH(s) != 1)
        error("the \"allowed\" function did not return a logical of length 1");

    return allowed;
}

// Return TRUE if the current iPred can be used in a term with iParent
// i.e. TRUE means no constraint.
//
// This calls the R function Allowed which was passed in as a parameter to
// ForwardPassR.  The fields of Allowed have been preallocated into
// AllowedFunc and so all we do here is fill in the values and call eval.

bool IsAllowed(
    const int iPred,        // in: candidate predictor
    const int iParent,      // in: candidate parent term
    const int Dirs[],       // in:
    const int nPreds,       // in:
    const int nMaxTerms)    // in:
{
    if (isNull(AllowedFunc))
       return TRUE;

    SEXP s = AllowedFunc;           // 1st element is the function
    s = CDR(s);                     // 2nd element is "degree"
    INTEGER(CADR(s))[0] = iPred+1;  // 3rd element is "pred"
    int *p = INTEGER(CADDR(s));     // 4th element is "parents"
    int i, nDegree = 1;
    for (i = 0; i < nPreds; i++) {
        p[i] = Dirs_(iParent, i);
        if (p[i])
            nDegree++;
    }
    INTEGER(CAR(s))[0] = nDegree;

    return EvalAllowedFunc();
}

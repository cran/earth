// allowed.c: routines for the "allowed" parameter of the R function earth().

#include "R.h"
#include "Rinternals.h"

#ifndef _MSC_VER // microsoft
#ifndef bool
    typedef int bool;
    #define false 0
    #define true  1
#endif
#endif

#define Dirs_(iTerm,iPred) Dirs[(iTerm) + (iPred)*(nMaxTerms)]

static SEXP AllowedFuncGlobal;
static SEXP AllowedEnvGlobal;
static int  nArgsGlobal;
static bool FirstGlobal;

// Initialize the R function AllowedFuncGlobal from the Allowed function
// argument which was passed into ForwardPassR.
// For efficiency, we initialize once here rather than in IsAllowed.
//
// The caller of ForwardPassR has already checked that Allowed is
// a function and has three arguments: degree, pred, parents.
//
// The "allowed" function has the following prototype, where
// namesx and first are optional.
//
//     allowed <- function(degree, pred, parents, namesx, first)
//     {
//         ...
//         TRUE   # return TRUE if allowed
//     }
//
// where "degree" is the MARS term degree, with pred in the term.
//       "pred" is column index in the input matrix x
//       "parents" is an integer vector of parent predictors
//                 (it's a copy of Dirs[iParent,]
//       "namesx" is optional and is the colnames of the x arg
//                to earth, after factor expansion
//       "first" is optional and is 1 the first time "allowed"
//               is invoked for the current model

void InitAllowedFunc(
        SEXP Allowed, // can be NULL
        int nAllowedArgs, SEXP Env,
        const char** sPredNames, int nPreds)
{
    if(Allowed == R_NilValue)
        AllowedFuncGlobal = NULL;
    else {
        if(nAllowedArgs < 3 || nAllowedArgs > 5)
            error("Bad nAllowedArgs %d", nAllowedArgs);

        AllowedEnvGlobal = Env;
        nArgsGlobal = nAllowedArgs;

        // the UNPROTECT for the PROTECT below is in FreeAllowedFunc()
        PROTECT(AllowedFuncGlobal = allocList(1 + nAllowedArgs));

        SEXP s = AllowedFuncGlobal; // 1st element is the function
        SETCAR(s, Allowed);
        SET_TYPEOF(s, LANGSXP);

        s = CDR(s);                 // 2nd element is "degree"
        SETCAR(s, allocVector(INTSXP, 1));

        s = CDR(s);                 // 3rd element is "pred"
        SETCAR(s, allocVector(INTSXP, 1));

        s = CDR(s);                 // 4th element is "parents"
        SETCAR(s, allocVector(INTSXP, nPreds));

        if(nAllowedArgs >= 4) {
            SEXP namesx;
            s = CDR(s);             // 5th element is "namesx"
            SETCAR(s, namesx = allocVector(STRSXP, nPreds));
            PROTECT(namesx);
            if(sPredNames == NULL)
                error("Bad sPredNames");
            for(int i = 0; i < nPreds; i++)
                SET_STRING_ELT(namesx, i, mkChar(sPredNames[i]));
            UNPROTECT(1);
        }
        if(nAllowedArgs >= 5) {
            s = CDR(s);             // 6th element is "first"
            SETCAR(s, allocVector(LGLSXP, 1));
        }
    }
    FirstGlobal = true;
}

void FreeAllowedFunc(void)
{
    if(AllowedFuncGlobal != NULL) {
        UNPROTECT(1);           // matches PROTECT in InitAllowedFunc
        AllowedFuncGlobal = NULL;
    }
}

static bool EvalAllowedFunc(void)
{
    if(AllowedFuncGlobal == NULL)
        error("EvalAllowedFunc: AllowedFuncGlobal == NULL");

    SEXP s = eval(AllowedFuncGlobal, AllowedEnvGlobal);

    bool allowed;
    switch(TYPEOF(s)) {         // be fairly permissive with return type
        case LGLSXP:
            allowed = (bool)(LOGICAL(s)[0] != 0);
            break;
        case INTSXP:
            allowed = INTEGER(s)[0] != 0;
            break;
        case REALSXP:
            allowed = (bool)(REAL(s)[0] != 0.);
            break;
        default:
            error("the \"allowed\" function returned a %s instead of a logical",
                  Rf_type2char(TYPEOF(s)));
            allowed = FALSE; // -Wall
            break;
    }
    if(LENGTH(s) != 1)
        error("the \"allowed\" function did not return a logical of length 1");

    return allowed;
}

// Return TRUE if the current iPred can be used in a term with iParent
// i.e. TRUE means no constraint.
//
// This calls the R function Allowed which was passed in as a parameter to
// ForwardPassR.  The fields of Allowed have been preallocated into
// AllowedFuncGlobal and so all we do here is fill in the values and call eval.

bool IsAllowed(
    const int iPred,        // in: candidate predictor
    const int iParent,      // in: candidate parent term
    const int Dirs[],       // in:
    const int nPreds,       // in:
    const int nMaxTerms)    // in:
{
    if(AllowedFuncGlobal == NULL)
       return TRUE;

    SEXP s = AllowedFuncGlobal;     // 1st element is the function
    s = CDR(s);                     // 2nd element is "degree"
    INTEGER(CADR(s))[0] = iPred+1;  // 3rd element is "pred"
    int* p = INTEGER(CADDR(s));     // 4th element is "parents"
    int i, nDegree = 1;
    for(i = 0; i < nPreds; i++) {
        p[i] = Dirs_(iParent, i);
        if(p[i])
            nDegree++;
    }
    INTEGER(CAR(s))[0] = nDegree;

    // optional 5th element already initialized to predictor names

    if(nArgsGlobal >= 5)            // optional 6th element is "first"
        *(LOGICAL(CAD4R(s))) = FirstGlobal;
    FirstGlobal = false;

    return EvalAllowedFunc();
}

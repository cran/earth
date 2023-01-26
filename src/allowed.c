// allowed.c: routines for the "allowed" parameter of the R function earth().

#include "R.h"
#include "Rinternals.h"
#include <stdbool.h> // defines bool, true, false

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
// a function and has 3...5 args: degree, pred, parents, namesx, first.
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
//
// This code is based on the "Writing R Extensions" manual
// Section 5.11 "Evaluating R expressions from C"

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

        // We use R_PreserveObject/R_ReleaseObject here instead of
        // PROTECT/UNPROTECT purely to avoid a false warning from CRAN rchk.
        // In the normal course of operation, FreeAllowedFunc() calls
        // R_ReleaseObject to undo the call below to R_PreserveObject.
        // But if there is a call to error in the earth C code (or
        // R_CheckUserInterrupt doesn't return), an on.exit in
        // the R code will call FreeEarth to callFreeAllowedFunc.
        AllowedFuncGlobal = allocList(1 + nAllowedArgs);
        R_PreserveObject(AllowedFuncGlobal);

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
            if(sPredNames == NULL)
                error("Bad sPredNames");
            PROTECT(namesx);
            for(int i = 0; i < nPreds; i++)
                SET_STRING_ELT(namesx, i, mkChar(sPredNames[i]));
            UNPROTECT(1);
        }
        if(nAllowedArgs >= 5) {
            s = CDR(s);             // 6th element is "first"
            SETCAR(s, allocVector(LGLSXP, 1));
        }
    }
    FirstGlobal = TRUE;
}

void FreeAllowedFunc(void)
{
    if(AllowedFuncGlobal != NULL) {
        // following matches R_PreserveObject in InitAllowedFunc
        R_ReleaseObject(AllowedFuncGlobal);
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

    // AllowedFuncGlobal has been protected by R_PreserveObject
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
        *(LOGICAL(CAD4R(s))) = FirstGlobal != 0;
    FirstGlobal = FALSE;

    return EvalAllowedFunc();
}

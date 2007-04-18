// earth.c
//
// This code is derived from code in mda:dmarss.r by Hastie and Tibshirani.
// Comments containing "$$" mark known issues
// Stephen Milborrow Feb 2007 Petaluma
//
// See the R earth documentation for descriptions of the principal data structures.
//
// This code uses a subset of C99.  To build the earth R library under gcc:
//   gcc -Wall -pedantic -Wextra -O3 -std=gnu99 -Ic:/a1/r/work/include -c earth.c -o earth.o
//
// For a standalone program with a main() routine under gcc (you will have to
// build Rddl.lib and Rblas.lib first):
//   gcc -DSTANDALONE -DMAIN -Wall -pedantic -Wextra -O3 -std=gnu99
//      -Ic:/a1/r/work/include -I../src/tests ../src/earth.c 
//      /a1/r/work/src/gnuwin32/Rdll.lib /a1/r/work/bin/Rblas.lib -o main.exe
//   main.exe
//
// For a standalone program with a main() routine under VisualC 6.0,
// treat earth.c as a C++ file because VisualC has no direct support for C99:
// I have not tried building under higher versions of VisualC.
//   md Debug
//   cl -nologo -DSTANDALONE -DMAIN -D_CRTDBG_MAP_ALLOC -TP -Zi -W3 -MLd
//       -I/a1/r/work/src/include -I../src/tests -Fp./Debug/vc60.PCH -Fo"./Debug/" 
//       -c ../src/earth.c
//   link -nologo -debug:full -out:main.exe ./Debug/earth.obj
//       \a1\r\work\src\gnuwin32\Rdll.lib \a1\r\work\bin\Rblas.lib
//   main.exe
//
// References:
//
// HastieTibs: Trevor Hastie and Robert Tibshirani
//      S library mda version 0.3.2 dmarss.r Ratfor code
//      Modifications for R by Kurt Hornik, Friedrich Leisch, Brian Ripley
//
// FriedmanMars: Multivariate Adaptive Regression Splines (with discussion)
//      Annals of Statistics 19/1, 1--141, 1991
//
// FriedmanFastMars: Friedman "Fast MARS"
//      Dep. of Stats. Stanford, Tech Report 110, May 1993

#if !STANDALONE
#define USING_R 1
#endif // STANDALONE

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#if _DEBUG && _MSC_VER
    #include <crtdbg.h> // microsoft malloc debugging library
#endif

#if _MSC_VER            // microsoft
    #define _C_ "C"
#else
    #define _C_
    #ifndef bool
        typedef int bool;
        #define false 0
        #define true  1
    #endif
#endif

#if USING_R             // R with gcc
    #include "R.h"
    #define printf Rprintf
    #define ASSERT(x)   \
        if (!(x)) error("assertion failed: %s [%s %d]", #x, __FILE__, __LINE__)
#else
    #include "earth.h"
    #define warning printf
    void error(const char *args, ...);
    #define R_FINITE(x) _finite(x)
    #define ASSERT(x)   \
        if (!(x)) error("assertion failed: %s [%s %d]\n", #x, __FILE__, __LINE__)
#endif

extern _C_ int dqrdc2_(double *x, int *ldx, int *n, int *p,
              double *tol, int *rank,
              double *qraux, int *pivot, double *work);

extern _C_ int dqrsl_(double *x, int *ldx, int *n, int *k,
             double *qraux, double *y,
             double *qy, double *qty, double *b,
             double *rsd, double *xb, int *job, int *info);

extern _C_ void dtrsl_(double *t, int *ldt, int *n, double *b, int *job, int *info);

extern _C_ void daxpy_(const int *n, const double *alpha,
                    const double *dx, const int *incx,
                    double *dy, const int *incy);

extern _C_ double ddot_(const int *n,
                    const double *dx, const int *incx,
                    const double *dy, const int *incy);

#define sq(x)       ((x) * (x))
#ifndef max
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif

#define INLINE      inline
#define USE_BLAS    1     // 1 is faster (tested on Windows XP Pentium with R BLAS)
                          // also, need USE_BLAS to use bxOrthCenteredT
#define BETA_CACHE  1     // 1 to cache betas in orthogResiduals, for speed
#define FAST_MARS   1     // 1 to use technique in FriedmanFastMars
#define IOFFSET     1     // 1 to convert 0-based indices to 1-based in printfs

static const char   *VERSION    = "version 0.1"; // change if you modify this file!
static const double BX_TOL      = 0.01;
static const double QR_TOL      = 0.01;
static const double MIN_GRSQ    = -10.0;
static const int    ONE         = 1;        // parameter for BLAS routines
#if _MSC_VER                                // microsoft compiler
static const double ZERO        = 0;
static const double POS_INF     = (1.0 / ZERO);
#else
static const double POS_INF     = (1.0 / 0.0);
#endif

// Poor man's array indexing -- not pretty, but pretty useful.
//
// Note that we use column major ordering. C programs usually use row major
// ordering but we don't here because the functions in this file are called
// by R and call Fortran routines which use column major ordering.

#define Dirs_(iTerm,iPred)      Dirs[(iTerm) + (iPred)*(nMaxTerms)]
#define Cuts_(iTerm,iPred)      Cuts[(iTerm) + (iPred)*(nMaxTerms)]

#define bx_(iCase,iTerm)                bx      [(iCase) + (iTerm)*(nCases)]
#define bxOrth_(iCase,iTerm)            bxOrth  [(iCase) + (iTerm)*(nCases)]
#define bxOrthCenteredT_(iTerm,iCase)   bxOrthCenteredT[(iTerm) + (iCase)*(nMaxTerms)]
#define x_(iCase,iPred)                 x       [(iCase) + (iPred)*(nCases)]
#define xOrder_(iCase,iPred)            xOrder  [(iCase) + (iPred)*(nCases)]

static int nTraceGlobal;        // copy of nTrace parameter
static int nMinSpanGlobal;      // copy of nMinSpan parameter

#if BETA_CACHE
static void FreeBetaCache(void);
#endif

//--------------------------------------------------------------------------------------------
#define free1(p) if (p) free(p); p = NULL   // is a macro so we can zero p

static void *malloc1(size_t size)
{
    void *p = malloc(size);
    if (!p)
        error("Out of memory (could not allocate %g MBytes)", ((double)size) / sq(1024.0));
    return p;
}

static void *calloc1(size_t num, size_t size)
{
    void *p = calloc(num, size);
    if (!p)
        error("Out of memory (could not allocate %g MBytes)", ((double)size) / sq(1024.0));
    return p;
}

//--------------------------------------------------------------------------------------------
// These are malloced blocks.  They unfortunately have to be declared globally so
// under R if the user interrupts we can free them using on.exit(.C("FreeR"))

static int  *xOrder;            // local to FindTerm
static bool *WorkingSet;        // local to FindTerm and backwardPass
static double *xbx;             // local to FindTerm
static double *CovSx;           // local to FindTerm
static double *CovCol;          // local to FindTerm
static double *CovSy;           // local to FindTerm
static double *bxOrth;          // local to ForwardPass

// Transposed and mean centered copy of bxOrth, for fast update in FindKnot.
// It's faster because there is better data locality as iTerm increases, so
// better L1 cache use.  Using this is measurably faster.
// This is used only if USE_BLAS is true.
static double *bxOrthCenteredT; // local to ForwardPass

static double *bxOrthMean;      // local to ForwardPass
static int  *nFactorsInTerm;    // local to earth or ForwardPassR
static int  *nUses;             // local to earth or ForwardPassR
#if USING_R
static int *iDirs;              // local to ForwardPassR
static bool *BoolFullSet;       // local to ForwardPassR
#endif
#if FAST_MARS
static void FreeQ(void);
#endif

#if USING_R
void FreeR(void)
{
    free1(WorkingSet);
    free1(CovSx);
    free1(CovCol);
    free1(CovSy);
    free1(xOrder);
    free1(bxOrthMean);
    free1(bxOrthCenteredT);
    free1(bxOrth);
    free1(BoolFullSet);
    free1(iDirs);
    free1(nUses);
    free1(nFactorsInTerm);
#if BETA_CACHE
    FreeBetaCache();
#endif
#if FAST_MARS
    FreeQ();
#endif
}
#endif

//--------------------------------------------------------------------------------------------
// Gets called periodically to service the R framework.
// Will not return if the user interrupts.

#if USING_R
static INLINE void ServiceR(int nCases)
{
    static unsigned nIters;
    static int nCasesPrev;
    static int nItersPerRCall;
    if (nCases != nCasesPrev) {
        // Init nItersPerRCall, knowing that FindKnot services nCases per call
        // to ServiceR.  We want to call R_CheckUserInterrupt often enough to be
        // responsive to interrupts but not so often that we slow computation.
        // The magic nbrs below seem ok in practice.

        nCasesPrev = nCases;
        nItersPerRCall = 100000/nCases;
        if (nItersPerRCall > 100)
            nItersPerRCall = 100;
        if (nItersPerRCall < 1)
            nItersPerRCall = 1;
    }
    if (nIters++ % nItersPerRCall == 0) {
        R_FlushConsole();
        R_CheckUserInterrupt();     // may never return
    }
}
#endif

//--------------------------------------------------------------------------------------------
// Order() gets the sorted indices of vector x, so x[sorted[i]] <= x[sorted[i+1]].
// Ties may be reordered. The returned indices are 0 based (as in C not as in R).
//
// This function is similar to the R library function rsort_with_index(),
// but is defined here to minimize R dependencies.

static const double *pxGlobal;

static int Compare(const void *p1, const void *p2)  // for qsort
{
    const int i1 = *(int *)p1;
    const int i2 = *(int *)p2;
    double Diff = pxGlobal[i1] - pxGlobal[i2];
    if (Diff < 0)
        return -1;
    else if (Diff > 0)
        return 1;
    else
        return 0;
}

static void Order(int sorted[],                     // out: vector with nx elements
                  const double x[], const int nx)   // in: x is a vector with nx elements
{
    for (int i = 0; i < nx; i++)
        sorted[i] = i;
    pxGlobal = x;
    qsort(sorted, nx, sizeof(int), Compare);
}

//--------------------------------------------------------------------------------------------
// Get order indices for an x array of dimensions nRows x nCols.
//
// Returns an nRows x nCols integer array of indices, where each column
// corresponds to a column of x.  See Order() for ordering details.
//
// Caller must free the returned array.

static int *OrderArray(const double x[], const int nRows, const int nCols)
{
    int *xOrder = (int *)malloc1(nRows * nCols * sizeof(int));

    for (int iCol = 0; iCol < nCols; iCol++)
        Order(xOrder + iCol*nRows, x + iCol*nRows, nRows);

    return xOrder;
}

//--------------------------------------------------------------------------------------------
// return the number of TRUEs in the boolean vector UsedCols

static int GetNbrUsedCols(const bool UsedCols[], const int nLen)
{
    ASSERT(UsedCols[0] == true);    // intercept term should always be set

    int nTrue = 0;

    for (int iCol = 0; iCol < nLen; iCol++)
        if (UsedCols[iCol])
            nTrue++;

    return nTrue;
}

//--------------------------------------------------------------------------------------------
// Copy used columns in x to *pxUsed and return the number of used columns
// UsedCols[i] is true for each each used column index in x
// Caller must free *pxUsed

static int CopyUsedCols(double **pxUsed,                // out
                    const double x[],                   // in: nCases x nCols
                    const int nCases, const int nCols,  // in
                    const bool UsedCols[])              // in
{
    int nUsedCols = GetNbrUsedCols(UsedCols, nCols);

    double *xUsed = (double *)malloc1(nCases * nUsedCols * sizeof(double));
    int iUsed = 0;
    for (int iCol = 0; iCol < nCols; iCol++)
        if (UsedCols[iCol]) {
            memcpy(xUsed + iUsed * nCases, x + iCol * nCases, nCases * sizeof(double));
            iUsed++;
        }
    *pxUsed = xUsed;
    return nUsedCols;
}

//--------------------------------------------------------------------------------------------
// Regress y on the used columns of x, in the standard way (using QR).
// UsedCols[i] is true for each each used col i in x; unused cols are ignored
//
// The returned Beta argument is computed from, and is indexed on,
// the compacted x vector, not on the original x.
//
// The returned Pivots should only be used if *pnRank != nUsedCols.
// The entries of Pivots refer to columns in the full x (and are 0 based).
// Entries in Pivots at *pnRank and above specify linearly dependent columns in x.
//
// To maximize compatibility we call the same routines as the R function lm.

static void Regress(
    double       Betas[],       // out
    double       Residuals[],   // out
    double       *pRss,         // out
    int          *pnRank,       // out: nbr of indep cols in x
    int          Pivots[],      // out:
    const double x[],           // in: nCases x nCols
    const double y[],           // in: nCases
    const int    nCases,        // in
    int          nCols,         // in: number of columns in x, some may not be used
    const bool   UsedCols[])    // in: specifies used columns in x
{
    double *xUsed;
    int nUsedCols = CopyUsedCols(&xUsed, x, nCases, nCols, UsedCols);
    int iCol;
    for (iCol = 0; iCol < nCols; iCol++)
        Pivots[iCol] = iCol+1;

    // compute Betas and yHat (use Residuals as a temporary buffer to store yHat)

    int job = 101;          // specify c=1 e=1 to compute qty, b, yHat
    int info;

    double *qraux = (double *)malloc1(nUsedCols * sizeof(double));
    double *work = (double *)malloc1(nCases * nUsedCols * sizeof(double));

    dqrdc2_(                // R function, based on LINPACK dqrdc
        xUsed,              // io:  x, on return upper tri of x is R of QR
        (int *)&nCases,     // in:  ldx (typecast discards const)
        (int *)&nCases,     // in:  n
        &nUsedCols,         // in:  p
        (double*)&QR_TOL,   // in:  tol
        pnRank,             // out: k, num of indep cols of x
        qraux,              // out: qraux
        Pivots,             // out: jpvt
        work);              // work

    dqrsl_(                 // LINPACK function
        xUsed,              // in:  x
        (int *)&nCases,     // in:  ldx (typecast discards const)
        (int *)&nCases,     // in:  n
        pnRank,             // in:  k
        qraux,              // in:  qraux
        (double *)y,        // in:  y
        work,               // out: qy, unused here
        work,               // out: qty, unused here
        Betas,              // out: b
        work,               // out: rsd, unused here
        Residuals,          // out: xb = yHat = ls approx of x*b
        &job,               // in:  job
        &info);             // in:  info

    ASSERT(info == 0);

    // compute Residuals and *pRss

    double Rss = 0;
    for (int iCase = 0; iCase < nCases; iCase++) {
        Residuals[iCase] = y[iCase] - Residuals[iCase];
        Rss += sq(Residuals[iCase]);
    }
    *pRss = Rss;

    if (*pnRank != nUsedCols) {
        // adjust pivots for missing cols in UsedCols and for 1 offset

        int *PivotOffset = (int *)malloc1(nCols * sizeof(int));
        int iCol, nOffset = 0, iOld = 0;
        for (iCol = 0; iCol < nCols; iCol++) {
            if (!UsedCols[iCol])
                nOffset++;
            else {
                PivotOffset[iOld] = nOffset;
                if (++iOld > nUsedCols)
                    break;
            }
        }
        for (iCol = 0; iCol < nUsedCols; iCol++)
            Pivots[iCol] = Pivots[iCol] - 1 + PivotOffset[Pivots[iCol] - 1];
        free1(PivotOffset);
    }
    free1(qraux);
    free1(xUsed);
    free1(work);
}

//--------------------------------------------------------------------------------------------
// Regress y on bx to get Residuals and Betas.
// If bx isn't of full rank, remove dependent cols, update UsedCols,
// and regress again on the bx with removed cols

static void RegressAndFix(
    double Betas[],         // out: nMaxTerms x 1
    double Residuals[],     // out: nCases x 1
    double *pRss,           // out:
    bool   UsedCols[],      // io:  will remove cols if necessary, nMaxTerms x 1
    const  double bx[],     // in:  nCases x nMaxTerms
    const double y[],       // in:  nCases x 1
    const int nCases,       // in
    const int nTerms)       // in
{
    int nRank;
    int *Pivots = (int *)malloc1(nTerms * sizeof(int));
    Regress(Betas, Residuals, pRss, &nRank, Pivots, bx, y, nCases, nTerms, UsedCols);
    int nUsedCols = GetNbrUsedCols(UsedCols, nTerms);
    int nDeficient = nUsedCols - nRank;
    if (nDeficient) {           // rank deficient?
        // Remove linearly dependent columns.
        // The lin dep columns are at index nRank and higher in Pivots.

        for (int iCol = nRank; iCol < nUsedCols; iCol++)
            UsedCols[Pivots[iCol]] = false;

        Regress(Betas, Residuals, pRss, &nRank, Pivots, bx, y, nCases, nTerms, UsedCols);
        nUsedCols = nUsedCols - nDeficient;
        if (nRank != nUsedCols)
            warning("Could not fix rank deficient bx: nUsedCols %d nRank %d",
                nUsedCols,  nRank);
        else if (nTraceGlobal >= 1)
            printf("Fixed rank deficient bx by removing %d term%s, %d term%s remain\n",
                nDeficient, ((nDeficient==1)? "": "s"),
                nUsedCols,  ((nUsedCols==1)? "": "s"));
        }
    free1(Pivots);
}

//--------------------------------------------------------------------------------------------
#if BETA_CACHE

// The BetaCache is used when searching for a new term pair, via FindTerm.
// Most of the calculation for the orthogonal regression betas is repeated
// with the same data, and thus we can save time by caching betas.
//
// iParent     is the term that forms the base for the new term
// iPred       is the predictor for the new term
// iOrthPred   is the column index in the xOrth array
//
// $$ could reduce the size of BetaCache, repeated values?

static double *BetaCacheGlobal; // [iOrthPred,iParent,iPred]
                                // dim nPreds x nMaxTerms x nMaxTerms

static void InitBetaCache(int nMaxTerms, int nPreds)
{
    int nCacheSize =  nMaxTerms * nMaxTerms * nPreds;

    if (nTraceGlobal >= 4)                  // print cache size in megabytes
        printf("BetaCache %.3f Mb\n", (nCacheSize * sizeof(double))/sq(1024.0));

    BetaCacheGlobal = (double *)malloc1(nCacheSize * sizeof(double));

    for (int i = 0; i < nCacheSize; i++)    // mark all cache entries as uninited
        BetaCacheGlobal[i] = POS_INF;
}

static void FreeBetaCache(void)
{
    if (BetaCacheGlobal)
        free1(BetaCacheGlobal);
}
#endif // BETA_CACHE

//--------------------------------------------------------------------------------------------
static INLINE double Mean(const double x[], int n)
{
    double mean = 0;
    for (int i = 0; i < n; i++)
        mean += x[i] / n;
    return mean;
}

//--------------------------------------------------------------------------------------------
// get mean centered sum of squares

static INLINE double SumOfSquares(const double x[], const double mean, int n)
{
    double ss = 0;
    for (int i = 0; i < n; i++)
        ss += sq(x[i] - mean);
    return ss;
}

//--------------------------------------------------------------------------------------------
static INLINE double GetGcv(const int nTerms, // nbr basis terms including intercept
                const int nCases, double Rss, const double Penalty)
{
    double cost;
    if (Penalty == -1)  // special case: terms and knots are free
        cost = 0;
    else
        cost = nTerms + (Penalty * (double)(nTerms-1) / 2);   // nKnots = (nTerms-1)/2

    return Rss / (nCases * sq(1 - cost/nCases));
}

//--------------------------------------------------------------------------------------------
// We only consider knots that are nMinSpan distance apart, to increase resistance
// to noise.  This function determines that distance.
// It implements eqn 43 FriedmanMars, but with an extension for nMinSpan.
// If bx==NULL then instead of counting valid entries in bx we use nCases,
// and ignore the term index iTerm.
//
// nMinSpan: if <0, add to internally calculated min span (i.e. decrease internal min span)
//           if =0, use internally calculated min span
//           if >0, use instead of internally calculated min span

static INLINE int GetMinSpan(int nCases, const double *bx,
                             const int iTerm, const int nMinSpan)
{
    if (nMinSpan > 0)                           // user specified a fixed span?
        return nMinSpan;

    int nUsed = 0;
    if (bx == NULL)
        nUsed = nCases;
    else for (int iCase = 0; iCase < nCases; iCase++)
        if (bx_(iCase,iTerm) > 0)
            nUsed++;

    static const double Term1 = 2.9702;         // -log(-log(0.95)
    static const double Term2 = 1.7329;         // 2.5 * log(2)

    double MinSpan = (Term1 + log((double)nCases * nUsed)) / Term2;

    MinSpan += nMinSpan;                        // possibly decrease
    if (MinSpan < 1)
        MinSpan = 1;

    return (int)MinSpan;
}

//--------------------------------------------------------------------------------------------
// We don't consider knots that are too close to the ends.
// This function determines how close to an end we can get.
// It implements eqn 45 FriedmanMars, re-expressed for efficient computation

static INLINE int GetEndSpan(const int nCases)
{
    static const double Log2  = 0.69315;            // log(2)
    static const double Term1 = 7.32193;            // 3 + log(20)/log(2);

    return (int)(Term1 + log(nCases) / Log2);
}

//--------------------------------------------------------------------------------------------
// Return true if model term type is not already in model
// i.e. if the hockey stick functions in a pre-existing term do not use the same
// predictors (ignoring knot values).
//
// In practice this nearly always returns true.

static bool GetNewFormFlag(const int iPred, const int iTerm,
                        const int Dirs[], const bool UsedCols[],
                        const int nTerms, const int nPreds, const int nMaxTerms)
{
    bool IsNewForm = true;
    for (int iTerm1 = 1; iTerm1 < nTerms; iTerm1++)     // start at 1 to skip intercept
        if (UsedCols[iTerm1]) {
            IsNewForm = false;
            if (Dirs_(iTerm1,iPred) == 0)
                return true;
            for (int iPred1 = 0; iPred1 < nPreds; iPred1++)
                if (iPred1 != iPred && Dirs_(iTerm1,iPred1) != Dirs_(iTerm,iPred1))
                    return true;
        }
    return IsNewForm;
}
//--------------------------------------------------------------------------------------------
// Return the Residuals from regressing y on the used columns of the
// orthogonal matrix bxOrth.  The mean of each column of bxOrth must be 0.
//
// In practice this function is called with the params shown in {braces}.
//
// This function must be fast.

static INLINE void OrthogResiduals(
    double Residuals[],     // out: vector with nCases elems { bxOrth[,nTerms] }
    const double y[],       // in:  vector with nCases elems { bx[,nTerms], xbx }
    const double bxOrth[],  // in:  nCases x nPreds          { bxOrth }
    const int nCases,       // in
    const int nPreds,       // in: number of cols in bxOrth   { nTerms in model }
    const bool UsedCols[],  // in: UsedCols[i] is TRUE if col is used, unused cols are ignored
                            //     Following parameters are only for the beta cache
    const int iParent,      // in: if >= 0, use BetaCacheGlobal {FindTerm iTerm, addTermP -1}
    const int iPred,        // in: predictor index
    const int nMaxTerms)    // in:
{
#if BETA_CACHE
    double *const pCache = BetaCacheGlobal + iParent*nMaxTerms + iPred*sq(nMaxTerms);
#else
    double Dummy = iParent; // prevent compiler warning: unused parameter
    Dummy = iPred;
    Dummy = nMaxTerms;
#endif
    memcpy(Residuals, y, nCases * sizeof(double));

    for (int iOrthPred = 0; iOrthPred < nPreds; iOrthPred++)
        if (UsedCols[iOrthPred]) {
            double Beta;
            const double *pbxOrth = &bxOrth_(0, iOrthPred);
#if BETA_CACHE
            if (iParent >= 0 && pCache[iOrthPred] != POS_INF)
                Beta = pCache[iOrthPred];
            else
#endif
            {
                double Xtx = 0, Xty = 0;
                for (int iCase = 0; iCase < nCases; iCase++) {
                    const double x1 = pbxOrth[iCase];
                    Xty += x1 * y[iCase];
                    Xtx += sq(x1);
                }
                Beta = Xty / Xtx;
                ASSERT(R_FINITE(Beta));
#if BETA_CACHE
                if (iParent >= 0)
                    pCache[iOrthPred] = Beta;
#endif
            }
#if USE_BLAS
            const double NegBeta = -Beta;
            daxpy_(&nCases, &NegBeta, pbxOrth, &ONE, Residuals, &ONE);
#else
            for (int iCase = 0; iCase < nCases; iCase++)
                Residuals[iCase] -= Beta * pbxOrth[iCase];
#endif
        }
}

//--------------------------------------------------------------------------------------------
// Init the rightmost column of bxOrth i.e the column indexed by nTerms.
// The new col is the normalized residuals from regressing y on the
// lower (i.e. already existing) cols of bxOrth.
// Also updates bxOrthCenteredT and bxOrthMean.
//
// In practice this function is called only with the params shown in {braces}

static INLINE void InitBxOrthCol(
    double bxOrth[],         // io: col nTerms is changed, other cols not touched
    double bxOrthCenteredT[],// io: kept in sync with bxOrth
    double bxOrthMean[],     // io: element at nTerms is updated
    bool   *pGoodCol,        // io: set to false if col sum-of-squares is under BX_TOL
    const double *y,         // in: { FindTerm y, addTermPair bx[,nTerms] }
    const int nTerms,        // in: column goes in at index nTerms, 0 is the intercept
    const bool WorkingSet[], // in
    const int nCases,        // in
    const int nMaxTerms,     // in
    const int iCacheTerm,    // in: if >= 0, use BetaCacheGlobal {FindTerm iTerm, addTermP -1}
                             //     if < 0 then recalc Betas from scratch
    const int iPred)         // in: predictor index, only used if iCacheTerm >= 0
{
    int iCase;
    *pGoodCol = true;

    if (nTerms == 0) {          // column 0, the intercept
        double root = 1 / sqrt(nCases);
        bxOrthMean[0] = root;
        for (iCase = 0; iCase < nCases; iCase++)
            bxOrth_(iCase,0) = root;
    } else if (nTerms == 1) {   // column 1, the first basis function

        double yMean = Mean(y, nCases);
        for (iCase = 0; iCase < nCases; iCase++)
            bxOrth_(iCase,1) = y[iCase] - yMean;
    } else
        OrthogResiduals(&bxOrth_(0,nTerms), // resids go in rightmost col of bxOrth at nTerms
            y, bxOrth, nCases, nTerms, WorkingSet, iCacheTerm, iPred, nMaxTerms);

    if (nTerms > 0) {
        // normalize the column to length 1 and init bxOrthMean[nTerms]

        bxOrthMean[nTerms] = Mean(&bxOrth_(0,nTerms), nCases);
        double bxOrthSS = SumOfSquares(&bxOrth_(0,nTerms), 0, nCases);
        const double Tol = (iCacheTerm < 0? 0: BX_TOL);
        if (bxOrthSS > Tol) {
            const double root = sqrt(bxOrthSS);
            for (iCase = 0; iCase < nCases; iCase++)
                bxOrth_(iCase,nTerms) /= root;
        } else {
            memset(&bxOrth_(0,nTerms), 0, nCases * sizeof(double));
            bxOrthMean[nTerms] = 0;
            *pGoodCol = false;
        }
    }
    for (iCase = 0; iCase < nCases; iCase++)        // keep bxOrthCenteredT in sync
        bxOrthCenteredT_(nTerms,iCase) = bxOrth_(iCase,nTerms) - bxOrthMean[nTerms];
}

//--------------------------------------------------------------------------------------------
static double GetCut(int iCase, const int iPred, const int nCases,
                        const double x[], const int xOrder[])
{
    if (iCase < 0 || iCase >= nCases)
        error("GetCut iCase %d: iCase < 0 || iCase >= nCases", iCase);
    const int ix = xOrder_(iCase,iPred);
    ASSERT(ix >= 0 && ix < nCases);
    return x_(ix,iPred);
}

//--------------------------------------------------------------------------------------------
// Add a new term pair to the arrays.
// Each term in the new term pair is a copy of an existing parent term but extended
// by multiplying it by a new hockey stick function at the selected knot.
// If the upper term in the term pair is invalid then we still add the upper
// term but mark it as FALSE in FullSet.

static void AddTermPair(
    int    Dirs[],              // io
    double Cuts[],              // io
    double bx[],                // io
    double bxOrth[],            // io
    double bxOrthCenteredT[],   // io
    double bxOrthMean[],        // io
    bool   FullSet[],           // io
    int    nFactorsInTerm[],    // io
    int    nUses[],             // io
    const int nTerms,           // in: new term pair goes in at index nTerms and nTerms1
    const int iBestParent,      // in: parent term
    const int iBestCase,        // in
    const int iBestPred,        // in
    const int nPreds,           // in
    const int nCases,           // in
    const int nMaxTerms,        // in
    const bool IsNewForm,       // in
    const bool IsLinPred,       // in
    const double x[],           // in
    const int xOrder[])         // in
{
    const double BestCut = GetCut(iBestCase, iBestPred, nCases, x, xOrder);
    ASSERT(IsLinPred || iBestCase != 0);
    const int nTerms1 = nTerms+1;

    // copy the parent term to the new term pair

    int iPred;
    for (iPred = 0; iPred < nPreds; iPred++) {
        Dirs_(nTerms, iPred) =
        Dirs_(nTerms1,iPred) = Dirs_(iBestParent,iPred);

        Cuts_(nTerms, iPred) =
        Cuts_(nTerms1,iPred) = Cuts_(iBestParent,iPred);

        if (nTraceGlobal >= 2)
            if (Dirs_(iBestParent,iPred) != 0)  // print parent term
                printf("%-3d ", iBestParent+IOFFSET);
    }
    // incoporate the new hockey stick function

    nFactorsInTerm[nTerms]  =
    nFactorsInTerm[nTerms1] = nFactorsInTerm[iBestParent] + 1;

    Dirs_(nTerms, iBestPred) = 1;
    Dirs_(nTerms1,iBestPred) = -1;

    Cuts_(nTerms, iBestPred) =
    Cuts_(nTerms1,iBestPred) = BestCut;

    FullSet[nTerms] = true;
    if (!IsLinPred && IsNewForm)
        FullSet[nTerms1] = true;

    // If the term is not valid, then we don't wan't to use it as the base for
    // a new term later (in FindTerm).  Enforce this by setting
    // nFactorsInTerm to a value greater than any posssible nMaxDegree.

    if (!FullSet[nTerms1])
        nFactorsInTerm[nTerms1] = 99;

    // fill in new columns of bx, at nTerms and nTerms+1 (left and right hockey sticks)

    int iCase;
    for (iCase = 0; iCase < nCases; iCase++)
        if (x_(iCase,iBestPred) - BestCut > 0)
            bx_(iCase,nTerms) = bx_(iCase,iBestParent) * (x_(iCase,iBestPred) - BestCut);
        else
            bx_(iCase,nTerms1) = bx_(iCase,iBestParent) * (BestCut - x_(iCase,iBestPred));

    nUses[iBestPred]++;

    // init the col in bxOrth at nTerms and init bxOrthMean[nTerms]

    bool GoodCol;
    InitBxOrthCol(bxOrth, bxOrthCenteredT, bxOrthMean, &GoodCol,
        &bx_(0,nTerms), nTerms, FullSet, nCases, nMaxTerms, -1, iPred);
                            // -1 means don't use BetaCacheGlobal, calc Betas afresh

    // init the col in bxOrth at nTerms1 and init bxOrthMean[nTerms1]

    if (FullSet[nTerms1]) {
        InitBxOrthCol(bxOrth, bxOrthCenteredT, bxOrthMean, &GoodCol,
            &bx_(0,nTerms1), nTerms1, FullSet, nCases, nMaxTerms, -1, iPred);
    } else {
        memset(&bxOrth_(0,nTerms1), 0, nCases * sizeof(double));
        bxOrthMean[nTerms1] = 0;
        for (iCase = 0; iCase < nCases; iCase++)    // keep bxOrthCenteredT in sync
            bxOrthCenteredT_(nTerms1,iCase) = 0;
    }
}

//--------------------------------------------------------------------------------------------
#if FAST_MARS

// For simplicty we always update the queue
// but we only use the queue if FastK is small enough.

typedef struct tQueue {
    int     iParent;        // parent term
    double  RssDelta;
    int     nTerms;         // number of terms when RssDelta was calculated
    double  AgedRank;
} tQueue;

static tQueue *Q;           // indexed on iTerm (this Q is used for updates)
static tQueue *SortedQ;     // indexed on iParent rank
static int    nQMax;        // number of elements used in queues

static void InitQ(const int nMaxTerms)
{
    int i;
    nQMax = 0;
    Q       = (tQueue *)malloc1(nMaxTerms * sizeof(tQueue));
    SortedQ = (tQueue *)malloc1(nMaxTerms * sizeof(tQueue));
    for (i = 0; i < nMaxTerms; i++) {
        Q[i].iParent = i;
        Q[i].nTerms = -1;   // not strictly needed, nice for debugging
        Q[i].RssDelta = -1;
        Q[i].AgedRank = -1;
    }
}

static void FreeQ(void)
{
    free1(SortedQ);
    free1(Q);
}

static void PrintSortedQ(int FastK)     // for debugging
{
    printf("\n\nQIndex Parent nTerms AgedRank  RssDelta\n");
    for (int i = 0; i < nQMax; i++) {
        printf("%c  %3d    %3d    %3d    %5.1f  %g\n",
            ((i == FastK-1)? 'K':' '),
            i+IOFFSET,
            SortedQ[i].iParent+IOFFSET,
            SortedQ[i].nTerms+IOFFSET,
            SortedQ[i].AgedRank,
            SortedQ[i].RssDelta);
    }
}

// Sort so highest RssDeltas are at low indices.
// Secondary sort key is iParent.  Not strictly needed, but removes
// possible differences in qsort implementations (which "sorts"
// identical keys unpredictably).

static int CompareQ(const void *p1, const void *p2)     // for qsort
{
    double Diff = ((tQueue*)p2)->RssDelta - ((tQueue*)p1)->RssDelta;
    if (Diff < 0)
        return -1;
    else if (Diff > 0)
        return 1;

    // Diff is 0, so sort now on iParent

    int iDiff = ((tQueue*)p1)->iParent - ((tQueue*)p2)->iParent;
    if (iDiff < 0)
        return -1;
    else if (iDiff > 0)
        return 1;
    return 0;
}

// Sort so lowest AgedRanks are at low indices.
// If AgedRanks are the same then sort on RssDelta and iParent.

static int CompareAgedQ(const void *p1, const void *p2) // for qsort
{
    double Diff = ((tQueue*)p1)->AgedRank - ((tQueue*)p2)->AgedRank;
    if (Diff < 0)
        return -1;
    else if (Diff > 0)
        return 1;

    // Diff is 0, so sort now on RssDelta

    Diff = ((tQueue*)p2)->RssDelta - ((tQueue*)p1)->RssDelta;
    if (Diff < 0)
        return -1;
    else if (Diff > 0)
        return 1;

    // Diff is still 0, so sort now on iParent

    int iDiff = ((tQueue*)p1)->iParent - ((tQueue*)p2)->iParent;
    if (iDiff < 0)
        return -1;
    else if (iDiff > 0)
        return 1;
    return 0;
}

static void AddTermToQ(
    const int iTerm,        // in
    const int nTerms,       // in
    const double RssDelta,  // in
    const bool Sort,        // in
    const int nMaxTerms,    // in
    const int FastK,        // in
    const double FastBeta)  // in: ageing Coef, 0 is no ageing, FastMARS recommends 1
{
    ASSERT(iTerm < nMaxTerms);
    Q[iTerm].nTerms = nTerms;
    Q[iTerm].RssDelta = max(Q[iTerm].RssDelta, RssDelta);
    nQMax = max(nQMax, iTerm+1);
    if (Sort) {
        memcpy(SortedQ, Q, nQMax * sizeof(tQueue));
        if (FastK < nMaxTerms) { // for (small gain in) efficiency, only sort if needed
            qsort(SortedQ, nQMax, sizeof(tQueue), CompareQ);            // sort on RssDelta
            if (FastBeta != 0) {
                for (int iRank = 0; iRank < nQMax; iRank++)
                    SortedQ[iRank].AgedRank =
                        iRank + FastBeta * (nTerms - SortedQ[iRank].nTerms);
                qsort(SortedQ, nQMax, sizeof(tQueue), CompareAgedQ);    // sort on aged rank
            }
        }
    }
}

static void UpdateRssDeltaInQ(const int iParent, const int nTerms, const double RssDelta)
{
    ASSERT(iParent == Q[iParent].iParent);
    ASSERT(iParent < nQMax);
    Q[iParent].nTerms = nTerms;
    Q[iParent].RssDelta = RssDelta;
}

static int GetNextParent(   // returns -1 if no more parents
    const int nTerms,       // use -1 to init
    const int FastK)
{
    static int iQ;          // index into sorted queue
    int iParent = -1;
    if (nTerms == -1) {     // init?
        iQ = -1;
        if (nTraceGlobal == 5)
            printf("\n|Considering parents ");
    } else {
        iQ++;
        int Max = min(nTerms, FastK);
        if (iQ < Max)
            iParent = SortedQ[iQ].iParent;
        if (nTraceGlobal == 5 && iParent >= 0)
            printf("%d [%g] ", iParent+IOFFSET, SortedQ[iQ].RssDelta);
    }
    return iParent;
}

#endif // FAST_MARS

//--------------------------------------------------------------------------------------------
// Print a summary of the model, for debug tracing

#if STANDALONE
static void PrintSummary(
        const int    nMaxTerms,     // in
        const int    nTerms,        // in: number of cols in bx, some may not be used
        const int    nPreds,        // in: number of predictors
        const bool   UsedCols[],    // in: specifies used colums in bx
        const int    *nFactorsInTerm, // in: number of hockey stick funcs in basis term
        const double Betas[],       // in: if NULL will print zeroes
        const int    Dirs[],        // in
        const double Cuts[])        // in
{
    printf("   nFacs       Beta\n");

    int iUsed = -1;
    for (int iTerm = 0; iTerm < nTerms; iTerm++) {
        if (UsedCols[iTerm]) {
            iUsed++;
            printf("%2.2d  %2d    %9.3g | ",
                iTerm, nFactorsInTerm[iTerm], (Betas? Betas[iUsed]: 0));
            }
        else
            printf("%2.2d  --    %9.3g | ", iTerm, 0.0);

        int iPred;
        for (iPred = 0; iPred < nPreds; iPred++)
            if (Dirs_(iTerm,iPred) == 0)
                printf(" . ");
            else
                printf("%2d ", Dirs_(iTerm,iPred));

        printf("|");

        for (iPred = 0; iPred < nPreds; iPred++)
            if (Dirs_(iTerm,iPred) == 0)
                printf("    .    ");
            else
                printf("%8.3g ", Cuts_(iTerm,iPred));

        printf("\n");
    }
    printf("\n");
}
#endif // STANDALONE

//--------------------------------------------------------------------------------------------
// The caller has selected a candidate predictor iPred and a candidate iParent.
// This function now selects a knot.
// The new term will be a copy of the parent term but extended by multiplying
// it by a new hockey stick function at the selected knot.
//
// The general idea: scan backwards through all (ordered) knots for the
// given predictor iPred, calculating RssDelta.
// If RssDelta > *pRssDeltaForThisTerm (and all else is ok), then select
// the knot (by updating *piBestCase and *pRssDeltaForThisTerm).
//
// Returns true if a new knot is found.  If so, it will have also
// updated *piBestCase and *pRssDeltaForThisTerm.
//
// There are currently nTerms in the model. We want to add a term pair
// at index nTerms and nTerms+1.
//
// This function must be fast.

static INLINE void FindKnot(
    int    *piBestCase,             // out: possibly updated, index into the ORDERED x's
    double *pRssDeltaForThisTerm,   // io:  possibly updated
    double CovCol[],                // io:  updated, nTerms x 1
    double CovSy[],                 // io:  initialized, nTerms x 1
    double CovSx[],                 // out: scratch buffer, overwritten, nTerms x 1
    const int nTerms,               // in
    const int iParent,              // in: parent term
    const int iPred,                // in: predictor index
    const int nCases,               // in
    const int nMaxTerms,            // in
    const double RssExistingModel,  // in: RSS before adding this new candidate term
    const double RssDeltaBase,      // in: change in RSS if predictor iPred enters linearly
    const double RssDeltaPrev,      // in: change in RSS when previous term was added
    const double bx[],              // in
    const double bxOrth[],          // in
    const double bxOrthCenteredT[], // in
    const double bxOrthMean[],      // in
    const double y[],               // in
    const double x[],               // in
    const int xOrder[],             // in
    const double yMean,             // in
    const double NewVarAdjust)      // in
{
    if (RssExistingModel <= 0)
        error("assertion failed: RssExistingModel <= 0 (y is all const?)");
    ASSERT(RssDeltaPrev > 0);
#if USE_BLAS
    double Dummy = bxOrth[0];       // prevent compiler warning: unused parameter
    Dummy = bxOrthMean[0];
#else
    double Dummy = bxOrthCenteredT[0];
    Dummy = nMaxTerms;
#endif
    const int nMinSpan = GetMinSpan(nCases, bx, iParent, nMinSpanGlobal);
    const int nEndSpan = GetEndSpan(nCases);

    CovSy[nTerms] = 0;
    memset(CovCol, 0, (nTerms+1) * sizeof(double));
    memset(CovSx,  0, (nTerms+1) * sizeof(double));
    double bxSum = 0, bx2Sum = 0, bx2xSum = 0, bxxSum = 0, ybxSum = 0, st = 0;

    for (int iCase = nCases - 2; iCase >= 0; iCase--) { // -2 allows for ix1
        const int    ix0 = xOrder_(iCase,  iPred);      // get the x's in descending order
        const double x0  = x_(ix0,iPred);
        const int    ix1 = xOrder_(iCase+1,iPred);
        const double x1  = x_(ix1,iPred);
        const double bx1 = bx_(ix1,iParent);
        const double xDelta = x1 - x0;
        const double bx2 = sq(bx1);
#if USE_BLAS
        daxpy_(&nTerms, &bx1, &bxOrthCenteredT_(0,ix1), &ONE, CovSx,  &ONE);
        daxpy_(&nTerms, &xDelta, CovSx, &ONE, CovCol, &ONE);
#else
        int it;
        for (it = 0; it < nTerms; it++) {
            CovSx[it]  += (bxOrth_(ix1,it) - bxOrthMean[it]) * bx1;
            CovCol[it] += xDelta * CovSx[it];
        }
#endif
        bxSum   += bx1;
        bx2Sum  += bx2;
        bx2xSum += bx2 * x1;
        bxxSum  += bx1 * x1;
        const double su = st;
        st = bxxSum - bxSum * x0;
        CovCol[nTerms] += xDelta * (2 * bx2xSum - bx2Sum * (x0 + x1)) +
                          (sq(su) - sq(st)) / nCases;
        ybxSum         += (y[ix1] - yMean) * bx1;
        CovSy[nTerms]  += xDelta * ybxSum;

        if (iCase % nMinSpan == 0 && CovCol[nTerms] > 0) {
            // calculate RssDelta and see if this knot beats the previous best
#if USE_BLAS
            double temp1 = CovSy[nTerms]  - ddot_(&nTerms, CovSy,  &ONE, CovCol, &ONE);
            double temp2 = CovCol[nTerms] - ddot_(&nTerms, CovCol, &ONE, CovCol, &ONE);
#else
            double temp1 = CovSy[nTerms];
            double temp2 = CovCol[nTerms];
            for (it = 0; it < nTerms; it++) {
                temp1 -= CovSy[it] * CovCol[it];
                temp2 -= sq(CovCol[it]);
            }
#endif
            if (temp2/CovCol[nTerms] > BX_TOL) {    // $$ HastieTibs code note says fix this
                                                    // but it seems to work well enough

                const double RssDelta = NewVarAdjust * (RssDeltaBase + sq(temp1) / temp2);

                // $$ HastieTibs code has an extra test here
                // !(iCase > 0 && x_(ix0,iPred) == x_(xOrder_(iCase-1,iPred),iPred))

                if (RssDelta > *pRssDeltaForThisTerm        &&
                        RssDelta < 1.01 * RssExistingModel  &&
                        RssDelta < 2 * RssDeltaPrev         &&
                        iCase >= nEndSpan                   &&
                        iCase < nCases - nEndSpan           &&
                        bx1 > 0) {
                    *piBestCase = iCase;
                    *pRssDeltaForThisTerm = RssDelta;
                }
            }
        }
    }
}

//--------------------------------------------------------------------------------------------
// The caller has selected a candidate parent term iParent.
// This function now selects a predictor.
//
// $$ These functions have a ridiculous number of parameters, I know.

static INLINE void FindPred(
    int    *piBestCase,             // out: return -1 if no new term available
                                    //      else return an index into the ORDERED x's
    int    *piBestPred,             // out
    int    *piBestParent,           // out: existing term on which we are basing the new term
    double *pBestRssDelta,          // out: adding new term reduces RSS this much
                                    //      set to 0 if no possible new term
    double *pBestRssDeltaForThisTerm, // out
    bool   *pIsNewForm,             // out
    bool   *pIsLinPred,             // out: true if knot is at min x val so x enters linearly
    double xbx[],                   // io
    double CovSx[],                 // io
    double CovCol[],                // io
    double CovSy[],                 // io
    double bxOrth[],                // io
    double bxOrthCenteredT[],       // io
    double bxOrthMean[],            // io
    const int iParent,              // in
    const double x[],               // in
    const double y[],               // in
    const int nCases,               // in
    const int nPreds,               // in
    const int nTerms,               // in
    const int nMaxTerms,            // in
    const double yMean,             // in
    const double RssExistingModel,  // in: RSS before adding this new candidate term
    const double RssDeltaPrev,      // in: change in RSS caused by adding previous term
    const double bx[],              // in
    const bool FullSet[],           // in
    const int xOrder[],             // in
    const int nUses[],              // in
    const int Dirs[],               // in
    const double NewVarPenalty)     // in: penalty for adding a new variable
{
#if USING_R
    ServiceR(nCases);
#endif
    for (int iPred = 0; iPred < nPreds; iPred++) {
        if (Dirs_(iParent,iPred) != 0) {    // predictor is in parent term?
            if (nTraceGlobal >= 6)
                printf("|Parent %-2d Pred %-2d                                   "
                    "                skip (pred is in parent)\n",
                    iParent+IOFFSET, iPred+IOFFSET);
        } else {
            const double NewVarAdjust = 1 + (nUses[iPred] == 0? 0: NewVarPenalty);
            if (nTraceGlobal >= 6)
                printf("|Parent %-2d Pred %-2d ", iParent+IOFFSET, iPred+IOFFSET);

            double RssDeltaBase = 0;    // change in RSS for iParent and new knot at iPred

            bool IsNewForm = GetNewFormFlag(iPred, iParent, Dirs,
                                        FullSet, nTerms, nPreds, nMaxTerms);
            if (IsNewForm) {
                // Add a candidate term at bx[,nTerms], with predictor iPred entering
                // linearly. Do this by setting the knot at the lowest value xMin of x,
                // since min(0,x-xMin)==x-xMin for all x.  The change in RSS caused by
                // adding this term forms the base RSS delta which we will try to beat
                // in the search in FindKnot.

                // set xbx to x * bx[,iParent]

                int iCase;
                for (iCase = 0; iCase < nCases; iCase++)
                    xbx[iCase] = x_(iCase,iPred) * bx_(iCase,iParent);

                // init bxOrth[,nTerms] and bxOrthMean[nTerms] for the candidate term
                // $$ look into IsNewForm handling here, it's confusing

                InitBxOrthCol(bxOrth, bxOrthCenteredT, bxOrthMean, &IsNewForm,
                    xbx, nTerms, FullSet, nCases, nMaxTerms, iParent, iPred);

                // init CovCol and CovSy[nTerms]

                memset(CovCol, 0, (nTerms-1) * sizeof(double));
                CovCol[nTerms] = 1;
                CovSy[nTerms] = 0;
                for (iCase = 0; iCase < nCases; iCase++)
                    CovSy[nTerms] += (y[iCase] - yMean) * bxOrth_(iCase,nTerms);

                // calculate change to RSS caused by adding candidate new term

                RssDeltaBase = 0;
                for (iCase = 0; iCase < nCases; iCase++)
                    RssDeltaBase += y[iCase] * bxOrth_(iCase,nTerms);
                RssDeltaBase = NewVarAdjust * sq(RssDeltaBase);

                if (nTraceGlobal >= 6)
                    printf("Case %4d Cut % 12.4g< "
                        "RssDelta %-12.5g ",
                        0+IOFFSET, GetCut(0, iPred, nCases, x, xOrder), RssDeltaBase);

                if (RssDeltaBase > *pBestRssDeltaForThisTerm)
                    *pBestRssDeltaForThisTerm = RssDeltaBase;
                if (RssDeltaBase > *pBestRssDelta) {
                    // The new term (with predictor entering linearly) beats other
                    // candidate terms so far.

                    *pBestRssDelta = RssDeltaBase;
                    *pIsLinPred    = true;
                    *piBestCase    = 0;         // knot is at the lowest value of x
                    *piBestPred    = iPred;
                    *piBestParent  = iParent;
                    if (nTraceGlobal >= 6)
                        printf("NEW");
                }
            }   // is NewForm
            if (nTraceGlobal >= 6)
                printf("\n");

            int iBestCase;
            FindKnot(&iBestCase, pBestRssDeltaForThisTerm, CovCol, CovSy, CovSx,
                    (IsNewForm? nTerms + 1: nTerms),
                    iParent, iPred, nCases, nMaxTerms,
                    RssExistingModel, RssDeltaBase, RssDeltaPrev,
                    bx, bxOrth, bxOrthCenteredT, bxOrthMean,
                    y, x, xOrder, yMean, NewVarAdjust);

            if (*pBestRssDeltaForThisTerm > *pBestRssDelta) {
                *pBestRssDelta = *pBestRssDeltaForThisTerm;
                *pIsLinPred    = false;
                *pIsNewForm    = IsNewForm;
                *piBestCase    = iBestCase;
                *piBestPred    = iPred;
                *piBestParent  = iParent;
                if (nTraceGlobal >= 6)
                    printf("|                  "
                        "Case %4d Cut % 12.4g  RssDelta %-12.5g NEW\n",
                        iBestCase+IOFFSET,
                        GetCut(iBestCase, iPred, nCases, x, xOrder),
                        *pBestRssDelta);
            }
        } // if Dirs
    } // for iPred
}

//--------------------------------------------------------------------------------------------
// Find a new term to add to the model, if possible, and return the
// selected case (i.e. knot), predictor, and parent term indices.
// The new term is a copy of an existing parent term but extended
// by multiplying the parent by a new hockey stick function at the selected knot.
//
// There are currently nTerms in the model. We want to add a term at index nTerms.

static void FindTerm(
    int    *piBestCase,             // out: return -1 if no new term available
                                    //      else return an index into the ORDERED x's
    int    *piBestPred,             // out
    int    *piBestParent,           // out: existing term on which we are basing the new term
    double *pBestRssDelta,          // out: adding new term reduces RSS this much
                                    //      set to 0 if no possible new term
    bool   *pIsNewForm,             // out
    bool   *pIsLinPred,             // out: true if knot is at min x val so x enters linearly
    double bxOrth[],                // io: column nTerms overwritten
    double bxOrthCenteredT[],       // io: kept in sync with bxOrth
    double bxOrthMean[],            // io: element nTerms overwritten
    const double x[],               // in
    const double y[],               // in
    const int nCases,               // in
    const int nPreds,               // in
    const int nTerms,               // in
    const int nMaxDegree,           // in
    const int nMaxTerms,            // in
    const double yMean,             // in
    const double RssExistingModel,  // in: RSS before adding this new candidate term
    const double RssDeltaPrev,      // in: change in RSS caused by adding previous term
    const double bx[],              // in
    const bool FullSet[],           // in
    const int xOrder[],             // in
    const int nFactorsInTerm[],     // in
    const int nUses[],              // in
    const int Dirs[],               // in
    const int FastK,                // in: Fast MARS K
    const double NewVarPenalty)     // in: penalty for adding a new variable
{
#if !FAST_MARS
    int Dummy = FastK;                  // prevent compiler warning: unused parameter
    Dummy = 0;
#endif
    if (nTraceGlobal >= 6)
        printf("\n|Searching for new term %-3d                    RssDelta 0\n",
            nTerms+IOFFSET);

    *piBestCase = -1;
    *pBestRssDelta = 0;
    *pIsLinPred = false;
    int iCase;

    xbx = (double *)malloc1(nCases * sizeof(double));
    CovSx  = (double *)malloc1(nMaxTerms * sizeof(double));
    CovCol = (double *)calloc1(nMaxTerms, sizeof(double));
    CovSy  = (double *)calloc1(nMaxTerms, sizeof(double));

    for (int iTerm = 0; iTerm < nTerms; iTerm++)
        for (iCase = 0; iCase < nCases; iCase++)
            CovSy[iTerm] += (y[iCase] - yMean) * bxOrth_(iCase,iTerm);

    int iParent;
#if FAST_MARS
    GetNextParent(-1, FastK);   // init queue iterator
    while ((iParent = GetNextParent(nTerms, FastK)) > -1) {
#else
    for (iParent = 0; iParent < nTerms; iParent++) {
#endif
        // Assume a bad RssDelta for iParent.
        // This pushes parents that can't be used to the bottom of the queue.
        double BestRssDeltaForThisTerm = -1;    // only used by FAST_MARS

        if (nFactorsInTerm[iParent] >= nMaxDegree) {
            if (nTraceGlobal >= 6)
                printf("|Parent %-2d                                                      "
                    "     skip (nFactorsInTerm %d)\n",
                    iParent+IOFFSET, nFactorsInTerm[iParent]);
        } else
            FindPred(piBestCase, piBestPred, piBestParent,
                pBestRssDelta, &BestRssDeltaForThisTerm, pIsNewForm, pIsLinPred,
                xbx, CovSx, CovCol, CovSy, bxOrth, bxOrthCenteredT, bxOrthMean,
                iParent, x, y, nCases, nPreds, nTerms, nMaxTerms, yMean, RssExistingModel,
                RssDeltaPrev, bx, FullSet, xOrder, nUses, Dirs, NewVarPenalty);
#if FAST_MARS
        UpdateRssDeltaInQ(iParent, nTerms, BestRssDeltaForThisTerm);
#endif
    } // iParent
#if FAST_MARS
    if (nTraceGlobal >= 5)
        printf("\n");
#else
    if (nTraceGlobal >= 6)
        printf("\n");
#endif
    free1(CovSy);
    free1(CovCol);
    free1(CovSx);
    free1(xbx);
}

//--------------------------------------------------------------------------------------------
static void PrintForwardProlog(const double RssNull,
                const int nCases,
                const int nPreds,
                const char *sPredNames[])       // in: predictor names, can be NULL
{
    if (nTraceGlobal == 1)
        printf("Term 1");
    else if (nTraceGlobal >= 2) {
        printf("Forward pass: model matrix %d x %d minspan %d endspan %d\n\n",
            nCases, nPreds,
            GetMinSpan(nCases, NULL, 0, nMinSpanGlobal), GetEndSpan(nCases));

        printf("         GRSq    RSq   DeltaRSq         RSS Pred ");
        if (sPredNames)
            printf("PredName  ");
        printf("       Cut  Terms   Parents\n");

        printf("1      0.0000 0.0000           %12.4g                  %s%d\n",
            RssNull, (sPredNames? "          ":""), IOFFSET);
    }
}

//--------------------------------------------------------------------------------------------
static void PrintForwardStep(
    const int nTerms,
    const int nUsedTerms,
    const int iBestCase,
    const int iBestPred,
    const double Rss,
    const double RSq,
    const double RSqDelta,
    const double Gcv,
    const double GcvNull,
    const int nCases,
    const int xOrder[],
    const double x[],
    const bool IsLinPred,
    const bool IsNewForm,
    const char *sPredNames[])   // in: predictor names, can be NULL
{
    if (nTraceGlobal == 1) {
        printf(", ");
        if (nTerms % 30 == 29)
            printf("\n");
        printf("%d", nTerms+IOFFSET);
    } else if (nTraceGlobal >= 2) {
        printf("%-4d%9.4g %6.4f %12.4g %9g  %3d",
            nTerms+IOFFSET, 1-Gcv/GcvNull, RSq, RSqDelta, Rss, iBestPred+IOFFSET);

        if (sPredNames) {
            if (sPredNames[iBestPred] && sPredNames[iBestPred][0])
                printf(" %8.8s ", sPredNames[iBestPred]);
            else
                printf(" %8.8s ", " ");
        }
        if (iBestPred < 0)
            printf("                   ");
        else {
            if (iBestCase == -1)
                printf("       none  ");
            else
                printf("% 11.5g%c ",
                    GetCut(iBestCase, iBestPred, nCases, x, xOrder), (IsLinPred? '<': ' '));
            if (!IsLinPred && IsNewForm)  // two new used terms?
                printf("%-3d %-3d ", nUsedTerms-2+IOFFSET, nUsedTerms-1+IOFFSET);
            else
                printf("%-3d     ", nUsedTerms-1+IOFFSET);
        }
        fflush(stdout);
    }
}

//--------------------------------------------------------------------------------------------
static void PrintForwardEpilog(
            const int nTerms,  int nMaxTerms,
            const double Thresh,
            const double RSq, const double RSqDelta,
            const double Gcv, const double GcvNull,
            const bool FullSet[])
{
    double GRSq = 1-Gcv/GcvNull;
    if (nTraceGlobal == 1)
        printf(" GRSq: %.4g RSq: %.4g ", GRSq, RSq);
    if (nTraceGlobal >= 1)
        printf("\n");
    if (nTraceGlobal >= 2) {
        // treat very low nMaxTerms as a special case
        // because RSDelta etc. not yet completely initialized

        if (nMaxTerms < 3) {
            printf("Reached max number of terms %d", nMaxTerms);
            if (nTerms < nMaxTerms)
                printf(" (no room for another term pair)");
            printf("\n");
        }
        // NOTE: this code must match the loop termination conditions in ForwardPass

        else if (GRSq < MIN_GRSQ)
            printf("Reached min GRSq (GRSq %g < %g)\n", GRSq, MIN_GRSQ);

        else if (RSqDelta < Thresh)
            printf("Reached delta RSq threshold (DeltaRSq %g < %g)\n",
                RSqDelta, Thresh);

        else if (RSq > 1-Thresh)
            printf("Reached max RSq (RSq %g > %g)\n", RSq, 1-Thresh);

        else {
            printf("Reached max number of terms %d", nMaxTerms);
            if (nTerms < nMaxTerms)
                printf(" (no room for another term pair)");
            printf("\n");
        }
        printf("Forward pass complete: %d terms", nTerms);
        int nUsed = GetNbrUsedCols(FullSet, nMaxTerms);
        if (nUsed != nTerms)
            printf(" (%d terms used)", nUsed);
        if (nTraceGlobal >= 2)
            printf("\n");
        if (nTraceGlobal >= 3)
            printf("\n");
    }
}

//--------------------------------------------------------------------------------------------
// Forward pass
//
// After initializing the intercept term, the main for loop adds terms in pairs.
// In the for loop, nTerms is the index of the potential new term; nTerms+1
// the index of its partner.
// The upper term in the term pair may not be useable.  If so we still
// increment nTerms by 2 but don't set the flag in FullSet.
//
// $$ feature: allow nMaxDegree to be a vector of length nPreds
// $$ feature: add option to prescale x and y

static void ForwardPass(
    int    *pnTerms,            // out: highest used term number in full model
    bool   FullSet[],           // out: 1 * nMaxTerms, indices of lin indep cols of bx
    double bx[],                // out: nCases * nMaxTerms
    int    Dirs[],              // out: nMaxTerms * nPreds, 1,0,-1 for iTerm, iPred
    double Cuts[],              // out: nMaxTerms * nPreds, cut for iTerm, iPred
    int    nFactorsInTerm[],    // out: number of hockey stick funcs in each MARS term
    int    nUses[],             // out:
    const double x[],           // in:  1 x nCases
    const double y[],           // in:  1 x nCases
    const int nCases,           // in:
    const int nPreds,           // in:
    const int nMaxDegree,       // in:
    const int nMaxTerms,        // in:
    const double Penalty,       // in:
    double Thresh,              // in: forward step threshold
    int FastK,                  // in: Fast MARS K
    const double FastBeta,      // in: Fast MARS ageing coef
    const double NewVarPenalty, // in: penalty for adding a new variable
    const char *sPredNames[])   // in: predictor names, can be NULL
{
    if (nTraceGlobal >= 5)
        printf("earth.c %s\n", VERSION);

    // The limits below are somewhat arbitrary and generous.
    // They are intended to catch gross errors on the part of the
    // callee, and to prevent crashes because of 0 sizes etc.
    // We use error rather than ASSERT because these are user settable params
    // and we want to be informative from the user's perspective.
    // The errors are reported using the variable names in the R code.

    // prevent possible minspan range problems, also prevent crash when nCases==0
    if (nCases < 10)
        error("need at least 10 rows in x, you have %d", nCases);
    if (nCases < nPreds)    // (this check may not actually be necessary)
        warning("Need as many rows as columns in x");
    if (nMaxDegree <= 0)
        error("degree %d <= 0", nMaxDegree);
    if (nMaxDegree >= 50)
        error("degree %d >= 50", nMaxDegree);
    if (nMaxTerms < 3)      // prevent internal misbehaviour
        error("nk %d < 3", nMaxTerms);
    if (nMaxTerms > 10000)
        error("nk %d > 10000", nMaxTerms);
    if (Penalty < 0 && Penalty != -1)
        error("penalty %g < 0, the only legal value less than 0 is -1, "
            "meaning terms and knots are free", Penalty);
    if (Penalty > 1000)
        error("penalty %g > 1000", Penalty);
    if (Thresh < 0)
        error("thresh %g < 0", Thresh);
    if (Thresh >= 1)
        error("thresh %g >= 1", Thresh);
    if (Thresh < 1e-10)     // needed for numerical stability, 1e-10 seems ok
        Thresh = 1e-10;
    if (nMinSpanGlobal < -1000)
        error("minspan %d < -1000", nMinSpanGlobal);
    if (nMinSpanGlobal > nCases/2)
        error("minspan %d > nrow(x)/2 %d", nMinSpanGlobal, nCases/2);
    if (FastK == -1)
        FastK = 10000;      // bigger than any nMaxTerms
    if (FastK < 3)
        error("fast.k %d < 3, the only legal value less than 3 is -1, ",
            "meaning no Fast MARS", FastK);
    if (FastBeta < 0)
        error("fast.beta %g < 0", FastBeta);
    if (FastBeta > 1000)
        error("fast.beta %g > 1000", FastBeta);
    if (nTraceGlobal < 0)
        warning("trace.earth %d < 0", nTraceGlobal);
    if (nTraceGlobal > 7)
        warning("trace.earth %d > 7", nTraceGlobal);
    if(NewVarPenalty < 0)
        warning("newvar.penalty %g < 0", NewVarPenalty);
    if(NewVarPenalty > 1)
        warning("newvar.penalty %g > 1", NewVarPenalty);

    bxOrth          = (double *)malloc1(nCases * nMaxTerms * sizeof(double));
    bxOrthCenteredT = (double *)malloc1(nMaxTerms * nCases * sizeof(double));
    bxOrthMean      = (double *)malloc1(nMaxTerms * sizeof(double));

    memset(FullSet,        0, nMaxTerms * sizeof(bool));
    memset(Dirs,           0, nMaxTerms * nPreds * sizeof(int));
    memset(Cuts,           0, nMaxTerms * nPreds * sizeof(double));
    memset(nFactorsInTerm, 0, nMaxTerms * sizeof(int));
    memset(nUses,          0, nPreds * sizeof(int));
    memset(bx, 0, nCases * nMaxTerms * sizeof(double));
    xOrder = OrderArray(x, nCases, nPreds);
#if BETA_CACHE
    InitBetaCache(nMaxTerms, nPreds);
#endif
    FullSet[0] = true;  // intercept
    bool GoodCol;
    InitBxOrthCol(bxOrth, bxOrthCenteredT, bxOrthMean, &GoodCol,   // intercept col 0
        &bx_(0,0), 0 /*nTerms*/, FullSet, nCases, nMaxTerms, -1, -1);
    ASSERT(GoodCol);
    for (int iCase = 0; iCase < nCases; iCase++)
        bx_(iCase,0) = 1;
    double yMean = Mean(y, nCases);
    double RssNull = SumOfSquares(y, yMean, nCases);        // intercept only model
    double Rss = RssNull, RssDelta = RssNull, RSq = 0, RSqDelta = 0;
    int nUsedTerms = 1;     // number of used basis terms including intercept, for GCV calc
    double Gcv = 0, GcvNull = GetGcv(1, nCases, RssNull, Penalty);
    PrintForwardProlog(RssNull, nCases, nPreds, sPredNames);
#if FAST_MARS
    InitQ(nMaxTerms);
    AddTermToQ(0, 1, RssNull, true, nMaxTerms, FastK, FastBeta);
#endif
    int nTerms, iBestCase;
    for (nTerms = 1;                                    // start after intercept
            nTerms < nMaxTerms-1 && RSq < 1-Thresh;     // -1 allows for upper term in pair
            nTerms += 2) {                              // add terms in pairs
        int iBestPred, iBestParent;
        bool IsNewForm, IsLinPred;

        FindTerm(&iBestCase, &iBestPred, &iBestParent, &RssDelta, &IsNewForm, &IsLinPred,
            bxOrth, bxOrthCenteredT, bxOrthMean, x, y,
            nCases, nPreds, nTerms, nMaxDegree, nMaxTerms,
            yMean, Rss, RssDelta, bx, FullSet, xOrder, nFactorsInTerm, nUses, Dirs,
            FastK, NewVarPenalty);

        nUsedTerms++;
        if (!IsLinPred && IsNewForm)    // add paired term too?
            nUsedTerms++;
        Rss -= RssDelta;
        if (Rss < 1e-10)                // RSS can go slightly neg due to rounding
            Rss = 0;
        Gcv = GetGcv(nUsedTerms, nCases, Rss, Penalty);
        const double OldRSq = RSq;
        RSq = 1-Rss/RssNull;
        RSqDelta = RSq - OldRSq;

        PrintForwardStep(nTerms, nUsedTerms, iBestCase, iBestPred, Rss, RSq, RSqDelta,
            Gcv, GcvNull, nCases, xOrder, x, IsLinPred, IsNewForm, sPredNames);

        // no explicit test for iBestCase<0, because if iBestCase<0 then RSqDelta==0

        if ((1-Gcv/GcvNull) < MIN_GRSQ || RSqDelta < Thresh) {
            if (nTraceGlobal >= 2)
                printf("reject term\n");
            break;                      // NOTE break
        }
        AddTermPair(Dirs, Cuts, bx, bxOrth, bxOrthCenteredT, bxOrthMean,
            FullSet, nFactorsInTerm, nUses,
            nTerms, iBestParent, iBestCase, iBestPred, nPreds, nCases,
            nMaxTerms, IsNewForm, IsLinPred, x, xOrder);
#if FAST_MARS
        AddTermToQ(nTerms,   nTerms, POS_INF, false, nMaxTerms, FastK, FastBeta);
        //$$ should only add 2nd term if it is a good term?
        AddTermToQ(nTerms+1, nTerms, POS_INF, true,  nMaxTerms, FastK, FastBeta);
        if (nTraceGlobal >= 7)
            PrintSortedQ(FastK);
#endif
        if (nTraceGlobal >= 2)
            printf("\n");
    }
    PrintForwardEpilog(nTerms, nMaxTerms, Thresh, RSq, RSqDelta, Gcv, GcvNull, FullSet);
    *pnTerms = nTerms;
#if BETA_CACHE
    FreeBetaCache();
#endif
#if FAST_MARS
    FreeQ();
#endif
    free1(xOrder);
    free1(bxOrthMean);
    free1(bxOrthCenteredT);
    free1(bxOrth);
}

//--------------------------------------------------------------------------------------------
// This is an interface from R to the C routine ForwardPass

#if USING_R
void ForwardPassR(              // for use by R
    int    FullSet[],           // out: 1 x nMaxTerms, bool vec of lin indep cols of bx
    double bx[],                // out: nCases x nMaxTerms
    double Dirs[],              // out: nMaxTerms x nPreds, elements are 1,0,-1
    double Cuts[],              // out: nMaxTerms x nPreds, cut for iTerm,iPred
    const double x[],           // in:  1 x nCases
    const double y[],           // in:  1 x nCases
    const int *pnCases,         // in:
    const int *pnPreds,         // in:
    const int *pnMaxDegree,     // in:
    const int *pnMaxTerms,      // in:
    const double *pPenalty,     // in:
    double *pThresh,            // in: forward step threshold
    const int *pnMinSpan,       // in:
    const int *pFastK,            // in: Fast MARS K
    const double *pFastBeta,      // in: Fast MARS ageing coef
    const double *pNewVarPenalty, // in: penalty for adding a new variable
    const int *pnTrace,           // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char *sPredNames[])     // in: predictor names in trace printfs, can be NULL
{
    nTraceGlobal = *pnTrace;
    nMinSpanGlobal = *pnMinSpan;

    const int nMaxTerms = *pnMaxTerms;
    const int nPreds = *pnPreds;
    const int nCases = *pnCases;

    // nUses is the number of time a predictor is used in the model
    nUses = (int *)malloc1(*pnPreds * sizeof(int));

    // nFactorsInTerm is number of hockey stick functions in basis term
    nFactorsInTerm = (int *)malloc1(nMaxTerms * sizeof(int));

    iDirs = (int *)calloc1(nMaxTerms * nPreds, sizeof(int));

    // convert int to bool (may be redundant, depending on compiler)
    BoolFullSet = (int *)malloc1(*pnMaxTerms * sizeof(bool));
    int iTerm;
    for (iTerm = 0; iTerm < *pnMaxTerms; iTerm++)
        BoolFullSet[iTerm] = FullSet[iTerm];

    int nTerms;
    ForwardPass(&nTerms, BoolFullSet, bx, iDirs, Cuts, nFactorsInTerm, nUses,
            x, y, nCases, nPreds, *pnMaxDegree, nMaxTerms,
            *pPenalty, *pThresh, *pFastK, *pFastBeta, *pNewVarPenalty, sPredNames);

    // remove linearly independent rows if necessary -- this updates BoolFullSet
    //
    // $$ memory optimization fix: memory peaks in the call to RegressAndFix I think.
    //    Should compact used cols of bx here so don't need extra copy created by
    //    CopyUsedCols in Regress(), will also be faster, less thrashing.

    double Rss;
    double *Betas = (double *)malloc1(nMaxTerms * sizeof(double));
    double *Residuals = (double *)malloc1(nCases * sizeof(double));
    RegressAndFix(Betas, Residuals, &Rss, BoolFullSet, bx, y, nCases, nMaxTerms);
    free1(Residuals);
    free1(Betas);

    for (int iTerm = 0; iTerm < nMaxTerms; iTerm++)     // convert int to double
        for (int iPred = 0; iPred < nPreds; iPred++)
            Dirs[iTerm + iPred * nMaxTerms] =
                iDirs[iTerm + iPred * nMaxTerms];

    for (iTerm = 0; iTerm < *pnMaxTerms; iTerm++)       // convert bool to int
        FullSet[iTerm] = BoolFullSet[iTerm];

    free1(iDirs);
    free1(nFactorsInTerm);
    free1(nUses);
}
#endif // USING_R

//--------------------------------------------------------------------------------------------
// Prune BestSet by deleting terms without increasing the GCV, if possible.
//
// This uses a brute force method: it does a full regression each time it
// temporarily removes a candidate term.  There are much quicker methods.
// $$ use Alan Miller's Fortran routines

#if STANDALONE
static void BackwardPass(
    bool   BestSet[],           // io: specifies which cols in bx are used
    double *pGcvBest,           // out
    double Residuals[],         // out
    double Betas[],             // out
    const int    nCases,        // in
    const int    nMaxTerms,     // in
    const double Penalty,       // in
    const double bx[],          // in
    const double y[])           // in
{
    double RssNull = SumOfSquares(y, Mean(y, nCases), nCases); // intercept only model
    double GcvNull = GetGcv(1, nCases, RssNull, Penalty);
    double Rss;
    RegressAndFix(Betas, Residuals, &Rss, BestSet, bx, y, nCases, nMaxTerms);
    int nUsedTerms = GetNbrUsedCols(BestSet, nMaxTerms);    // includes intercept
    double GcvBest = GetGcv(nUsedTerms, nCases, Rss, Penalty);
    double RssBest = Rss;
    if (nTraceGlobal >= 3) {
        printf("Backward pass:\nDelTerm     RSq      GRSq\n");
        printf("        %7.4f %9.4f max\n", 1- Rss/RssNull, 1-GcvBest/GcvNull);
    }
    WorkingSet = (bool *)malloc1(nMaxTerms * sizeof(bool));
    memcpy(WorkingSet, BestSet, nMaxTerms * sizeof(bool));
    int *Pivots = (int *)malloc1(nMaxTerms * sizeof(int));
    int nRank;

    for (int iTerm = nUsedTerms-1; iTerm > 0; iTerm--) {
        // set iDelete to the best term for deletion

        int iDelete = -1;
        double RssTestBest = POS_INF;

        int iUsedTerm = 1;
        for (int iTestTerm = 1; iTestTerm < nMaxTerms; iTestTerm++)
            if (WorkingSet[iTestTerm]) {
                WorkingSet[iTestTerm] = false;
                Regress(Betas, Residuals, &Rss, &nRank, Pivots, bx, y,
                    nCases, nMaxTerms, WorkingSet);
                ASSERT(nRank == iTerm);
                WorkingSet[iTestTerm] = true;
                iUsedTerm++;
                if (Rss < RssTestBest) {        // new min?
                    iDelete = iTestTerm;
                    RssTestBest = Rss;
                }
            }
        ASSERT(iDelete > 0);

        // delete term at iDelete from the WorkingSet and possibly update BestSet

        WorkingSet[iDelete] = false;
        char *sMax = "";
        double GcvTestBest = GetGcv(iTerm, nCases, RssTestBest, Penalty);
        if (GcvTestBest <= GcvBest) {           // new min?
            sMax = " max";
            GcvBest = GcvTestBest;
            RssBest = RssTestBest;
            memcpy(BestSet, WorkingSet, nMaxTerms * sizeof(bool));
        }
        if (nTraceGlobal >= 3)
            printf("   %4d %7.4f %9.4f%s\n",
                iDelete+IOFFSET, 1-RssTestBest/RssNull, 1-GcvTestBest/GcvNull, sMax);

        Regress(Betas, Residuals, &Rss, &nRank, Pivots, bx, y, nCases, nMaxTerms, WorkingSet);
    }
    int nUsedTerms1 = GetNbrUsedCols(BestSet, nMaxTerms);
    if (nTraceGlobal >= 3)
        printf("\n");
    if (nTraceGlobal >= 1)
        printf(
            "Backward pass complete: %d terms (pruned %d term%s), GRSq: %.4g RSq: %.4g\n\n",
            nUsedTerms1, nUsedTerms-nUsedTerms1, (nUsedTerms-nUsedTerms1 == 1? "": "s"),
            1-GcvBest/GcvNull, 1-RssBest/RssNull);

    // get Betas etc. for the final set BestSet of terms

    Regress(Betas, Residuals, &Rss, &nRank, Pivots, bx, y, nCases, nMaxTerms, BestSet);
    ASSERT(nRank == nUsedTerms1);
    *pGcvBest = GcvBest;
    free1(WorkingSet);
    free1(Pivots);
}
#endif // STANDALONE

//--------------------------------------------------------------------------------------------
#if STANDALONE
static int DiscardUnusedTerms(
    double bx[],             // io: nCases x nMaxTerms
    int    Dirs[],           // io: nMaxTerms x nPreds
    double Cuts[],           // io: nMaxTerms x nPreds
    bool   FullSet[],        // io
    int    nFactorsInTerm[], // io
    const int nMaxTerms,
    const int nPreds,
    const int nCases)
{
    int nUsed = 0, iTerm;
    for (iTerm = 0; iTerm < nMaxTerms; iTerm++)
        if (FullSet[iTerm]) {
            memcpy(bx + nUsed * nCases, bx + iTerm * nCases, nCases * sizeof(double));
            for (int iPred = 0; iPred < nPreds; iPred++) {
                Dirs_(nUsed, iPred) = Dirs_(iTerm, iPred);
                Cuts_(nUsed, iPred) = Cuts_(iTerm, iPred);
            }
            nFactorsInTerm[nUsed] = nFactorsInTerm[iTerm];
            nUsed++;
        }
    memset(FullSet, 0, nMaxTerms * sizeof(bool));
    for (iTerm = 0; iTerm < nUsed; iTerm++)
        FullSet[iTerm] = true;
    return nUsed;
}
#endif

//--------------------------------------------------------------------------------------------
#if STANDALONE
void Earth(
    double bx[],            // out: nCases x nMaxTerms
    double *pGcvBest,       // out:
    bool   BestSet[],       // out: nMaxTerms x 1, indices of best set of cols of bx
    int    *pnTerms,        // out: max term nbr in final model, after removing lin dep terms
    int    Dirs[],          // out: nMaxTerms x nPreds, 1,0,-1 for term iTerm, predictor iPred
    double Cuts[],          // out: nMaxTerms x nPreds, cut for term iTerm, predictor iPred
    double Residuals[],     // out: nCases x 1
    double Betas[],         // out: nMaxTerms x 1
    const double x[],       // in:  nCases x nPreds
    const double y[],       // in:  nCases x 1
    const int nCases,       // in:
    const int nPreds,       // in:
    const int nMaxDegree,   // in: Friedman's mi
    const int nMaxTerms,    // in: includes the intercept term
    const double Penalty,   // in:
    double Thresh,          // in: forward step threshold
    const int nMinSpan,     // in: set to non zero to override internal calculation
    const bool Prune,       // in: do backward pass (current implementation is slow)
    const int FastK,        // in: Fast MARS K
    const double FastBeta,  // in: Fast MARS ageing coef
    const double NewVarPenalty, // in: penalty for adding a new variable
    const int nTrace,           // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char *sPredNames[])   // in: predictor names in trace printfs, can be NULL
{
    nTraceGlobal = nTrace;
    nMinSpanGlobal = nMinSpan;

    // nUses is the number of time a predictor is used in the model
    nUses = (int *)malloc1(nPreds * sizeof(int));

    // nFactorsInTerm is number of hockey stick functions in basis term
    nFactorsInTerm = (int *)malloc1(nMaxTerms * sizeof(int));

    int nTerms;
    ForwardPass(&nTerms, BestSet, bx, Dirs, Cuts, nFactorsInTerm, nUses,
        x, y, nCases, nPreds, nMaxDegree, nMaxTerms,
        Penalty, Thresh, FastK, FastBeta, NewVarPenalty, sPredNames);

    if (!Prune || nTraceGlobal >= 6) {   // nTraceGlobal check because PrintSummary below
        // get Residuals and Betas, and ensure bx is full rank by updating BestSet
        double Rss;
        RegressAndFix(Betas, Residuals, &Rss, BestSet, bx, y, nCases, nMaxTerms);
    }
    nTerms = DiscardUnusedTerms(bx, Dirs, Cuts, BestSet, nFactorsInTerm,
                                nMaxTerms, nPreds, nCases);
    if (nTraceGlobal >= 6)
        PrintSummary(nMaxTerms, nTerms, nPreds, BestSet, nFactorsInTerm, Betas, Dirs, Cuts);
    if (Prune) {
        BackwardPass(BestSet, pGcvBest, Residuals, Betas, nCases, nTerms, Penalty, bx, y);
        if (nTraceGlobal > 6)
            PrintSummary(nMaxTerms, nTerms, nPreds, BestSet, nFactorsInTerm, Betas, Dirs, Cuts);
    }
    *pnTerms = nTerms;
    free1(nFactorsInTerm);
    free1(nUses);
}
#endif // STANDALONE

//--------------------------------------------------------------------------------------------
// return the max number of knots in any term

#if STANDALONE
static int GetMaxKnotsPerTerm(
    const bool   UsedCols[],    // in
    const int    Dirs[],        // in
    const int    nPreds,        // in
    const int    nTerms,        // in
    const int    nMaxTerms)     // in
{
    int nKnotsMax = 0;
    for (int iTerm = 1; iTerm < nTerms; iTerm++)
        if (UsedCols[iTerm]) {
            int nKnots = 0; // number of knots in this term
            for (int iPred = 0; iPred < nPreds; iPred++)
                if (Dirs_(iTerm, iPred) != 0)
                    nKnots++;
            if (nKnots > nKnotsMax)
                nKnotsMax = nKnots;
        }
    return nKnotsMax;
}
#endif // STANDALONE

//--------------------------------------------------------------------------------------------
// print a string representing the earth expresssion, one term per line
// $$ spacing is not quite right and is overly complicated

#if STANDALONE
void FormatEarth(
    const bool   UsedCols[],// in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],    // in: nMaxTerms x nPreds, 1,0,-1 for term iTerm, predictor iPred
    const double Cuts[],    // in: nMaxTerms x nPreds, cut for term iTerm, predictor iPred
    const double Betas[],   // in: nMaxTerms x 1
    const int    nPreds,
    const int    nTerms,
    const int    nMaxTerms,
    const int    nDigits)   // number of significant digits to print
{
    int iBestTerm = 0;
    int nKnotsMax = GetMaxKnotsPerTerm(UsedCols, Dirs, nPreds, nTerms, nMaxTerms);
    int nKnots = 0;
    char s[1000];
    ASSERT(nDigits >= 0);
    char sFormat[50];  sprintf(sFormat,  "%%-%d.%dg", nDigits+6, nDigits);
    char sFormat1[50]; sprintf(sFormat1, "%%%d.%dg",  nDigits+6, nDigits);
    int nPredWidth;
    if (nPreds > 100)
        nPredWidth = 3;
    else if (nPreds > 10)
        nPredWidth = 2;
    else
        nPredWidth = 1;
    char sPredFormat[20]; sprintf(sPredFormat, "%%%dd", nPredWidth);
    char sPad[500]; sprintf(sPad, "%*s", 28+nDigits+nPredWidth, " ");    // comment pad

    printf(sFormat, Betas[0]);      // intercept
    while (nKnots++ < nKnotsMax)
        printf(sPad);
    printf(" // 0\n");

    for (int iTerm = 1; iTerm < nTerms; iTerm++)
        if (UsedCols[iTerm]) {
            iBestTerm++;
            printf("%+-9.3g", Betas[iBestTerm]);
            nKnots = 0;
            for (int iPred = 0; iPred < nPreds; iPred++) {
                switch(Dirs_(iTerm, iPred)) {
                    case  0:
                        break;
                    case -1:
                        sprintf(s, " * max(0, %s - %*sx[%s])",
                            sFormat, nDigits+2, " ", sPredFormat);
                        printf(s, Cuts_(iTerm, iPred), iPred);
                        nKnots++;
                        break;
                    case  1:
                        sprintf(s, " * max(0, x[%s]%*s-  %s)",
                            sPredFormat,  nDigits+2, " ", sFormat1);
                        printf(s, iPred, Cuts_(iTerm, iPred));
                        nKnots++;
                        break;
                    default:
                        ASSERT(false);
                        break;
                }
            }
            while (nKnots++ < nKnotsMax)
                printf(sPad);
            printf(" // %d\n", iBestTerm);
        }
}
#endif // STANDALONE

//--------------------------------------------------------------------------------------------
// return the value predicted by an earth model, given  a vector of inputs x

#if STANDALONE
double PredictEarth(
    const double x[],           // in: vector nPreds x 1 of input values
    const bool   UsedCols[],    // in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],        // in: nMaxTerms x nPreds, 1,0,-1 for iTerm iPred
    const double Cuts[],        // in: nMaxTerms x nPreds, cut for term iTerm predictor iPred
    const double Betas[],       // in: nMaxTerms x 1
    const int    nPreds,
    const int    nTerms,
    const int    nMaxTerms)
{
    double yHat = Betas[0];     // intercept
    int iTerm11 = 0;
    for (int iTerm = 1; iTerm < nTerms; iTerm++)
        if (UsedCols[iTerm]) {
            iTerm11++;
            double Term = Betas[iTerm11];
            for (int iPred = 0; iPred < nPreds; iPred++)
                switch(Dirs_(iTerm, iPred)) {
                    case  0: break;
                    case -1: Term *= max(0, Cuts_(iTerm, iPred) - x[iPred]); break;
                    case  1: Term *= max(0, x[iPred] - Cuts_(iTerm, iPred)); break;
                    default: ASSERT("bad direction" == NULL); break;
                }
            yHat += Term;
        }
    return yHat;
}
#endif // STANDALONE

//--------------------------------------------------------------------------------------------
// example main routine

#if STANDALONE && MAIN
void error(const char *args, ...)
{
    char s[1000];
    va_list p;
    va_start(p, args);
    vsprintf(s, args, p);
    va_end(p);
    printf("\nError: %s\n", s);
    exit(-1);
}

void xerbla_(char *srname, int *info)
{
    char buf[7];
    strncpy(buf, srname, 6);
    buf[6] = '\0';
    error("BLAS/LAPACK routine %6s gave error code %d", buf, -(*info));
}

int main(void)
{
    const int nTrace = 3;
    const int nMaxTerms = 21;
    const int nCases = 100;
    const int nPreds = 1;

    double *x = (double *)malloc1(nCases * nPreds * sizeof(double));
    double *y = (double *)malloc1(nCases *          sizeof(double));
    for (int i = 0; i < nCases; i++) {
        double x0 = (double)i / nCases;
        x[i] = x0;
        y[i] = sin(4 * x0);
    }
    double *bx        = (double *)malloc1(nCases    * nMaxTerms * sizeof(double));
    bool   *BestSet   = (bool *)  malloc1(nMaxTerms *             sizeof(bool));
    int    *Dirs      = (int *)   malloc1(nMaxTerms * nPreds *    sizeof(int));
    double *Cuts      = (double *)malloc1(nMaxTerms * nPreds *    sizeof(double));
    double *Residuals = (double *)malloc1(nCases    *             sizeof(double));
    double *Betas     = (double *)malloc1(nMaxTerms *             sizeof(double));
    double GcvBest;
    int    nTerms;

    Earth(bx, &GcvBest, BestSet, &nTerms, Dirs, Cuts, Residuals, Betas,
        x, y, nCases, nPreds, 1, nMaxTerms, 2, 0.001, 0, true, 20, 0, 0, nTrace, NULL);

    printf("Expression:\n");
    FormatEarth(BestSet, Dirs, Cuts, Betas, nPreds, nTerms, nMaxTerms, 3);

    free1(x);
    free1(y);
    free1(bx);
    free1(BestSet);
    free1(Dirs);
    free1(Cuts);
    free1(Residuals);
    free1(Betas);

    return 0;
}
#endif

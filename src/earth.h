// earth.h: externs for earth.c
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// A copy of the GNU General Public License is available at
// http://www.r-project.org/Licenses

#if !defined(EARTH_H)
#define EARTH_H

void FreeEarth(void);

#if USING_R

SEXP ForwardPassR(              // for use by R
    SEXP SEXP_FullSet,          // out: nMaxTerms x 1, bool vec of lin indep cols of bx
    SEXP SEXP_bx,               // out: MARS basis matrix, nCases x nMaxTerms
    SEXP SEXP_Dirs,             // out: nMaxTerms x nPreds, elements are -1,0,1,2
    SEXP SEXP_Cuts,             // out: nMaxTerms x nPreds, cut for iTerm,iPred
    SEXP SEXP_iTermCond,        // out: reason we terminated the forward pass
    SEXP SEXP_x,                // in: nCases x nPreds, unweighted x
    SEXP SEXP_y,                // in: nCases x nResp, unweighted but scaled y
    SEXP SEXP_yw,               // in: nCases x nResp, weighted and scaled y
    SEXP SEXP_WeightsArg,       // in: nCases x 1, never R_NilValue
    SEXP SEXP_nCases,           // in: number of rows in x and elements in y
    SEXP SEXP_nResp,            // in: number of cols in y
    SEXP SEXP_nPreds,           // in: number of cols in x
    SEXP SEXP_nMaxDegree,       // in:
    SEXP SEXP_Penalty,          // in:
    SEXP SEXP_nMaxTerms,        // in:
    SEXP SEXP_Thresh,           // in: forward step threshold
    SEXP SEXP_nMinSpan,         // in:
    SEXP SEXP_nEndSpan,         // in:
    SEXP SEXP_nFastK,           // in: Fast MARS K
    SEXP SEXP_FastBeta,         // in: Fast MARS ageing coef
    SEXP SEXP_NewVarPenalty,    // in: penalty for adding a new variable (default is 0)
    SEXP SEXP_LinPreds,         // in: nPreds x 1, 1 if predictor must enter linearly
    SEXP SEXP_Allowed,          // in: constraints function, can be R NULL
    SEXP SEXP_nAllowedArgs,     // in: number of arguments to Allowed function, 3...5
    SEXP SEXP_Env,              // in: environment for Allowed function
    SEXP SEXP_AdjustEndSpan,    // in:
    SEXP SEXP_nAutoLinPreds,    // in: assume predictor linear if knot is min predictor value
    SEXP SEXP_nUseBetaCache,    // in: 1 to use the beta cache, for speed
    SEXP SEXP_Trace,            // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    SEXP SEXP_sPredNames);      // in: predictor names in trace printfs

void EvalSubsetsUsingXtxR(      // for use by R
    double        PruneTerms[], // out: specifies which cols in bx are in best set
    double        RssVec[],     // out: nTerms x 1
    const int*    pnCases,      // in
    const int*    pnResp,       // in: number of cols in y
    const int*    pnMaxTerms,   // in
    const double  bx[],         // in: MARS basis matrix, all cols must be indep
    const double  y[],          // in: nCases * nResp (possibly weighted)
    const double* pTrace);      // in

void RegressR(                  // for testing earth routine Regress from R
    double       Betas[],       // out: (nUsedCols+1) * nResp, +1 is for intercept
    double       Residuals[],   // out: nCases * nResp
    double       Rss[],         // out: RSS, summed over all nResp
    double       Diags[],       // out: diags of inv(transpose(x) * x)
    int*         pnRank,        // out: nbr of indep cols in x
    int          iPivots[],     // out: nCols
    const double x[],           // in: nCases x nCols
    const double y[],           // in: nCases x nResp
    const int*   pnCases,       // in: number of rows in x and in y
    const int*   pnResp,        // in: number of cols in y
    int*         pnCols,        // in: number of columns in x, some may not be used
    const int    UsedColsR[]);  // in: specifies used columns in x (assume R LOGICAL is stored as int)
#endif // USING_R

#if STANDALONE

void Earth(
    double* pBestGcv,            // out: GCV of the best model i.e. BestSet columns of bx
    int*    pnTerms,             // out: max term nbr in final model, after removing lin dep terms
    int*    piTermCond,          // out: reason we terminated the foward pass
    bool    BestSet[],           // out: nMaxTerms x 1, indices of best set of cols of bx
    double  bx[],                // out: nCases x nMaxTerms
    int     Dirs[],              // out: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
    double  Cuts[],              // out: nMaxTerms x nPreds, cut for iTerm, iPred
    double  Residuals[],         // out: nCases x nResp
    double  Betas[],             // out: nMaxTerms x nResp
    const double x[],            // in: nCases x nPreds
    const double y[],            // in: nCases x nResp
    const double WeightsArg[],   // in: nCases x 1, can be NULL, not yet supported
    const size_t nCases,         // in: number of rows in x and elements in y
    const int nResp,             // in: number of cols in y
    const int nPreds,            // in: number of cols in x
    const int nMaxDegree,        // in: Friedman's mi
    const int nMaxTerms,         // in: includes the intercept term
    const double Penalty,        // in: GCV penalty per knot
    const double Thresh,         // in: forward step threshold
    const int nMinSpan,          // in: set to non zero to override internal calculation
    const int nEndSpan,          // in: set to non zero to override internal calculation
    const bool Prune,            // in: do backward pass
    const int nFastK,            // in: Fast MARS K
    const double FastBeta,       // in: Fast MARS ageing coef
    const double NewVarPenalty,  // in: penalty for adding a new variable
    const int LinPreds[],        // in: nPreds x 1, 1 if predictor must enter linearly
    const double AdjustEndSpan,  // in: for adjusting endspan for interaction terms
    const bool AutoLinPreds,     // in: assume predictor linear if knot is max predictor value
    const bool UseBetaCache,     // in: 1 to use the beta cache, for speed
    const double Trace,          // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char* sPredNames[]);   // in: predictor names in trace printfs, can be NULL

void FormatEarth(
    const bool   UsedCols[], // in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],     // in: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
    const double Cuts[],     // in: nMaxTerms x nPreds, cut for iTerm, iPred
    const double Betas[],    // in: nMaxTerms x nResp
    const int    nPreds,
    const int    nResp,      // in: number of cols in y
    const int    nTerms,
    const int    nMaxTerms,
    const int    nDigits,    // number of significant digits to print
    const double MinBeta);   // terms with fabs(betas) less than this are not printed, 0 for all

void PredictEarth(
    double       y[],        // out: vector nResp
    const double x[],        // in: vector nPreds x 1 of input values
    const bool   UsedCols[], // in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],     // in: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
    const double Cuts[],     // in: nMaxTerms x nPreds, cut for iTerm, iPred
    const double Betas[],    // in: nMaxTerms x nResp
    const int    nPreds,     // in: number of cols in x
    const int    nResp,      // in: number of cols in y
    const int    nTerms,
    const int    nMaxTerms);

#endif // STANDALONE

#endif // EARTH_H

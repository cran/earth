// earth.h: externs for earth.c

#if !defined(EARTH_H)
#define EARTH_H

#if USING_R
void FreeR(void);
#endif

#if USING_R
void ForwardPassR(
    int    FullSet[],           // out: 1 x nMaxTerms, bool vec of lin indep cols of bx
    double bx[],                // out: MARS basis matrix, nCases x nMaxTerms
    double Dirs[],              // out: nMaxTerms x nPreds, elements are 1,0,-1
    double Cuts[],              // out: nMaxTerms x nPreds, cut for iTerm,iPred
    const double x[],           // in: nCases x nPreds
    const double y[],           // in: nCases x nResponses
    const int *pnCases,         // in: number of rows in x and elements in y
    const int *pnResponses,     // in: number of cols in y
    const int *pnPreds,         // in: number of cols in x
    const int *pnMaxDegree,     // in:
    const int *pnMaxTerms,      // in:
    const double *pPenalty,     // in:
    double *pThresh,            // in: forward step threshold
    const int *pnMinSpan,       // in:
    const int *pFastK,            // in: Fast MARS K
    const double *pFastBeta,      // in: Fast MARS ageing coef
    const double *pNewVarPenalty, // in: penalty for adding a new variable
    const int *pnTrace,           // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char *sPredNames[]);    // in: predictor names in trace printfs, can be NULL
#endif

#if USING_R
void EvalSubsetsUsingXtxR(
    double       PruneTerms[],  // out: specifies which cols in bx are in best set
    double       RssVec[],      // out: nTerms x 1
    const int    *pnCases,      // in
    const int    *pnResponses,  // in: number of cols in y
    const int    *pnMaxTerms,   // in
    const double bx[],          // in: MARS basis matrix, all cols must be independent
    const double y[]);          // in: nCases * nResponses
#endif

#if STANDALONE
void Earth(
    double *pBestGcv,       // out: GCV of the best model i.e. BestSet columns of bx
    int    *pnTerms,        // out: max term nbr in final model, after removing lin dep terms
    bool   BestSet[],       // out: nMaxTerms x 1, indices of best set of cols of bx
    double bx[],            // out: nCases x nMaxTerms
    int    Dirs[],          // out: nMaxTerms x nPreds, 1,0,-1 for term iTerm, predictor iPred
    double Cuts[],          // out: nMaxTerms x nPreds, cut for term iTerm, predictor iPred
    double Residuals[],     // out: nCases x nResponses
    double Betas[],         // out: nMaxTerms x nResponses
    const double x[],       // in: nCases x nPreds
    const double y[],       // in: nCases x nResponses
    const int nCases,       // in: number of rows in x and elements in y
    const int nResponses,   // in: number of cols in y
    const int nPreds,       // in: number of cols in x
    const int nMaxDegree,   // in: Friedman's mi
    const int nMaxTerms,    // in: includes the intercept term
    const double Penalty,   // in: GCV penalty per knot
    double Thresh,          // in: forward step threshold
    const int nMinSpan,     // in: set to non zero to override internal calculation
    const bool Prune,       // in: do backward pass (current implementation is slow)
    const int FastK,        // in: Fast MARS K
    const double FastBeta,  // in: Fast MARS ageing coef
    const double NewVarPenalty, // in: penalty for adding a new variable
    const int nTrace,           // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char *sPredNames[]);  // in: predictor names in trace printfs, can be NULL

void FormatEarth(
    const bool   UsedCols[],// in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],    // in: nMaxTerms x nPreds, 1,0,-1 for term iTerm, predictor iPred
    const double Cuts[],    // in: nMaxTerms x nPreds, cut for term iTerm, predictor iPred
    const double Betas[],   // in: nMaxTerms x nResponses
    const int    nPreds,
    const int    nResponses,// in: number of cols in y
    const int    nTerms,
    const int    nMaxTerms,
    const int    nDigits,   // number of significant digits to print
    const double MinBeta);  // terms with abs(beta) less than this are not printed, 0 for all

void PredictEarth(
    double       y[],           // out: vector nResponses
    const double x[],           // in: vector nPreds x 1 of input values
    const bool   UsedCols[],    // in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],        // in: nMaxTerms x nPreds, 1,0,-1 for iTerm iPred
    const double Cuts[],        // in: nMaxTerms x nPreds, cut for term iTerm predictor iPred
    const double Betas[],       // in: nMaxTerms x nResponses
    const int    nPreds,
    const int    nResponses,    // in: number of cols in y
    const int    nTerms,
    const int    nMaxTerms);
#endif

#endif // EARTH_H

// earth.h: externs for earth.c

#if !defined(EARTH_H)
#define EARTH_H

#if USING_R
void FreeR(void);               // for use by R
#endif

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
    const char *sPredNames[]);    // in: predictor names in trace printfs, can be NULL
#endif

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
    const char *sPredNames[]);  // in: predictor names in trace printfs, can be NULL

void FormatEarth(
    const bool   UsedCols[],// in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],    // in: nMaxTerms x nPreds, 1,0,-1 for term iTerm, predictor iPred
    const double Cuts[],    // in: nMaxTerms x nPreds, cut for term iTerm, predictor iPred
    const double Betas[],   // in: nMaxTerms x 1
    const int    nPreds,
    const int    nTerms,
    const int    nMaxTerms,
    const int    nDigits);  // number of significant digits to print

double PredictEarth(
    const double x[],           // in: vector nPreds x 1 of input values
    const bool   UsedCols[],    // in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],        // in: nMaxTerms x nPreds, 1,0,-1 for iTerm iPred
    const double Cuts[],        // in: nMaxTerms x nPreds, cut for term iTerm predictor iPred
    const double Betas[],       // in: nMaxTerms x 1
    const int    nPreds,
    const int    nTerms,
    const int    nMaxTerms);
#endif

#endif // EARTH_H

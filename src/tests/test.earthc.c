// test.earth.c: main() for testing earth c routines
// Comments containing "$$" mark known issues

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <crtdbg.h>
#include "earth.h"

#define PRINT_TIME 0

#define sq(x)   ((x) * (x))

//-----------------------------------------------------------------------------
void error(const char *args, ...)
{
    char *s = (char *)malloc(1000);
    va_list p;
    va_start(p, args);
    vsprintf(s, args, p);
    va_end(p);
    printf("\nError: %s\n", s);
    free(s);
    exit(-1);
}

//-----------------------------------------------------------------------------
static double getGcv(const int nTerms, // nbr basis terms including intercept
                const int nCases, double Rss, const double penalty)
{
    double cost;
    if (penalty < 0)    // special case: terms and knots are free
        cost = 0;
    else
        cost = nTerms + (penalty * (double)(nTerms-1) / 2);   // nKnots = (nTerms-1)/2

    return Rss / (nCases * sq(1 - cost/nCases));
}

//-----------------------------------------------------------------------------
static double RandUniform(void) // uniform rand number from -1 to +1
{
    return (double)((rand() % 20000) - 10000) / 10000;
}

//-----------------------------------------------------------------------------
static double RandGauss(void)   // standard normal random number
{
    double r = 0;
    for (int i = 0; i < 12; i++)    // by central limit theorem sum of uniforms is gaussian
        r += RandUniform();
    return r / 2;
}

//-----------------------------------------------------------------------------
static double funcNoise(const double x[]) { return RandGauss(); }

static double func0(const double x[]) { return x[0]; }

static double func1(const double x[]) { return x[0] + x[1] + .1 * RandGauss(); }

static double func2(const double x[]) { return x[0] + x[1] + x[0]*x[1]; }

static double func3(const double x[]) { return cos(x[0]) + x[1]; }

static double func4(const double x[]) { return sin(2 * x[0]) + 2*x[1] + 0.5*x[0]*x[1]; }

static double func5(const double x[]) { return x[0] + x[1] + x[3] + x[4] + x[5] + x[1]*x[2] + (x[3]+1)*(x[4]+1)*x[5]; }

static double func4lin(const double x[]) { return x[0] + x[1] + x[3] + x[4]; }

static double func6(const double x[]) {                 // 5 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5]
            + .1 * RandGauss();
}

static double func6clean(const double x[]) {            // 5 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5];
}

static double func7(const double x[]) {                 // 10 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +x[6] +x[7] +x[8] +x[9] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5] + x[6]*x[7] + x[8]*x[9];
}

static double func8(const double x[]) {                 // 20 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +x[6] +x[7] +x[8] +x[9] +
           x[10]+x[11]+x[12]+x[13]+x[14]+x[15]+x[16]+x[17]+x[18]+x[19] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5] + x[6]*x[7] + x[8]*x[9] +
            + .1 * RandGauss();
}

static double func56(const double x[]) {    // Friedman MARS paper eqn 56
    return 0.1 * exp(4*x[0]) + 4 / (1 + exp(-20*(x[1]-0.5)) + 3*x[2] + 2*x[3] + x[4] + RandGauss());
}

//-----------------------------------------------------------------------------
static void testEarth(char sTestName[],
                double (__cdecl *func)(const double xVec[]),    // function to be modeled
                const int nCases, const int nPreds,
                const int nMaxDegree, const int nMaxTerms,
                const int nTrace, const bool Format,
                const double ForwardStepThresh,
                const int K, const double AgeCoeff, const double NewVarPenalty,
                const int seed)
{
    double *x          = (double *)malloc(nCases * nPreds * sizeof(double));
    double *y          = (double *)malloc(nCases * sizeof(double));
    bool   *BestSet    = (bool *)  malloc(nMaxTerms * sizeof(bool));
    int    *Dirs       = (int *)   malloc(nMaxTerms * nPreds * sizeof(int));
    double *Cuts       = (double *)malloc(nMaxTerms * nPreds * sizeof(double));
    double *Residuals  = (double *)malloc(nCases * sizeof(double));
    double *Betas      = (double *)malloc(nMaxTerms * sizeof(double));
    double *bx         = (double *)malloc(nCases * nMaxTerms * sizeof(double));

    static int nTest;
    nTest++;

    printf("=============================================================================\n");
    printf("TEST %d: %s\n", nTest, sTestName);

    // init x and y

    srand(seed);
    double *xVec = (double *)malloc(nPreds * sizeof(double));
    int iCase;
    static int xxx;
    for (iCase = 0; iCase < nCases; iCase++) {
        for (int iPred = 0; iPred < nPreds; iPred++) {
            double x1;
            x1 = (double)((rand() % 20000) - 10000) / 10000;    // rand number from -1 to +1
            x[iCase + iPred * nCases] = x1;
            xVec[iPred] = x[iCase + iPred * nCases];
        }
        y[iCase] = func(xVec);
    }
    xxx++;
    free(xVec);

    double BestGcv;
    int nTerms;
    const double Penalty = ((nMaxDegree>1)? 3:2);
    clock_t Time = clock();

    Earth(bx, &BestGcv, BestSet, &nTerms, Dirs, Cuts, Residuals, Betas,
        x, y, nCases, nPreds, nMaxDegree, nMaxTerms, Penalty, ForwardStepThresh,
        0, true, // MinSpan Prune
        K, AgeCoeff, NewVarPenalty, nTrace, NULL);

    // calc nUsedTerms

    int nUsedTerms = 0;
    for (int iTerm = 0; iTerm < nTerms; iTerm++)
        if (BestSet[iTerm])
            nUsedTerms++;

    // calc RSquared, GRSquared

    double Rss = 0, Tss = 0, meanY = 0;
    for (iCase = 0; iCase < nCases; iCase++)
        meanY += y[iCase];
    meanY /= nCases;
    xVec = (double *)malloc(nPreds * sizeof(double));
    for (iCase = 0; iCase < nCases; iCase++) {
        for (int iPred = 0; iPred < nPreds; iPred++)
            xVec[iPred] = x[iCase + iPred * nCases];
        double yHat = PredictEarth(xVec, BestSet, Dirs, Cuts, Betas, nPreds, nTerms, nMaxTerms);
        double Residual = y[iCase] - yHat;
        Rss += sq(Residual);
        Tss += sq(y[iCase] - meanY);
    }
    free(xVec);
    const double RSq =  1 - Rss/Tss;
    const double GcvNull =  getGcv(1, nCases, Tss, Penalty);
    const double GRSq =  1 - getGcv(nUsedTerms, nCases, Rss, Penalty) / GcvNull;

    // show results

    printf("RESULT %d: GRSq %g RSq %g nTerms %d of %d of %d [%.2f secs]\n", nTest,
        GRSq, RSq, nUsedTerms, nTerms, nMaxTerms,
#if PRINT_TIME
        (double)(clock() - Time) / CLOCKS_PER_SEC);
#else
        99.99);
#endif
    if (Format && nTrace) {
        printf("\nTEST %d: FUNCTION\n", nTest);
        FormatEarth(BestSet, Dirs, Cuts, Betas, nPreds, nTerms, nMaxTerms, 3);
        printf("\n");
    }
    free(x);
    free(y);
    free(BestSet);
    free(Dirs);
    free(Cuts);
    free(Residuals);
    free(Betas);
    free(bx);
}

//-----------------------------------------------------------------------------
int main(void)
{                                               // nCases  nPreds  nMaxDegree nMaxTerms nTrace  Format   Thresh K AgeCoeff seed
  testEarth("noise N=1000",            funcNoise,     1000,      1,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("x0 N=10",                     func0,       10,      1,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("x0 N=1000",                   func0,     1000,      1,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("x0 + noise N=1000",           func0,     1000,    1+1,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("x0 + x1 N=1000",              func1,     1000,      2,          2,       11,     6,  true,   0.001, 20, 0, 0, 99);
  testEarth("x0 + x1 + noise N=1000",      func1,     1000,    2+8,          2,       51,     0,  true,   0.001, 20, 0, 0, 99);
  testEarth("x0 + x1 + x0*x1 N=30",        func2,       30,      2,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("x0 + x1 + x0*x1 N=1000",      func2,     1000,      2,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("cos(x0) + x1 N=1000",         func3,     1000,      2,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("sin(2*x0)+2*x1*.5*x0*x1",     func4,     1000,      2,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("sin(2*x0)+2*x1*.5*x0*x1",     func4,     1000,      3,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("3rd order, mi=2 ni=11",       func5,     1000,      6,          2,       11,     1,  true,   0.001, 20, 0, 0, 99);
  testEarth("3rd order, mi=2 ni=51",       func5,     1000,      6,          2,       51,     2,  true,   0.001, 20, 0, 0, 99);
  testEarth("3rd order, mi=3",             func5,     1000,      6,          3,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("5 preds + noise",             func6,      200,   5+10,          2,      101,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("5 preds clean",          func6clean,      200,   5+10,          2,      101,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("10 preds + noise",            func7,      200,  10+40,          2,      101,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("20 preds + noise, N=100",     func8,      100,  20+10,          2,      101,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("20 preds + noise, N=400",     func8,      400,  20+10,          2,      101,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("3rd order, mi=3, +noise",     func5,     1000,     10,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("eqn56 mi=1 N=100",            func56,     100,     10,          1,      101,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("eqn56 mi=2 N=100",            func56,     100,     10,          2,       51,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("eqn56 mi=10 N=100",           func56,     100,     10,         10,       51,     3,  true,   0.001, 20, 0, 0, 99);
//testEarth("eqn56 mi=10 N=1000",          func56,    1000,     10,         10,      101,     3,  true,   0.001, 20, 0, 0, 99);
//testEarth("eqn56 mi=10 N=5000",          func56,    5000,     10,         10,      101,     3,  true,   0.001, 20, 0, 0, 99);
  testEarth("x0 + x1 + x0*x1 N=30",        func2,       30,      2,          2,       51,     7,  true,   0.001, 99, 0, 0, 99);
  testEarth("x0 + x1 + x0*x1 N=30",        func2,       30,      2,          2,       51,     7,  true,   0.001,  4, 0, 0, 99);
  testEarth("x0 + x1 + x0*x1 N=30",        func2,       30,      2,          2,       51,     7,  true,   0.001,  4, 1, 0, 99);
  return 0;
}

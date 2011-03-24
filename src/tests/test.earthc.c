// test.earth.c: main() for testing earth c routines
// Comments containing "TODO" mark known issues

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <crtdbg.h>
#include "../earth.h"

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
static double funcNoise(const double x[], const int iResponse) { return RandGauss(); }

static double func0(const double x[], const int iResponse) { return x[0]; }

static double func1(const double x[], const int iResponse) { return x[0] + x[1] + .1 * RandGauss(); }

static double func2(const double x[], const int iResponse) { return x[0] + x[1] + x[0]*x[1]; }

static double func3(const double x[], const int iResponse) { return cos(x[0]) + x[1]; }

static double func4(const double x[], const int iResponse) { return sin(2 * x[0]) + 2*x[1] + 0.5*x[0]*x[1]; }

static double func5(const double x[], const int iResponse) { return x[0] + x[1] + x[3] + x[4] + x[5] + x[1]*x[2] + (x[3]+1)*(x[4]+1)*x[5]; }

static double func4lin(const double x[], const int iResponse) { return x[0] + x[1] + x[3] + x[4]; }

static double func6(const double x[], const int iResponse) {                 // 5 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5]
            + .1 * RandGauss();
}

static double func6clean(const double x[], const int iResponse) {            // 5 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5];
}

static double func7(const double x[], const int iResponse) {                 // 10 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +x[6] +x[7] +x[8] +x[9] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5] + x[6]*x[7] + x[8]*x[9];
}

static double func8(const double x[], const int iResponse) {                 // 20 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +x[6] +x[7] +x[8] +x[9] +
           x[10]+x[11]+x[12]+x[13]+x[14]+x[15]+x[16]+x[17]+x[18]+x[19] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5] + x[6]*x[7] + x[8]*x[9] +
            + .1 * RandGauss();
}

static double func9(const double x[], const int iResponse) { return x[1]; }

static double func56(const double x[], const int iResponse) {    // Friedman MARS paper eqn 56
    return 0.1 * exp(4*x[0]) + 4 / (1 + exp(-20*(x[1]-0.5)) + 3*x[2] + 2*x[3] + x[4] + RandGauss());
}

// functions for testing multiple responses

static double func0_1(const double x[], const int iResponse) 
{
    if (iResponse == 0)
        return func0(x, iResponse);
    else if (iResponse == 1)
        return func1(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func2_2(const double x[], const int iResponse) 
{
    if (iResponse == 0)
        return func2(x, iResponse);
    else if (iResponse == 1)
        return func2(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func0_4(const double x[], const int iResponse) 
{
    if (iResponse == 0)
        return func0(x, iResponse);
    else  if (iResponse == 1)
        return func4(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func0_2_4(const double x[], const int iResponse) 
{
    if (iResponse == 0)
        return func0(x, iResponse);
    else  if (iResponse == 1)
        return func2(x, iResponse);
    else  if (iResponse == 2)
        return func4(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func2_4_0(const double x[], const int iResponse) 
{
    if (iResponse == 0)
        return func2(x, iResponse);
    else  if (iResponse == 1)
        return func4(x, iResponse);
    else  if (iResponse == 2)
        return func0(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func4_2_0(const double x[], const int iResponse) 
{
    if (iResponse == 0)
        return func4(x, iResponse);
    else  if (iResponse == 1)
        return func2(x, iResponse);
    else  if (iResponse == 2)
        return func0(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func4_6(const double x[], const int iResponse) 
{
    if (iResponse == 0)
        return func4(x, iResponse);
    else  if (iResponse == 1)
        return func6(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

// functions for testing NewVarPenalty

static double func1collinear(const double x[], const int iResponse) 
    { return x[0] + x[1] + .1 * RandGauss(); }

static double func2collinear(const double x[], const int iResponse) 
    { return cos(x[0]) + cos(x[1]) + .1 * RandGauss(); }

//-----------------------------------------------------------------------------
static void TestEarth(char sTestName[],
                double (__cdecl *func)(const double xVec[], const int iResponse),
                const int nCases, const int nresponses, const int nPreds,
                const int nMaxDegree, const int nMaxTerms,
                const int nTrace, const bool Format,
                const double ForwardStepThresh,
                const int K, const double FastBeta, const double NewVarPenalty,
                const int seed, 
                const double Collinear = 0)
{
    #define y_(iCase,iResponse) y[(iCase) + (iResponse)*(nCases)]

    int *LinPreds  = (int *)calloc(nPreds, sizeof(int));

    double *x         = (double *)malloc(nCases    * nPreds     * sizeof(double));
    double *y         = (double *)malloc(nCases    * nresponses * sizeof(double));
    double *bx        = (double *)malloc(nCases    * nMaxTerms  * sizeof(double));
    bool   *BestSet   = (bool *)  malloc(nMaxTerms *              sizeof(bool));
    int    *Dirs      = (int *)   malloc(nMaxTerms * nPreds     * sizeof(int));
    double *Cuts      = (double *)malloc(nMaxTerms * nPreds     * sizeof(double));
    double *Residuals = (double *)malloc(nCases    * nresponses * sizeof(double));
    double *Betas     = (double *)malloc(nMaxTerms * nresponses * sizeof(double));

    static int nTest;
    nTest++;

    printf("=============================================================================\n");
    printf("TEST %d: %s\n", nTest, sTestName);

    // init x and y

    srand(seed);
    double *xVec = (double *)malloc(nPreds * sizeof(double));
    int iCase;
    for (iCase = 0; iCase < nCases; iCase++) {
        for (int iPred = 0; iPred < nPreds; iPred++) {
            double xtemp;
            xtemp = (double)((rand() % 20000) - 10000) / 10000;    // rand number from -1 to +1
            x[iCase + iPred * nCases] = xtemp;
            xVec[iPred] = x[iCase + iPred * nCases];
        }
        if (Collinear > 0) {
            // copy column 0 to 1 with added noise
            double xtemp =  x[iCase] + Collinear * RandGauss();
            x[iCase + 1 * nCases] = xtemp;
            xVec[1] = xtemp;
        }
        for (int iResponse = 0; iResponse < nresponses; iResponse++)
            y_(iCase, iResponse) = func(xVec, iResponse);
    }
    free(xVec);

    double BestGcv;
    int nTerms;
    const double Penalty = ((nMaxDegree>1)? 3:2);
    clock_t Time = clock();

    Earth(&BestGcv, &nTerms, BestSet, bx, Dirs, Cuts, Residuals, Betas,
        x, y, NULL, // weights are NULL
        nCases, nresponses, nPreds, nMaxDegree, nMaxTerms, Penalty, ForwardStepThresh,
        0,      // MinSpan
        true,   // Prune
        K, FastBeta, NewVarPenalty, LinPreds, true, nTrace, NULL);

    // calc nUsedTerms

    int nUsedTerms = 0;
    for (int iTerm = 0; iTerm < nTerms; iTerm++)
        if (BestSet[iTerm])
            nUsedTerms++;

    // calc RSquared, GRSquared

    for (int iResponse = 0; iResponse < nresponses; iResponse++) {
        double Rss = 0, Tss = 0, meanY = 0;
        for (iCase = 0; iCase < nCases; iCase++)
            meanY += y_(iCase, iResponse);
        meanY /= nCases;
        xVec = (double *)malloc(nPreds * sizeof(double));
        double *yHat = (double *)malloc(nresponses * sizeof(double));
        for (iCase = 0; iCase < nCases; iCase++) {
            for (int iPred = 0; iPred < nPreds; iPred++)
                xVec[iPred] = x[iCase + iPred * nCases];
            PredictEarth(yHat, xVec, BestSet, Dirs, Cuts, Betas, nPreds, nresponses, nTerms, nMaxTerms);
            double Residual = y_(iCase,iResponse) - yHat[iResponse];
            Rss += sq(Residual);
            Tss += sq(y_(iCase,iResponse) - meanY);
        }
        free(yHat);
        free(xVec);
        const double RSq =  1 - Rss/Tss;
        const double GcvNull =  getGcv(1, nCases, Tss, Penalty);
        const double GRSq =  1 - getGcv(nUsedTerms, nCases, Rss, Penalty) / GcvNull;

#if PRINT_TIME
        double TimeDelta = (double)(clock() - Time) / CLOCKS_PER_SEC;
#else
        double TimeDelta = 99.99;
#endif
        // show results
        if (nresponses > 1) {
            printf("RESULT %d Response %d: GRSq %g RSq %g nTerms %d of %d of %d",
                nTest, iResponse+1, GRSq, RSq, nUsedTerms, nTerms, nMaxTerms, TimeDelta);
            if (iResponse == 0)
                printf(" [%.2f secs]", TimeDelta);
            printf("\n");
            
        }
        else
            printf("RESULT %d: GRSq %g RSq %g nTerms %d of %d of %d [%.2f secs]\n", 
                nTest, GRSq, RSq, nUsedTerms, nTerms, nMaxTerms, TimeDelta);
    }
    if (Format && nTrace) {
        printf("\nTEST %d: FUNCTION %s\n", nTest, sTestName);
        FormatEarth(BestSet, Dirs, Cuts, Betas, nPreds, nresponses, nTerms, nMaxTerms, 3, 1e-6);
        printf("\n");
    }
    free(LinPreds);
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
{                                        // func    nCases     nResp nPreds  nMaxDegree nMaxTerms nTrace Form Thresh  K B N s
  TestEarth("noise N=1000",            funcNoise,     1000,        1,     1,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("x0 N=10",                     func0,       10,        1,     1,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("x0 N=1000",                   func0,     1000,        1,     1,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("x0 + noise N=1000",           func0,     1000,        1,   1+1,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("x0 + x1 N=1000",              func1,     1000,        1,     2,          2,       11,     7,true,0.001,20,0,0,99);
  TestEarth("x0 + x1 + noise N=1000",      func1,     1000,        1,   2+8,          2,       51,     0,true,0.001,20,0,0,99);
  TestEarth("x0 + x1 + x0*x1 N=30",        func2,       30,        1,     2,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("x0 + x1 + x0*x1 N=1000",      func2,     1000,        1,     2,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("cos(x0) + x1 N=1000",         func3,     1000,        1,     2,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("sin(2*x0)+2*x1*.5*x0*x1",     func4,     1000,        1,     2,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("sin(2*x0)+2*x1*.5*x0*x1",     func4,     1000,        1,     3,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("3rd order, mi=2 ni=11",       func5,     1000,        1,     6,          2,       11,     1,true,0.001,20,0,0,99);
  TestEarth("3rd order, mi=2 ni=51",       func5,     1000,        1,     6,          2,       51,     2,true,0.001,20,0,0,99);
  TestEarth("3rd order, mi=3",             func5,     1000,        1,     6,          3,       51,     3,true,0.001,20,0,0,99);
  TestEarth("5 preds + noise",             func6,      200,        1,  5+10,          2,      101,     3,true,0.001,20,0,0,99);
  TestEarth("5 preds clean",          func6clean,      200,        1,  5+10,          2,      101,     3,true,0.001,20,0,0,99);
  TestEarth("10 preds + noise",            func7,      200,        1, 10+40,          2,      101,     3,true,0.001,20,0,0,99);
  TestEarth("20 preds + noise, N=100",     func8,      100,        1, 20+10,          2,      101,     3,true,0.001,20,0,0,99);
  TestEarth("20 preds + noise, N=400",     func8,      400,        1, 20+10,          2,      101,     3,true,0.001,20,0,0,99);
  TestEarth("3rd order, mi=3, +noise",     func5,     1000,        1,    10,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("eqn56 mi=1 N=100",            func56,     100,        1,    10,          1,      101,     3,true,0.001,20,0,0,99);
  TestEarth("eqn56 mi=2 N=100",            func56,     100,        1,    10,          2,       51,     3,true,0.001,20,0,0,99);
  TestEarth("eqn56 mi=10 N=100",           func56,     100,        1,    10,         10,       51,     3,true,0.001,20,0,0,99);
  // Following two tests are slow so are commented out
//TestEarth("eqn56 mi=10 N=1000",          func56,    1000,        1,    10,         10,      101,     3,true,0.001,20,0,0,99);
//TestEarth("eqn56 mi=10 N=5000",          func56,    5000,        1,    10,         10,      101,     3,true,0.001,20,0,0,99);
  TestEarth("x0 + x1 + x0*x1 N=30",        func2,       30,        1,     2,          2,       51,     7,true,0.001,99,0,0,99);
  TestEarth("x0 + x1 + x0*x1 N=30",        func2,       30,        1,     2,          2,       51,     7,true,0.001, 4,0,0,99);
  TestEarth("x0 + x1 + x0*x1 N=30",        func2,       30,        1,     2,          2,       51,     7,true,0.001, 4,1,0,99);

  // test multiple responses                func    nCases     nResp nPreds  nMaxDegree nMaxTerms nTrace Format   Thresh  K B N s

  TestEarth("x0|x+x1+noise N=30",        func0_1,      100,        2,     2,          1,       51,     4, true,0.001,20,0,0,99);

  TestEarth("x0+x1+x0*x1|x0+x1+x0*x1 degree=1, N=100",        
                                         func2_2,      100,        2,     2,          1,       51,     4, true,0.001,20,0,0,99);

  TestEarth("x0+x1+x0*x1|x0+x1+x0*x1 degree=2 N=100",        
                                         func2_2,      100,        2,     2,          2,       51,     7, true,0.001,20,0,0,99);

  TestEarth("x0|sin(2*x0) + 2*x1 + 0.5*x0*x1 + 8 noise preds, N=50",
                                         func0_4,       50,        2,    10,          2,      101,     3, true,0.001,20,0,0,99);

  TestEarth("x0|x0+x1+x0*x1|sin(2*x0) + 2*x1 + 0.5*x0*x1  + 8 noise preds N=50",
                                         func0_2_4,     50,        3,   3+8,          2,      101,     3, true,0.001,20,0,0,99);

  TestEarth("|x0+x1+x0*x1|sin(2*x0) + 2*x1 + 0.5*x0*x1|x0  + 8 noise preds N=50",
                                         func2_4_0,     50,        3,   3+8,          2,      101,     3, true,0.001,20,0,0,99);

  TestEarth("sin(2*x0) + 2*x1 + 0.5*x0*x1|x0+x1+x0*x1|x0  + 8 noise preds N=50",
                                         func4_2_0,     50,        3,   3+8,          2,      101,     3, true,0.001,20,0,0,99);

  //TODO following gives lousy GRSq for Response 2, investigate
  TestEarth("sin(2*x0) + 2*x1 + 0.5*x0*x1|2nd order 6 preds + noise N=50",
                                         func4_6,     1000,         2,     6,         2,      101,     3, true,0.001,20,0,0,99);

  // test NewVarPenalty                     func    nCases     nResp nPreds  nMaxDegree nMaxTerms nTrace Form Thresh  K B  NP  s Colin

  TestEarth("cos(x1) + cos(x2), x1 and x2 xcollinear, NewVarPenalty=0",
                                   func2collinear,      100,        1,     2,          1,       51,     3,true, 0.001,20,0, 0.0,99,.005);
  TestEarth("cos(x1) + cos(x), x1 and x2 xcollinear, NewVarPenalty=0.05",
                                   func2collinear,      100,        1,     2,          1,       51,     3,true, 0.001,20,0,0.05,99,.005);
  return 0;
}

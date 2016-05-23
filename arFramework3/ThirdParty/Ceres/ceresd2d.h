// Header for ceresd2d.cpp
#include "mex.h"
#include <string.h> // For memcpy // For NULL
#include <typeinfo>
#include "matrix.h"
#include <vector>

#include "ceres/ceres.h"
#include "ceres/loss_function.h"
#include <math.h>
#include <stdio.h>

// Replace Mex inputs with meaningful names
#define parFUN              prhs[0]
#define parX0               prhs[1]
#define parLB               prhs[2]
#define parUB               prhs[3]
#define parOPTS             prhs[4]
#define parPRINTLEVEL       prhs[5]

using namespace std;

using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::LossFunction;
using ceres::TrivialLoss;
using ceres::HuberLoss;
using ceres::SoftLOneLoss;
using ceres::CauchyLoss;
using ceres::ArctanLoss;

void validateInput( const mxArray *prhs[], int nrhs, int *npars, int *bounded );
void loadOptions( const mxArray *prhs[], int* TrustRegionStrategyType, int* DoglegType, int* LossFunctionType, double* LossFunctionVar, double* TolFun, 
                double* TolX, double* TolGradient, int* MaxIter, bool* useNonmonotonicSteps,
                int* maxConsecutiveNonmonotonicSteps, double* maxSolverTimeInSeconds,  int* NumThreads, int* NumLinearSolverThreads,  
                double* InitialTrustRegionRadius, double* MaxTrustRegionRadius, double* MinTrustRegionRadius, double* MinRelativeDecrease, 
                double* MinLMDiagonal, double* MaxLMDiagonal, int* MaxNumConsecutiveInvalidSteps, bool* JacobiScaling, bool* useInnerIterations, 
                double* InnerIterationTolerance, int* LinearSolverType, int* printLevel);
void CallFunctionCheck(const mxArray *prhs[], int p, int* nresidals, int printLevel);
void PrintFunctionInformation(int MaxIter, double TolFun, double TolX, double TolGradient, int bounded, int n, int p);
bool getValueFromStruct( const mxArray *prhs[], const char* fieldName, bool oldValue );
double getValueFromStruct( const mxArray *prhs[], const char* fieldName, double oldValue );
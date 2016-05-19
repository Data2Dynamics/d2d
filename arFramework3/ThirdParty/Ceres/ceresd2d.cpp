// MEX interface for Ceres Solver                        //
// (C) Franz-Georg Wieland 2016                          //
// uses parts of NL2SOL MEX interface by Joep Vanlier    //
// Contact: franz-georg.wieland at mars.uni-freiburg.de  //

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
void PrintFunctionInformation(int MaxIter, int TolFun, int TolX, int TolGradient, int bounded, int n, int p);

class UData {
public:
    UData(int p, const mxArray *prhs[])
    {
        call2Data[1] = mxCreateDoubleMatrix(p, 1, mxREAL);
        call2Data[0] = (mxArray*) parFUN;
    }
    ~UData()
    {
        mxDestroyArray(call2Data[1]);
    }
    
    mxArray *getParameterArray() { return call2Data[1]; };
    mxArray *getFunctionHandle() { return call2Data[0]; };
    mxArray** getCallData() { return call2Data; };    
    int numrhsarguments = 2;     // Number of arguments put in mexCallMatlab (function and parameters)
    int getNumberofArguments() { return numrhsarguments; };
    
    mxArray *call2Data[2]; 
};

///////////// Test CostFunction for simple problem ////////////////////////
// class QuadraticCostFunction
//   : public CostFunction
// {
//  public:
//   QuadraticCostFunction()
//   {
//     mutable_parameter_block_sizes()->push_back(2);            
//     set_num_residuals(2);    
//   }
//   
//   virtual ~QuadraticCostFunction() 
//   {
//   }
// 
//   virtual bool Evaluate(double const* const* parameters,
//                         double* residuals,
//                         double** jacobians) const 
//   {      
//     
//    double q1[] = { parameters[0][0], parameters[0][1] };
// 
//     mexPrintf("\npara1 %e",parameters[0][0]);
//     mexPrintf("\npara2 %e\n",parameters[0][1]);
// 
// 
//     residuals[0] = 10 - q1[0];
//     residuals[1] = -50 - q1[1];
// 
//     if (jacobians != NULL && jacobians[0] != NULL) {
//       jacobians[0][0] = -1;
//       jacobians[0][1] = 0;
//       jacobians[0][2] = 0;
//       jacobians[0][3] = -1;
//     }
// 
//     return true;
//   }
// };
//////////END Test CostFunction for simple problem ////////////////////////


// Ceres CostFuntion class with information to residuals and jacobian
class D2DCostFunction : public CostFunction
{
 public:
      D2DCostFunction(int n,
                      int p,
                      const mxArray *prhs[])     // Constructor of class
      {
            N = n;      // Number of Residuals
            P = p;      // Number of Parameters
            userData = new UData(P, prhs);

            // Set number of residuals according to D2D input
            set_num_residuals(N);

            // Set number of parameters in parameter block 1 according to D2D input
            mutable_parameter_block_sizes()->push_back(P);                                                  
      }

      virtual ~D2DCostFunction()        // Destructor of class
      {
            delete userData;
      }

      virtual bool Evaluate(double const* const* parameters,
                            double* residuals,
                            double** jacobians) const
      {      
            int i;                                      // (Local) iterator
            int l;                                      // (Local) iterator 
            int result;                                 // Temporary result variable         
      
            double *returnValRes;                       // pointer to returned residual 
            double *returnValJac;                       // pointer to returned jacobian 
            mxArray *lhs[2];                            // Temporary varibale for output 
           
            double pm[P];                               // Temporary store for parameter variables
            
            // Read out parameters to temporary store
            for ( i = 0; i < P; i++ )
            {
                pm[i] = parameters[0][i];
            }
            // Insert parameters 
            memcpy( mxGetPr(userData->getParameterArray()), pm, P*sizeof(double) );
            
            // Evaluate function to get resdiuals and jacobian entries
            mexCallMATLAB( 2, lhs, userData->getNumberofArguments(), userData->getCallData(), "feval" );
            
            // Fetch results for residual and jac
            returnValRes = mxGetPr(lhs[0]);
            returnValJac = mxGetPr(lhs[1]);

            // Assign residuals 
            for ( i = 0; i < (N); i++ )
            {
                residuals[i] = returnValRes[i];
            }
            
            // Assign Jacobian elements
            for ( i = 0; i < P; i++ )
            {
                for ( l = 0; l < N; l++ )
                {
                    if (jacobians != NULL && jacobians[0] != NULL)
                    {
                        jacobians[0][l*P+i] = returnValJac[l+(N)*(i)];
                    }
                    
                }
            }

            // Clean up result 
            mxDestroyArray(lhs[0]);
            mxDestroyArray(lhs[1]);   
            return true;
      }
     
  private:  
       int N;                       // Placeholder for number residuals
       int P;                       // Placeholder for number parameters
       UData *userData;             // Class for function handle
};

// The MEX gateway function
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
    
    int     p;                          // Length of parameter vector
    int     n;                          // Length of residuals
    int     bounded;                    // Bounded status (1 = Bounded)                     

    // Data from Options Struct with Default values
    // For detailed documentation see http://ceres-solver.org/
    int        TrustRegionStrategyType = 1;
    int        DoglegType = 1;
    int        LossFunctionType = 1;
    double     LossFunctionVar = 1;
    double     TolFun = 0;
    double     TolX = 1e-6;
    double     TolGradient = 0;
    int        MaxIter = 1000;
    bool       useNonmonotonicSteps = false;
    int        maxConsecutiveNonmonotonicSteps = 5;       
    double     maxSolverTimeInSeconds = 1e6;        
    int        NumThreads = 1;
    int        NumLinearSolverThreads = 1;       
    double     InitialTrustRegionRadius = 1e4;
    double     MaxTrustRegionRadius = 1e16;
    double     MinTrustRegionRadius = 1e-32;
    double     MinRelativeDecrease = 1e-3;
    double     MinLMDiagonal = 1e-8;
    double     MaxLMDiagonal = 1e32;
    int        MaxNumConsecutiveInvalidSteps = 5;
    bool       JacobiScaling = true;
    bool       useInnerIterations = false;
    double     InnerIterationTolerance = 1e-3;
    int        LinearSolverType = 1;
    int        printLevel = 0;                 
                  
       
    double  *tempstor;               // Temporary storage container
    int 	i;                       // Iterator
    
    // Solver output
    double  *x;                      // Output handle 
    double  *xtemp;                  // Temporary storage for inital guess
    double  *finCost;                // Temporary storage for final cost
    double  *iter;                   // Temporary storage for iteration count
    double  *lbtemp;                 // Temporary storage for Lower bound       
    double  *ubtemp;                 // Temporary storage for Upper bound  
    double  *texitflag;
    int     resize;                  // Temporary Checking variable for inital trust region radius
    double  minparasize;             // Temporary variable for minimal 1D-parameter space size
     
   
    // Validate input, check whether we have bounds and determine number of parameters 
    validateInput( prhs, nrhs, &p, &bounded );
    
    
    // Read out Options struct
    if ( nrhs > 3 )
        loadOptions( prhs, &TrustRegionStrategyType, &DoglegType, &LossFunctionType, &LossFunctionVar,           
                     &TolFun, &TolX, &TolGradient, &MaxIter, &useNonmonotonicSteps,
					 &maxConsecutiveNonmonotonicSteps, &maxSolverTimeInSeconds, &NumThreads,
					 &NumLinearSolverThreads, &InitialTrustRegionRadius, &MaxTrustRegionRadius,
					 &MinTrustRegionRadius, &MinRelativeDecrease, &MinLMDiagonal, &MaxLMDiagonal,
					 &MaxNumConsecutiveInvalidSteps, &JacobiScaling, &useInnerIterations, 
					 &InnerIterationTolerance, &LinearSolverType, &printLevel);
        
       
    // Read out printLevel if specifically assigned (overwrites printLevel from options struct)
    if ( nrhs > 5 )
    {
        tempstor = mxGetPr(parPRINTLEVEL);
        printLevel = (int) tempstor[0];
    }
        
    // Bounds vector initialization (needs parameter vector length)
    double lb[p];       // Lower Bounds
    double ub[p];       // Upper Bounds
    
    // Fetch bounds
    if ( bounded == 1)
    {   
        if(printLevel > 0)
        {
            mexPrintf( "Fetching bounds... \n" );
        }
        // Fetch Pointers to Bounds
        lbtemp      = mxGetPr(parLB);
        ubtemp      = mxGetPr(parUB);
        
        // Assign bounds to variables
        for ( i=0; i < p; i++ )
        {
            lb[i] = *(lbtemp+i); 
            ub[i] = *(ubtemp+i);
        }
        if(printLevel > 0)
        {
            mexPrintf( "Bounds fetched\n" ); 
        }   
    }
    
        
    // Call Function to check if it gives out Jacobian
    CallFunctionCheck(prhs, p, &n, printLevel);          
    
    // Create Parameter output handle
    plhs[0]     = mxCreateDoubleMatrix(p,1, mxREAL);
    x           = mxGetPr( plhs[0] );    
    
    // Fetch inital parameter values
    xtemp       = mxGetPr(parX0);
    // Copy initial guess in case everything fails
    memcpy(x, xtemp, p*sizeof(double));    
    
    if ( printLevel > 1 )
        PrintFunctionInformation(MaxIter, TolFun, TolX, TolGradient, bounded, n, p);
    
    
    //////////////////////////////////////////////////////////////////////////////////////
    ///////////////////// CERES PART OF PROGRAM //////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

        
    // Build the problem.
    Problem problem;          
    if(printLevel > 0)
    {
        mexPrintf("Problem build\n");
    }
    
    
    // Set up the residual function (Cost Function)
    CostFunction* cost_function = new D2DCostFunction(n,p, prhs);
    // CostFunction* cost_function = new QuadraticCostFunction;     // Test CostFunction for simple problem
    if(printLevel > 0)
    {
        mexPrintf("CostFunction created\n"); 
    }
    
    LossFunction* loss = NULL;
    // Add loss function to problem
    switch(LossFunctionType)
    {
        case 1:
            break;
        case 2:
            loss = new TrivialLoss();
            break;
        case 3:
            loss = new HuberLoss(LossFunctionVar);
            break;
        case 4:
            loss = new SoftLOneLoss(LossFunctionVar);
            break; 
        case 5:
            loss = new CauchyLoss(LossFunctionVar);
            break; 
        case 6:
            loss = new ArctanLoss(LossFunctionVar);
            break;  
        default:
            mexErrMsgTxt("Error: invalid choice for LossFunction. Please choose int from 1 to 6"); 
            break;       
    }
    
 
    // Add Residual function to problem
    problem.AddResidualBlock(cost_function, loss, x);

    // Set bounds on parameters
    for ( i=0; i < p; i++ )
    {
        problem.SetParameterLowerBound(x,i,lb[i]);
        problem.SetParameterUpperBound(x,i,ub[i]);
    }
    if(printLevel > 0)
    {
        mexPrintf("CostFunction added\n"); 
    }
    
    
    
        
    // Initalize Options struct for solver
    Solver::Options options;
    
    /////////////////////////////////////////////////////////////////////////////////
    //////////////// OPTIONS READOUT ////////////////////////////////////////////////
    
    
    // Trust-region strategy
    // Options are 
    // 1: Dogleg
    // 2: Levenberg-Marquardt   
    if(TrustRegionStrategyType == 1)
    {
        options.trust_region_strategy_type = ceres::DOGLEG;
    }
    else if(TrustRegionStrategyType == 2)
    {
        options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    }
    
    // Dogleg Type
    // 1: Traditional Dogleg
    // 2: Subspace Dogleg      
    if(DoglegType == 1)
    {
        options.dogleg_type = ceres::SUBSPACE_DOGLEG;
    }
    else if(DoglegType == 2)
    {
        options.dogleg_type = ceres::TRADITIONAL_DOGLEG;
    }
    else
    {
        mexPrintf("Dogleg Type Setup invalid. Standard is used!\n");
    };
    
       
    
    // Maximum number of iterations
    options.max_num_iterations = MaxIter;  

    // Abortion Tolerances
    // to function change with each step, Default: 1e-6
    options.function_tolerance = TolFun;
    // to gradient change , Default: 1e-10
    options.gradient_tolerance = TolGradient;
    // to parameter change , Default: 1e-8
    options.parameter_tolerance = TolX;     
    
    

    // "Relax the requirement that the trust-region algorithm take strictly decreasing steps." Default: false
    options.use_nonmonotonic_steps = useNonmonotonicSteps;
    // The window size used by the step selection algorithm to accept non-monotonic steps. Default: 5
    options.max_consecutive_nonmonotonic_steps = maxConsecutiveNonmonotonicSteps;
    
    
    // Maximum time solver runs in seconds
    options.max_solver_time_in_seconds = maxSolverTimeInSeconds;

    // Number of threads to evaluate Jacobian, Default: 1
    options.num_threads = NumThreads;

    // Number of threads used by lin solver, Default: 1
    options.num_linear_solver_threads = NumLinearSolverThreads;

    // Inital trus-region radius, Default: 1e4
    options.initial_trust_region_radius = InitialTrustRegionRadius;
    
   

    // Maximum/Minimum trust region radius, Default: 1e16, 1e-3
    options.max_trust_region_radius = MaxTrustRegionRadius;  
    options.min_trust_region_radius = MinTrustRegionRadius;

    // Minimum relative decrease before trust-region step is accepted
    options.min_relative_decrease = MinRelativeDecrease;

    // Minimal/Maximal values for diagonal regulazation matrix of Levenber-Marquardt method, Default: 1e6, 1e32
    options.min_lm_diagonal = MinLMDiagonal;
    options.max_lm_diagonal = MaxLMDiagonal;

    // Maximum number of retries (with smaller trust-region/different conditioning) for solver if given trust-region leads to invalid results, Default: 5
    options.max_num_consecutive_invalid_steps = MaxNumConsecutiveInvalidSteps;


    // Linear solver type to solve linear least square problem, 
    // Options:
    // 1  DENSE_NORMAL_CHOLESKY  !!Dependency on EIGEN!!  
    // 2  DENSE_QR (for small problems)
    // 3  SPARSE_NORMAL_CHOLESKY !!Dependency on SuiteSparse and EIGEN!!
    // 4  CGNR for general sparse problems, !!inexact step algorithm used!!
    // for bundle adjustment problems:
    // 5  DENSE_SCHUR
    // 6  SPARSE_SCHUR
    // 7  ITERATIVE_SCHUR
    switch(LinearSolverType)
    {
        case 1:
            options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
            break;
        case 2:
            options.linear_solver_type = ceres::DENSE_QR;
            break;
        case 3:
            options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
            break;
        case 4:
            options.linear_solver_type = ceres::CGNR;
            break;
        case 5:
             options.linear_solver_type = ceres::DENSE_SCHUR;
            break;
        case 6:
            options.linear_solver_type = ceres::SPARSE_SCHUR;
            break;
        case 7:
            options.linear_solver_type = ceres::ITERATIVE_SCHUR;
            break;
    }



    // Determines weather Jacobian is scaled by norm of its columns before being passed to solver,
    // improves numerical conditioning of equations, Default: true
    options.jacobi_scaling = JacobiScaling;


    // Inner iterations (extra optimization on each trust-region step), Default: false, 1e-3
    options.use_inner_iterations = useInnerIterations;
    options.inner_iteration_tolerance = InnerIterationTolerance;


    printLevel = 0;
    if(printLevel >0)
    {
       // If set to true, logging output is sent to STDOUT, Default: false
      options.minimizer_progress_to_stdout = true;

      // Logging type of solver
      // Options:
      // SILENT
      // PER_MINIMIZER_ITERATION
      options.logging_type = ceres::PER_MINIMIZER_ITERATION;
    }
    if(printLevel = 0)
    {
      options.minimizer_progress_to_stdout = false;
      options.logging_type = ceres::SILENT;
    }

  /////////////////////// END OPTIONS READOUT //////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
    
    // Summary of solver output
    Solver::Summary summary;
    if(printLevel > 0)
    {
        mexPrintf("Solver initialized\n");
    }
    
    Solve(options, &problem, &summary);
    if(printLevel > 0)
    {
        mexPrintf("Solver done!\n"); 
    }
       
    
    // Print summary of ceres
    if(printLevel > 0)
    {
        cout << summary.FullReport() << endl; 
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    //////////// END OF CERES PART OF PROGRAM ////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    
    // Handle output of mex File
        
    // Resnorm
    if ( nlhs > 1 ) {
        plhs[1]     = mxCreateDoubleMatrix(1,1, mxREAL);
        finCost        = mxGetPr( plhs[1] );
        *finCost       = summary.final_cost;             
    }   
    
    
    // Residual and Jacobian (if requested)
    if ( ( nlhs > 2 )) 
    {
        UData *temp = new UData( p, prhs );
        
        int numrhsarguments=2;
        mxArray *lhs[2];            // Temporary varibale for output   

        memcpy( mxGetPr(temp->getParameterArray()), x, p*sizeof(double) );        

        mexCallMATLAB( 2, lhs, temp->getNumberofArguments(), temp->getCallData(), "feval" );

        // Point memory to the residual
        plhs[2] = lhs[0];

        // If Jacobian requested, point memory to it
        if((nlhs > 5))      
        {
                 plhs[5] = lhs[1];
        }    
        
        delete temp;
    }    
    
    
    
    // Exitflag
      // Iterations
    if ( nlhs > 3 ) {
        plhs[3]     = mxCreateDoubleMatrix(1,1, mxREAL);
        texitflag        = mxGetPr( plhs[3] );
        *texitflag       = summary.IsSolutionUsable();
    }      
    
    // Iterations
    if ( nlhs > 4 ) {
        plhs[4]     = mxCreateDoubleMatrix(1,1, mxREAL);
        iter        = mxGetPr( plhs[4] );
        *iter       = summary.num_successful_steps+summary.num_unsuccessful_steps;
    }        
}


double getValueFromStruct( const mxArray *prhs[], const char* fieldName, double oldValue )
{
    double* value;
    bool*   valuebool;
    int fieldNumber;
    mxArray *field;
    
    fieldNumber = mxGetFieldNumber(parOPTS, fieldName);
    
    /* Does the field exist? */
    if ( fieldNumber != -1 )
    {
        field = mxGetFieldByNumber(parOPTS, 0, fieldNumber);
        /* Not a valid value? Return -1.0f */
        if(mxIsLogical(field))
        {
            // Fetch logical values
            valuebool = mxGetLogicals(field);
            return valuebool[0];
        }
        else if( ( field == NULL ) || !mxIsDouble(field) || mxIsComplex(field) || mxIsEmpty(field) )
        {
            if (!mxIsEmpty(field))
                mexPrintf( "Invalid value for field %s\n", fieldName );
            return oldValue;
        } 
        else
        {
            /* Fetch the value */
            value = mxGetPr(field);
            return value[0];
        }
    } else
    {
        /* Not a valid value? Return -1.0f */
        return oldValue;
    }
    
    
    
}



// Load options from (custom) ceres options struct
void loadOptions( 
     const mxArray *prhs[], 
     int*        TrustRegionStrategyType,
     int*        DoglegType,
     int*        LossFunctionType, 
     double*     LossFunctionVar,
     double*     TolFun,
     double*     TolX,
     double*     TolGradient,
     int*        MaxIter,
     bool*       useNonmonotonicSteps,
     int*        maxConsecutiveNonmonotonicSteps,
     double*     maxSolverTimeInSeconds,        
     int*        NumThreads,
     int*        NumLinearSolverThreads,       
     double*     InitialTrustRegionRadius,
     double*     MaxTrustRegionRadius,
     double*     MinTrustRegionRadius,
     double*     MinRelativeDecrease,
     double*     MinLMDiagonal,
     double*     MaxLMDiagonal,
     int*        MaxNumConsecutiveInvalidSteps,
     bool*       JacobiScaling,
     bool*       useInnerIterations,
     double*     InnerIterationTolerance,
     int*        LinearSolverType,
     int*        printLevel)
{
//    mexPrintf("Loading Options...\n");
    char* stringBuffer;
    int fieldNumber;
    mxArray *field;
    
    
     *TrustRegionStrategyType           = (int) getValueFromStruct( prhs, "TrustRegionStrategyType", (double) *TrustRegionStrategyType );
     *DoglegType                        = (int) getValueFromStruct( prhs, "DoglegType", (double) *DoglegType );
     *LossFunctionType                  = (int) getValueFromStruct( prhs, "LossFunctionType", (double) *LossFunctionType );
     *LossFunctionVar                   = getValueFromStruct( prhs, "LossFunctionVar", (double) *LossFunctionVar );
     *TolFun                            = getValueFromStruct( prhs, "TolFun", (double) *TolFun );
     *TolX                              = getValueFromStruct( prhs, "TolX", (double) *TolX );
     *TolGradient                       = getValueFromStruct( prhs, "TolGradient", (double) *TolGradient );
     *MaxIter                           = (int) getValueFromStruct( prhs, "MaxIter", (double) *MaxIter );
     *useNonmonotonicSteps              = (bool) getValueFromStruct( prhs, "useNonmonotonicSteps", (double) *useNonmonotonicSteps );
     *maxConsecutiveNonmonotonicSteps   = (int) getValueFromStruct( prhs, "maxConsecutiveNonmonotonicSteps", (double) *maxConsecutiveNonmonotonicSteps );
     *maxSolverTimeInSeconds            = getValueFromStruct( prhs, "maxSolverTimeInSeconds", (double) *maxSolverTimeInSeconds );
     *NumThreads                        = (int) getValueFromStruct( prhs, "NumThreads", (double) *NumThreads );
     *NumLinearSolverThreads            = (int) getValueFromStruct( prhs, "NumLinearSolverThreads", (double) *NumLinearSolverThreads );
     *InitialTrustRegionRadius          = getValueFromStruct( prhs, "InitialTrustRegionRadius", (double) *InitialTrustRegionRadius );
     *MaxTrustRegionRadius              = getValueFromStruct( prhs, "MaxTrustRegionRadius", (double) *MaxTrustRegionRadius );
     *MinTrustRegionRadius              = getValueFromStruct( prhs, "MinTrustRegionRadius", (double) *MinTrustRegionRadius );
     *MinRelativeDecrease               = getValueFromStruct( prhs, "MinRelativeDecrease", (double) *MinRelativeDecrease );
     *MinLMDiagonal                     = getValueFromStruct( prhs, "MinLMDiagonal", (double) *MinLMDiagonal );
     *MaxLMDiagonal                     = getValueFromStruct( prhs, "MaxLMDiagonal", (double) *MaxLMDiagonal );
     *MaxNumConsecutiveInvalidSteps     = (int) getValueFromStruct( prhs, "MaxNumConsecutiveInvalidSteps", (double) *MaxNumConsecutiveInvalidSteps );
     *JacobiScaling                     = (bool) getValueFromStruct( prhs, "JacobiScaling", (double) *JacobiScaling );
     *useInnerIterations                = (bool) getValueFromStruct( prhs, "useInnerIterations", (double) *useInnerIterations );
     *InnerIterationTolerance           = getValueFromStruct( prhs, "InnerIterationTolerance", (double) *InnerIterationTolerance );
     *LinearSolverType                  = (int) getValueFromStruct( prhs, "LinearSolverType", (double) *LinearSolverType );
     *printLevel                        = (int) getValueFromStruct( prhs, "printLevel", (double) *printLevel );
}


// Call Function to check if it gives out Jacobian
void CallFunctionCheck(const mxArray *prhs[], int p, int* nresidals, int printLevel)
{
    if(printLevel > 0)
    {
        mexPrintf("Checking function output...\n");
    }
    
    UData *tempUData = new UData(p,prhs);

    int result;                 // Temporary result variable
    mxArray *lhs[2];            // Temporary varibale for output   
    
    double  *x0;                 	// Intial parameter values 
    
    // Fetch inital parameter values
    x0 = mxGetPr(parX0);

    // Insert initial parameters 
    memcpy( mxGetPr(tempUData->getParameterArray()), x0, p*sizeof(double) );        
    

    // Evaluate function to get length of residuals and see if it is sensitivity capable
    result = mexCallMATLAB( 2, lhs, tempUData->getNumberofArguments(), tempUData->getCallData(), "feval" );
    if ( result )
        mexErrMsgTxt("Error calling objective function");    
    
    // Fetch residual (to check the size)
    if(!mxIsDouble(lhs[0]) || mxIsComplex(lhs[0]) || mxIsEmpty(lhs[0]))
        mexErrMsgTxt("Error: Objective fumexCallMATLABnction did not return residual vector");
    else
        *nresidals = (int) mxGetNumberOfElements( lhs[0] );
    


    // Check existence and size of Jacobian
    if(!mxIsDouble(lhs[1]) || mxIsComplex(lhs[1]) || mxIsEmpty(lhs[1]))
    {
        mexErrMsgTxt("Error: Objective function did not return Jacobian.");
    } 
    else
    {
        // Check if Jacobian is the correct size before we segfault all over the place
        if ( !( (mxGetN( lhs[1] ) == p) && ( mxGetM( lhs[1] ) == *nresidals ) ) )
        {
            mexErrMsgTxt("Error: Jacobian wrong size.");
        }
    }        
    mexPrintf("Function output check successful\n");
    
    // Clean up
    delete tempUData;
}




// Check Input arguments
// Similar format as lsqnonlin:    lsqnonlin := fit(FUN,X0,LB,UB,OPTIONS)
//                                 ceresd2d  := fit(FUN,X0,LB,UB,OPTIONS)
void validateInput( const mxArray *prhs[], int nrhs, int *npars, int *bounded )
{
//    mexPrintf("Validating Input...\n");
    int pars;
    
    *bounded = 0;
    *npars = 0;
    
    if(nrhs < 2)
    {
        mexPrintf("Ceres Solver 1.11 (by Sameer Agarwal and Keir Mierle and Others)\n");
        mexPrintf("MATLAB wrapper by Franz-Georg Wieland (contact: franz-georg.wieland at mars.uni-freiburg.de)\n");
        mexErrMsgTxt("You must supply at least 2 arguments to ceres!\n\n[X,RESNORM,RESIDUAL,EXITFLAG,ITERATIONS,FEVALS,JACOBIAN] = mexnl2sol(fun,x0,(lb),(ub),(opts),(printlevel))\n");
    }
        
    /* Check Types */
    if(!mxIsFunctionHandle(parFUN) && !mxIsChar(parFUN))
        mexErrMsgTxt("First argument must be a function handle that computes the RHS!");
    if(!mxIsDouble(parX0) || mxIsComplex(parX0) || mxIsEmpty(parX0))
        mexErrMsgTxt("Second argument must be a real vector of initial parameters");

    /* Get size of the parameter vector */
    pars = (int) mxGetNumberOfElements(parX0);
    *npars = pars;
    
    /* Check whether we have bounds */
    if ( nrhs == 3 )
        mexErrMsgTxt("Specify both lower and upper bounds or none at all.");
    
    if ( nrhs > 3 )
    {
        /* If both lb and ub are empty ==> no bounds */
        if( mxIsEmpty(parLB) && mxIsEmpty(parUB) )
            *bounded = 0;
        else {
            if (mxIsEmpty(parLB) || mxIsEmpty(parUB))
                mexErrMsgTxt("Specify both lower and upper bounds or none at all.");
 
            if( !mxIsDouble(parLB) || mxIsComplex(parLB) || ( pars != (int) mxGetNumberOfElements( parLB ) ) )
                mexErrMsgTxt("Lower bound must be specified as a numeric vector of the same length as p0");
    
            if( !mxIsDouble(parUB) || mxIsComplex(parUB) || ( pars != (int) mxGetNumberOfElements( parUB ) ) )
                mexErrMsgTxt("Upper bound must be specified as a numeric vector of the same length as p0");
            
            *bounded = 1;
        }
                
        if ( nrhs > 4 )
            if ( !mxIsStruct( parOPTS ) )
                mexErrMsgTxt("Fifth argument (when specified) must be an OPTIONS struct");
        
        if ( nrhs > 5 )
            if ( mxIsEmpty( parPRINTLEVEL ) || !mxIsDouble( parPRINTLEVEL ) || mxIsComplex( parPRINTLEVEL ) )
                mexErrMsgTxt("Last argument (when specified) must be a numeric value indicating a print level");        
    }
}

void PrintFunctionInformation(int MaxIter, int TolFun, int TolX, int TolGradient, int bounded, int n, int p)
{
    mexPrintf("Ceres v1.11 - Settings\n");
    mexPrintf("  Maximum iterations:           %d\n", MaxIter );
    mexPrintf("  Absolute function tolerance:  %e\n", TolFun );
    mexPrintf("  X tolerance:                  %e\n", TolX );    
    mexPrintf("  Gradient tolerance:           %e\n", TolGradient );
    if ( bounded == 1 )
        mexPrintf("  Bounds:                       Enabled\n" );
    else
        mexPrintf("  Bounds:                       Not specified\n" );    
    mexPrintf("  Num. parameters               %d\n", p );
    mexPrintf("  Num. residuals                %d\n\n", n );
}  



#include"math.h"
# include "mex.h"
#include <omp.h>
#include "C:\Users\liuzhenye\Desktop\test2\ÄâÅ£¶Ù·¨\2\Eigen\Eigen\Dense"
//#include <Eigen/Dense>
#include "C:\Users\liuzhenye\Desktop\test2\ÄâÅ£¶Ù·¨\2\Eigen\Eigen\Eigen"

using namespace std;
using namespace Eigen;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize actual_number_of_non_zeros;
	double *ptr;
	double *value;
#if MX_HAS_INTERLEAVED_COMPLEX
	mxComplexDouble *complexPtr;
#endif

	/* Check for proper number of input and output arguments */
	if (nrhs != 2 && nrhs != 3)
	{
		mexErrMsgIdAndTxt("MATLAB:updatevalue:invalidNumInputs",
			"Two or three input argument required.");
	}

	if (!mxIsSparse(prhs[0]))
	{
		mexErrMsgIdAndTxt("MATLAB:updatevalue:invalidInputSparisty",
			"Input argument must be sparse\n");
	}

	//mexPrintf("%d\n", mxGetNumberOfElements(prhs[1]));
	//mexErrMsgTxt("Three input required.");

	if (nrhs == 2)
	{
		ptr = mxGetPr(prhs[0]);
		value = mxGetPr(prhs[1]);
		memcpy(ptr, value, mxGetNumberOfElements(prhs[1]) * sizeof(double));
	/*	for (int i = 0; i < mxGetNumberOfElements(prhs[1]); i++)
		{
			ptr[i] = value[i];
		}*/
	}
	else if(nrhs == 3)
	{
		double *index;
		ptr = mxGetPr(prhs[0]);
		value = mxGetPr(prhs[1]);
		index = mxGetPr(prhs[2]);
		#pragma omp parallel for
	/*	mexPrintf("%d\n", mxGetNumberOfElements(prhs[2]));*/
		for (int i = 0; i < mxGetNumberOfElements(prhs[2]); i++)
		{
	///*		mxprintf("%ld %lf\n", index[i], index[i]);*/			
	//	/*	if (index[i] >= mxGetNzmax(prhs[0]))
	//		{
	//			break;
	//		}*/
			//mexPrintf("%d\n", (int)index[i]);
			ptr[i] = value[(int)index[i]-1];
		}
	}

}

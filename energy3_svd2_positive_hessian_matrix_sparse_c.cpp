#include"math.h"
# include "mex.h"
#include <omp.h>
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Dense"
//#include <Eigen/Dense>
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Eigen"

void hessian_matrix_triplet(const int face_number, const int point_number, const double * points, const double * faces, const double * l_target, double * col, double *row, double *value)
{
	int i, j, k, l;
	#pragma omp parallel for
	for (int i = 0; i < face_number; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int xi_index = faces[j*face_number + i];
			int yi_index = xi_index + point_number;
			int zi_index = xi_index + 2 * point_number;
			int xj_index = faces[((j + 1) % 3)*face_number + i];
			int yj_index = faces[((j + 1) % 3)*face_number + i] + point_number;
			int zj_index = faces[((j + 1) % 3)*face_number + i] + 2 * point_number;

			/*	int xi_index = faces[i*3+j];
				int yi_index = xi_index + 1;
				int zi_index = xi_index + 2 ;
				int xj_index = faces[((j + 1) % 3) + i*3];
				int yj_index = xj_index + 1;
				int zj_index = zj_index + 2 ;*/

			double xi = points[xi_index - 1];
			double yi = points[yi_index - 1];
			double zi = points[zi_index - 1];
			double xj = points[xj_index - 1];
			double yj = points[yj_index - 1];
			double zj = points[zj_index - 1];
			double l = l_target[j*face_number + i];
			/*	row[(i-1)*108+(j-1)*36] = xi_index;
				row[i * 108 + j * 36+1] = xi_index;
				row[i * 108 + j * 36+2] = xi_index;
				row[i * 108 + j * 36+3] = xi_index;
				row[i * 108 + j * 36+4] = xi_index;
				row[i * 108 + j * 36+5] = xi_index;*/
			double ld=4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) - 4 * l*l;
			double ld2= 12 * pow((xi - xj), 2.0) + 12 * pow((yi - yj), 2.0) + 12 * pow((zi - zj), 2.0) - 4 * l*l;
			double norm = sqrt(pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0));
			double x = (xi - xj) / norm;
			double y = (yi - yj) / norm;
			double z = (zi - zj) / norm;
			double lambda = (norm - l) / norm * 2;

			row[i * 108 + j * 36] = xi_index; col[i * 108 + j * 36] = xi_index; value[i * 108 + j * 36] = 4*(xi-xj)*(xi-xj) / (2 * (norm*norm)) - (2 * (l - norm)) / norm + (4*(xi-xj)*(xi-xj) * (l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*(1 - x *x);
			row[i * 108 + j * 36 + 1] = yi_index; col[i * 108 + j * 36 + 1] = yi_index; value[i * 108 + j * 36 + 1] = 4*(yi-yj)*(yi-yj) / (2 * (norm*norm)) - (2 * (l - norm)) / norm + (4*(yi-yj)*(yi-yj) * (l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*(1 - y *y);
			row[i * 108 + j * 36 + 2] = zi_index; col[i * 108 + j * 36 + 2] = zi_index; value[i * 108 + j * 36 + 2] = 4*(zi-zj)*(zi-zj) / (2 * (norm*norm)) - (2 * (l - norm)) / norm + (4*(zi-zj)*(zi-zj) * (l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*(1 - z*z);
			row[i * 108 + j * 36 + 3] = yi_index; col[i * 108 + j * 36 + 3] = xi_index; value[i * 108 + j * 36 + 3] = ((2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (2 * (norm*norm)) + ((2 * xi - 2 * xj)*(2 * yi - 2 * yj)*(l - norm)) / (2 * norm*norm*norm) - (lambda < 0)*abs(lambda)*x*y;
			row[i * 108 + j * 36 + 4] = zi_index; col[i * 108 + j * 36 + 4] = xi_index; value[i * 108 + j * 36 + 4] = ((2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) + ((2 * xi - 2 * xj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) - (lambda < 0)*abs(lambda)*x*z;
			row[i * 108 + j * 36 + 5] = zi_index; col[i * 108 + j * 36 + 5] = yi_index; value[i * 108 + j * 36 + 5] = ((2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) + ((2 * yi - 2 * yj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) - (lambda < 0)*abs(lambda)*z*y;
			memcpy(value + i * 108 + j * 36 + 6, value + i * 108 + j * 36+3 , 3 * sizeof(double));
			row[i * 108 + j * 36 + 6] = xi_index; col[i * 108 + j * 36 + 6] = yi_index; //value[i * 108 + j * 36 + 6] = -2 * x * y*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
			row[i * 108 + j * 36 + 7] = xi_index; col[i * 108 + j * 36 + 7] = zi_index; //value[i * 108 + j * 36 + 9] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
			row[i * 108 + j * 36 + 8] = yi_index; col[i * 108 + j * 36 + 8] = zi_index;// value[i * 108 + j * 36 + 10] = -2 * z * y*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
			memcpy(value + i * 108 + j * 36 + 9, value + i * 108 + j * 36 , 9 * sizeof(double));
			row[i * 108 + j * 36 + 9] = xj_index; col[i * 108 + j * 36 + 9] = xj_index;//value[i * 108 + j * 36 + 1] = ld + 2*(ld < 0)*abs(ld)*(y*y + z*z) + 2 * pow(2 * xi - 2 * xj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*x*x;			
			row[i * 108 + j * 36 + 10] = yj_index; col[i * 108 + j * 36 + 10] = yj_index;//value[i * 108 + j * 36 + 3] = ld + 2 * (ld < 0)*abs(ld)*(x*x + z*z) + 2 * pow(2 * yi - 2 * yj, 2.0)+2 * (ld2 < 0)*abs(ld2)*y*y;	
			row[i * 108 + j * 36 + 11] = zj_index; col[i * 108 + j * 36 + 11] = zj_index;// value[i * 108 + j * 36 + 5] = ld + 2 * (ld < 0)*abs(ld)*(y*y + x*x) + 2 * pow(2 * zi - 2 * zj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*z *z;
			row[i * 108 + j * 36 + 12] = yj_index; col[i * 108 + j * 36 + 12] = xj_index; //value[i * 108 + j * 36 + 31] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
			row[i * 108 + j * 36 + 13] = zj_index; col[i * 108 + j * 36 + 13] = xj_index; //value[i * 108 + j * 36 + 33] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
			row[i * 108 + j * 36 + 14] = zj_index; col[i * 108 + j * 36 + 14] = yj_index; //value[i * 108 + j * 36 + 34] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
			row[i * 108 + j * 36 + 15] = xj_index; col[i * 108 + j * 36 + 15] = yj_index; //value[i * 108 + j * 36 + 30] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
			row[i * 108 + j * 36 + 16] = xj_index; col[i * 108 + j * 36 + 16] = zj_index; //value[i * 108 + j * 36 + 32] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;			
			row[i * 108 + j * 36 + 17] = yj_index; col[i * 108 + j * 36 + 17] = zj_index; //value[i * 108 + j * 36 + 35] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;

			row[i * 108 + j * 36 + 18] = xi_index; col[i * 108 + j * 36 + 18] = xj_index; value[i * 108 + j * 36 + 18] = (2 * (l - norm)) / norm - 4*(xi-xj)*(xi-xj) / (2 * (norm*norm)) - (4*(xi-xj)*(xi-xj) * (l - norm)) / (2 * norm*norm*norm)- (lambda < 0)*abs(lambda)*(1 - x * x);
			row[i * 108 + j * 36 + 19] = yi_index; col[i * 108 + j * 36 + 19] = yj_index; value[i * 108 + j * 36 + 19] = (2 * (l - norm)) / norm - 4*(yi-yj)*(yi-yj) / (2 * (norm*norm)) - (4*(yi-yj)*(yi-yj) * (l - norm)) / (2 * norm*norm*norm)- (lambda < 0)*abs(lambda)*(1 - y * y);
			row[i * 108 + j * 36 + 20] = zi_index; col[i * 108 + j * 36 + 20] = zj_index; value[i * 108 + j * 36 + 20] = (2 * (l - norm)) / norm - 4*(zi-zj)*(zi-zj) / (2 * (norm*norm)) - (4*(zi-zj)*(zi-zj) * (l - norm)) / (2 * norm*norm*norm)- (lambda < 0)*abs(lambda)*(1 - z * z);
			row[i * 108 + j * 36 + 21] = yi_index; col[i * 108 + j * 36 + 21] = xj_index; value[i * 108 + j * 36 + 21] = -((2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (2 * (norm*norm)) - ((2 * xi - 2 * xj)*(2 * yi - 2 * yj)*(l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*x*y;
			row[i * 108 + j * 36 + 22] = zi_index; col[i * 108 + j * 36 + 22] = xj_index; value[i * 108 + j * 36 + 22] = -((2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) - ((2 * xi - 2 * xj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*x*z;
			row[i * 108 + j * 36 + 23] = zi_index; col[i * 108 + j * 36 + 23] = yj_index; value[i * 108 + j * 36 + 23] = -((2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) - ((2 * yi - 2 * yj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*z*y;
			memcpy(value + i * 108 + j * 36 + 24, value + i * 108 + j * 36 + 21, 3 * sizeof(double));
			row[i * 108 + j * 36 + 24] = xi_index; col[i * 108 + j * 36 + 24] = yj_index; //value[i * 108 + j * 36 + 18] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 25] = xi_index; col[i * 108 + j * 36 + 25] = zj_index; //value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 26] = yi_index; col[i * 108 + j * 36 + 26] = zj_index; //value[i * 108 + j * 36 + 24] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
			memcpy(value + i * 108 + j * 36 + 27, value + i * 108 + j * 36+18, 9 * sizeof(double));
			row[i * 108 + j * 36 + 27] = xj_index; col[i * 108 + j * 36 + 27] = xi_index; //value[i * 108 + j * 36 + 13] = -ld - 2 * pow(2 * xi - 2 * xj, 2);			
			row[i * 108 + j * 36 + 28] = yj_index; col[i * 108 + j * 36 + 28] = yi_index; //value[i * 108 + j * 36 + 15] = -ld - 2 * pow(2 * yi - 2 * yj, 2);			
			row[i * 108 + j * 36 + 29] = zj_index; col[i * 108 + j * 36 + 29] = zi_index; //value[i * 108 + j * 36 + 17] = -ld - 2 * pow(2 * zi - 2 * zj, 2);
			row[i * 108 + j * 36 + 30] = yj_index; col[i * 108 + j * 36 + 30] = xi_index; //value[i * 108 + j * 36 + 19] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 31] = zj_index; col[i * 108 + j * 36 + 31] = xi_index; //value[i * 108 + j * 36 + 20] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 32] = zj_index; col[i * 108 + j * 36 + 32] = yi_index; //value[i * 108 + j * 36 + 25] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);											
			row[i * 108 + j * 36 + 33] = xj_index; col[i * 108 + j * 36 + 33] = yi_index; //value[i * 108 + j * 36 + 23] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);						
			row[i * 108 + j * 36 + 34] = xj_index; col[i * 108 + j * 36 + 34] = zi_index; //value[i * 108 + j * 36 + 26] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);			
			row[i * 108 + j * 36 + 35] = yj_index; col[i * 108 + j * 36 + 35] = zi_index; //value[i * 108 + j * 36 + 28] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
			

		
	/*		if (ld2 < 0)
			{
				for (int k = 0; k < 6; k++)
				{
					value[i * 108 + j * 36 + k] = 1e-19;
				}
				for (int k = 6; k < 36; k++)
				{
					value[i * 108 + j * 36 + k] = 1e-20;
				}
			}*/

				
		}
	}

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *points, *l_target, *value, *col, *row, *faces;
	//int ;
	bool transpose;
	int RowA, ColA, RowB, ColB;
	if (nlhs != 4) {
		mexErrMsgTxt("One output required.");
	}
	else if (nrhs != 3) {
		mexErrMsgTxt("Three input required.");
	}
	//if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[0]) || mxIsComplex(prhs[1])) {
	//	mexErrMsgTxt("Input Array must be double.");
	//}
	//获得矩阵的行数
	size_t point_number = mxGetM(prhs[0]);
	//if (point_number != 3180)
	//{
	//	mexErrMsgTxt("point_number");
	//}

	//获得矩阵的列数
	ColA = mxGetN(prhs[0]);
	int face_number = mxGetM(prhs[1]);
	ColB = mxGetN(prhs[1]);
	//判断行列是否相等
	//if (RowA != RowB || ColA != ColB) {
	//	mexErrMsgTxt("Rows and Cols must be same.");
	//}
	//获取指向输入参数的指针
	points = mxGetPr(prhs[0]);
	faces = mxGetPr(prhs[1]);
	l_target = mxGetPr(prhs[2]);
	//生成输出参量的mxArray
	plhs[0] = mxCreateDoubleMatrix(face_number, 108, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(face_number, 108, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(face_number, 108, mxREAL);
	int cmplx;
	cmplx = 0;
	
	//获取指向输出参数的指针
	col = mxGetPr(plhs[0]);
	row = mxGetPr(plhs[1]);
	value = mxGetPr(plhs[2]);
	//H= mxGetPr(plhs[3]);


	
	//mexWarnMsgIdAndTxt("1");
	hessian_matrix_triplet(face_number, point_number, points, faces, l_target, col, row, value);
	/*mexErrMsgTxt("2");*/
	Eigen::SparseMatrix < double >  A(3*point_number, 3*point_number);
	std::vector < Eigen::Triplet < double > > triplets;
	for (size_t i = 0; i < 108 * face_number; i++)
	{
		triplets.emplace_back(col[i]-1, row[i]-1, value[i]);
	}
	/*mexErrMsgTxt("3");*/
	A.setFromTriplets(triplets.begin(), triplets.end());
	//mexErrMsgTxt("4");
	int *innerIndex = A.innerIndexPtr();
	double *valueptr = A.valuePtr();
	/*mexErrMsgTxt("5");*/

	/*mexErrMsgTxt("8");*/
	//jcs[point_number] = A.nonZeros();
	/*mexErrMsgTxt("9");*/
	plhs[3] = mxCreateSparse(3 * point_number, 3 * point_number, A.nonZeros(), mxREAL);
	double *pr, *pi, *si, *sr;
	mwIndex *irs, *jcs;
	sr = mxGetPr(plhs[3]);
	/*si = mxGetPi(plhs[3]);*/
	irs = mxGetIr(plhs[3]);
	jcs = mxGetJc(plhs[3]);
	for (int i = 0; i < A.nonZeros(); i++)
	{
		sr[i] = valueptr[i];
	}
	/*mexErrMsgTxt("6");*/
	for (int i = 0; i < A.nonZeros(); i++)
	{
		irs[i] = innerIndex[i];
	}
	//mexErrMsgTxt("7");
	int *OuterStarts = A.outerIndexPtr();
	for (int i = 0; i < 3*point_number + 1; i++)
	{
		jcs[i] = OuterStarts[i];
	}


	//mxSetIr(plhs[3], irs);
	//mxSetJc(plhs[3], jcs);
	//mxSetPr(plhs[3], sr);

}


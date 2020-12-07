#include"math.h"
# include "mex.h"
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Dense"
//#include <Eigen/Dense>
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Eigen"

void hessian_matrix_triplet(const int face_number, const int point_number, const double * points, const double * faces, const double * l_target, double * col, double *row, double *value)
{
	int i, j, k, l;
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
			row[i * 108 + j * 36] = xi_index; col[i * 108 + j * 36] = xi_index; value[i * 108 + j * 36] = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) + 2 * pow(2 * xi - 2 * xj, 2.0) - 4 * l*l;
			row[i * 108 + j * 36 + 1] = xj_index; col[i * 108 + j * 36 + 1] = xj_index; value[i * 108 + j * 36 + 1] = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) + 2 * pow(2 * xi - 2 * xj, 2.0) - 4 * l*l;
			row[i * 108 + j * 36 + 2] = yi_index; col[i * 108 + j * 36 + 2] = yi_index; value[i * 108 + j * 36 + 2] = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) + 2 * pow(2 * yi - 2 * yj, 2.0) - 4 * l*l;
			row[i * 108 + j * 36 + 3] = yj_index; col[i * 108 + j * 36 + 3] = yj_index; value[i * 108 + j * 36 + 3] = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) + 2 * pow(2 * yi - 2 * yj, 2.0) - 4 * l*l;
			row[i * 108 + j * 36 + 4] = zi_index; col[i * 108 + j * 36 + 4] = zi_index; value[i * 108 + j * 36 + 4] = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) + 2 * pow(2 * zi - 2 * zj, 2.0) - 4 * l*l;
			row[i * 108 + j * 36 + 5] = zj_index; col[i * 108 + j * 36 + 5] = zj_index; value[i * 108 + j * 36 + 5] = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) + 2 * pow(2 * zi - 2 * zj, 2.0) - 4 * l*l;

			row[i * 108 + j * 36 + 6] = xi_index; col[i * 108 + j * 36 + 6] = yi_index; value[i * 108 + j * 36 + 6] = 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 7] = yi_index; col[i * 108 + j * 36 + 7] = xi_index; value[i * 108 + j * 36 + 7] = 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 8] = zi_index; col[i * 108 + j * 36 + 8] = xi_index; value[i * 108 + j * 36 + 8] = 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 9] = xi_index; col[i * 108 + j * 36 + 9] = zi_index; value[i * 108 + j * 36 + 9] = 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 10] = yi_index; col[i * 108 + j * 36 + 10] = zi_index; value[i * 108 + j * 36 + 10] = 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 11] = zi_index; col[i * 108 + j * 36 + 11] = yi_index; value[i * 108 + j * 36 + 11] = 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);

			row[i * 108 + j * 36 + 12] = xi_index; col[i * 108 + j * 36 + 12] = xj_index; value[i * 108 + j * 36 + 12] = 4 * l*l - 4 * pow(yi - yj, 2.0) - 4 * pow(zi - zj, 2.0) - 2 * pow(2 * xi - 2 * xj, 2.0) - 4 * pow(xi - xj, 2.0);
			row[i * 108 + j * 36 + 13] = xj_index; col[i * 108 + j * 36 + 13] = xi_index; value[i * 108 + j * 36 + 13] = 4 * l*l - 4 * pow(yi - yj, 2) - 4 * pow(zi - zj, 2) - 2 * pow(2 * xi - 2 * xj, 2) - 4 * pow(xi - xj, 2);
			row[i * 108 + j * 36 + 14] = yi_index; col[i * 108 + j * 36 + 14] = yj_index; value[i * 108 + j * 36 + 14] = 4 * l*l - 4 * pow(yi - yj, 2) - 4 * pow(zi - zj, 2) - 2 * pow(2 * yi - 2 * yj, 2) - 4 * pow(xi - xj, 2);
			row[i * 108 + j * 36 + 15] = yj_index; col[i * 108 + j * 36 + 15] = yi_index; value[i * 108 + j * 36 + 15] = 4 * l*l - 4 * pow(yi - yj, 2) - 4 * pow(zi - zj, 2) - 2 * pow(2 * yi - 2 * yj, 2) - 4 * pow(xi - xj, 2);
			row[i * 108 + j * 36 + 16] = zi_index; col[i * 108 + j * 36 + 16] = zj_index; value[i * 108 + j * 36 + 16] = 4 * l*l - 4 * pow(yi - yj, 2) - 4 * pow(zi - zj, 2) - 2 * pow(2 * zi - 2 * zj, 2) - 4 * pow(xi - xj, 2);
			row[i * 108 + j * 36 + 17] = zj_index; col[i * 108 + j * 36 + 17] = zi_index; value[i * 108 + j * 36 + 17] = 4 * l*l - 4 * pow(yi - yj, 2) - 4 * pow(zi - zj, 2) - 2 * pow(2 * zi - 2 * zj, 2) - 4 * pow(xi - xj, 2);

			row[i * 108 + j * 36 + 18] = xi_index; col[i * 108 + j * 36 + 18] = yj_index; value[i * 108 + j * 36 + 18] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 19] = yj_index; col[i * 108 + j * 36 + 19] = xi_index; value[i * 108 + j * 36 + 19] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 20] = zj_index; col[i * 108 + j * 36 + 20] = xi_index; value[i * 108 + j * 36 + 20] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 21] = xi_index; col[i * 108 + j * 36 + 21] = zj_index; value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);

			row[i * 108 + j * 36 + 22] = yi_index; col[i * 108 + j * 36 + 22] = xj_index; value[i * 108 + j * 36 + 22] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 23] = xj_index; col[i * 108 + j * 36 + 23] = yi_index; value[i * 108 + j * 36 + 23] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 24] = yi_index; col[i * 108 + j * 36 + 24] = zj_index; value[i * 108 + j * 36 + 24] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 25] = zj_index; col[i * 108 + j * 36 + 25] = yi_index; value[i * 108 + j * 36 + 25] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);

			row[i * 108 + j * 36 + 26] = xj_index; col[i * 108 + j * 36 + 26] = zi_index; value[i * 108 + j * 36 + 26] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 27] = zi_index; col[i * 108 + j * 36 + 27] = xj_index; value[i * 108 + j * 36 + 27] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 28] = yj_index; col[i * 108 + j * 36 + 28] = zi_index; value[i * 108 + j * 36 + 28] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 29] = zi_index; col[i * 108 + j * 36 + 29] = yj_index; value[i * 108 + j * 36 + 29] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);

			row[i * 108 + j * 36 + 30] = xj_index; col[i * 108 + j * 36 + 30] = yj_index; value[i * 108 + j * 36 + 30] = 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 31] = yj_index; col[i * 108 + j * 36 + 31] = xj_index; value[i * 108 + j * 36 + 31] = 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			row[i * 108 + j * 36 + 32] = xj_index; col[i * 108 + j * 36 + 32] = zj_index; value[i * 108 + j * 36 + 32] = 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 33] = zj_index; col[i * 108 + j * 36 + 33] = xj_index; value[i * 108 + j * 36 + 33] = 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 34] = zj_index; col[i * 108 + j * 36 + 34] = yj_index; value[i * 108 + j * 36 + 34] = 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
			row[i * 108 + j * 36 + 35] = yj_index; col[i * 108 + j * 36 + 35] = zj_index; value[i * 108 + j * 36 + 35] = 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
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


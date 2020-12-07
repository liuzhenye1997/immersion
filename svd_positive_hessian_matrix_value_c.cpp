#include"math.h"
# include "mex.h"
#include <omp.h>
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Dense"
//#include <Eigen/Dense>
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Eigen"

using namespace std;
using namespace Eigen;

void hessian_matrix_value(const int face_number, const int point_number, const double * points, const double * faces, const double * l_target, double *value)
{
	int k, l;
	//omp_set_num_threads(16);
#pragma omp parallel for

	for (int i = 0; i < face_number; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			//mexPrintf("%d\n", i*3+j);
			int xi_index = faces[j*face_number + i];
			int yi_index = xi_index + point_number;
			int zi_index = xi_index + 2 * point_number;
			int xj_index = faces[((j + 1) % 3)*face_number + i];
			int yj_index = faces[((j + 1) % 3)*face_number + i] + point_number;
			int zj_index = faces[((j + 1) % 3)*face_number + i] + 2 * point_number;



			double xi = points[xi_index - 1];
			double yi = points[yi_index - 1];
			double zi = points[zi_index - 1];
			double xj = points[xj_index - 1];
			double yj = points[yj_index - 1];
			double zj = points[zj_index - 1];
			double l = l_target[j*face_number + i];

			double ld = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) - 4 * l*l;
			double ld2 = 12 * pow((xi - xj), 2.0) + 12 * pow((yi - yj), 2.0) + 12 * pow((zi - zj), 2.0) - 4 * l*l;
			double norm = sqrt(pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0));
			double x = (xi - xj) / norm;
			double y = (yi - yj) / norm;
			double z = (zi - zj) / norm;
			int ld0 = ld < 0;
			int ld20 = ld2 < 0;
			double abs_ld = abs(ld);
			double abs_ld2 = abs(ld2);
			double x2 = pow(2 * xi - 2 * xj, 2.0);
			double y2 = pow(2 * yi - 2 * yj, 2.0);
			double z2 = pow(2 * zi - 2 * zj, 2.0);
			double xy = (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			double xz = (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			double yz = (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
		
	
			value[i * 108 + j * 36] = ld + 2 * (ld < 0)*abs(ld)*(y*y + z * z) + 2 * pow(2 * xi - 2 * xj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*x *x;
			value[i * 108 + j * 36 + 1] = ld + 2 * (ld < 0)*abs(ld)*(x*x + z * z) + 2 * pow(2 * yi - 2 * yj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*y*y;
			value[i * 108 + j * 36 + 2] = ld + 2 * (ld < 0)*abs(ld)*(y*y + x * x) + 2 * pow(2 * zi - 2 * zj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*z *z;
			value[i * 108 + j * 36 + 3] = -2 * x * y*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
			value[i * 108 + j * 36 + 4] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
			value[i * 108 + j * 36 + 5] = -2 * z * y*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
			memcpy(value + i * 108 + j * 36 + 6, value + i * 108 + j * 36 + 3, 3 * sizeof(double));
			memcpy(value + i * 108 + j * 36 + 9, value + i * 108 + j * 36, 9 * sizeof(double));
			value[i * 108 + j * 36 + 18] = -ld - 2 * pow(2 * xi - 2 * xj, 2.0);
			value[i * 108 + j * 36 + 19] = -ld - 2 * pow(2 * yi - 2 * yj, 2);
			value[i * 108 + j * 36 + 20] = -ld - 2 * pow(2 * zi - 2 * zj, 2);
			value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
			value[i * 108 + j * 36 + 22] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
			value[i * 108 + j * 36 + 23] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
			memcpy(value + i * 108 + j * 36 + 24, value + i * 108 + j * 36 + 21, 3 * sizeof(double));
			memcpy(value + i * 108 + j * 36 + 27, value + i * 108 + j * 36 + 18, 9 * sizeof(double));

		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *points, *l_target, *value, *col, *row, *faces;
	//int ;
	bool transpose;
	int RowA, ColA, RowB, ColB;
	if (nlhs != 1) {
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

	plhs[0] = mxCreateDoubleMatrix(face_number * 108, 1, mxREAL);
	//获取指向输出参数的指针

	value = mxGetPr(plhs[0]);
	hessian_matrix_value(face_number, point_number, points, faces, l_target, value);

}


#include"math.h"
# include "mex.h"
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Dense"
//#include <Eigen/Dense>
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Eigen"

void best_step(const int face_number, const int point_number, const double * points, const double * faces, const double * l_target, double * d,double *coefficient)
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



			double xi = points[xi_index - 1];
			double yi = points[yi_index - 1];
			double zi = points[zi_index - 1];
			double xj = points[xj_index - 1];
			double yj = points[yj_index - 1];
			double zj = points[zj_index - 1];

			double dxi = d[xi_index - 1];
			double dxj = d[xj_index - 1];
			double dyi = d[yi_index - 1];
			double dyj = d[yj_index - 1];
			double dzi = d[zi_index - 1];
			double dzj = d[zj_index - 1];

			double l = l_target[j*face_number + i]* l_target[j*face_number + i];

			coefficient[0] = coefficient[0] + 2 * (pow(dxi - dxj, 2) + pow(dyi - dyj, 2) + pow(dzi - dzj, 2))*(2 * pow(dxi - dxj, 2) + 2 * pow(dyi - dyj, 2) + 2 * pow(dzi - dzj, 2));
			coefficient[1] = coefficient[1] + 2 * (pow(dxi - dxj , 2) + pow(dyi - dyj, 2) + pow(dzi - dzj, 2))*(2 * (dxi - dxj)*(xi - xj) + 2 * (dyi - dyj)*(yi - yj) + 2 * (dzi - dzj)*(zi - zj)) + 2 * (2 * (dxi - dxj)*(xi - xj) + 2 * (dyi - dyj)*(yi - yj) + 2 * (dzi - dzj)*(zi - zj))*(2 * pow(dxi - dxj, 2) + 2 * pow(dyi - dyj, 2) + 2 * pow(dzi - dzj, 2));
			coefficient[2] = coefficient[2] + (2 * pow(2 * (dxi - dxj)*(xi - xj) + 2 * (dyi - dyj)*(yi - yj) + 2 * (dzi - dzj)*(zi - zj),2) + 2 * (2 * pow(dxi - dxj, 2) + 2 * pow(dyi - dyj, 2) + 2 * pow(dzi - dzj, 2))*(pow(xi - xj, 2) - l + pow(yi - yj, 2) + pow(zi - zj, 2)));
			coefficient[3] = coefficient[3] + 2 * (2 * (dxi - dxj)*(xi - xj) + 2 * (dyi - dyj)*(yi - yj) + 2 * (dzi - dzj)*(zi - zj))*(pow(xi - xj, 2) - l + pow(yi - yj,2) + pow(zi - zj, 2));
		}
	}

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *points, *l_target, *value, *col, *row, *faces,*d;
	//int ;
	bool transpose;
	int RowA, ColA, RowB, ColB;
	if (nlhs != 1) {
		mexErrMsgTxt("One output required.");
	}
	else if (nrhs != 4) {
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
	d = mxGetPr(prhs[2]);
	l_target = mxGetPr(prhs[3]);
	//生成输出参量的mxArray
	plhs[0] = mxCreateDoubleMatrix(4, 1, mxREAL);
	double *coefficient;
	coefficient = mxGetPr(plhs[0]);
	
	best_step(face_number, point_number, points, faces, l_target,d, coefficient);
}


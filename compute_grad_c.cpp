#include"math.h"
# include "mex.h"
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Dense"
//#include <Eigen/Dense>
#include "C:\Users\liuzhenye\Desktop\test2\拟牛顿法\2\Eigen\Eigen\Eigen"

void compute_grad(const int face_number, const int point_number, const double * points_temp, const double * faces, const double * l_target, double * l_temp, double *grad)
{
	#pragma omp_set_num_threads(3)
	#pragma omp parallel for 
	for (int k = 0; k < 3; k++)
	{
		for (int i = 0; i < face_number; i++)
		{
			for (int j = 0; j < 3; j++)
			{

				int ij = faces[i + j * face_number] - 1;
				int ij1 = faces[i + (j + 1) % 3 * face_number] - 1;
				grad[ij + k * point_number] = grad[ij + k * point_number] + 2 * (pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2))*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]);
				grad[ij1 + k * point_number] = grad[ij1 + k * point_number] - 2 * (pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2))*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]);
			}
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *points, *l_target, *value, *col, *row, *faces,*d,*l_temp;
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
	l_target = mxGetPr(prhs[2]);
	l_temp = mxGetPr(prhs[3]);
	//生成输出参量的mxArray
	plhs[0] = mxCreateDoubleMatrix(point_number, 3, mxREAL);
	double *grad;
	grad = mxGetPr(plhs[0]);
	
	compute_grad(face_number, point_number, points, faces, l_target, l_temp,grad);
}


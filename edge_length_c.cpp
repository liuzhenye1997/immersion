# include "mex.h"
#include<cmath>
using namespace std;


//计算三角形网格边的长度
void edge_length(const int face_number, const int point_number, const double * points, const double * faces, double * length)
{
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

			double xi = points[xi_index - 1];
			double yi = points[yi_index - 1];
			double zi = points[zi_index - 1];
			double xj = points[xj_index - 1];
			double yj = points[yj_index - 1];
			double zj = points[zj_index - 1];
			length[j*face_number + i] = 0;
			length[j*face_number + i] += (xi - xj)*(xi - xj);
			length[j*face_number + i] += (yi - yj)*(yi - yj);
			length[j*face_number + i] += (zi - zj)*(zi - zj);
			length[j*face_number + i] = sqrt(length[j*face_number + i]);
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *points, *l_target, *value, *col, *row, *faces, *d;
	//int ;
	bool transpose;
	int RowA, ColA, RowB, ColB;
	if (nlhs != 1) {
		mexErrMsgTxt("One output required.");
	}
	else if (nrhs != 2) {
		mexErrMsgTxt("Two input required.");
	}

	//获得矩阵的行数
	size_t point_number = mxGetM(prhs[1]);


	//获得矩阵的列数
	ColA = mxGetN(prhs[0]);
	int face_number = mxGetM(prhs[0]);
	ColB = mxGetN(prhs[1]);

	//获取指向输入参数的指针
	points = mxGetPr(prhs[1]);
	faces = mxGetPr(prhs[0]);
	//生成输出参量的mxArray
	plhs[0] = mxCreateDoubleMatrix(face_number, 3, mxREAL);
	double *length;
	length = mxGetPr(plhs[0]);

	edge_length(face_number, point_number, points, faces, length);


}


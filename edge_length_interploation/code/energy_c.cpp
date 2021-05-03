
# include "mex.h"
#include "D:\test2\拟牛顿法\2\Eigen\Eigen\Dense"
//#include <Eigen/Dense>
#include "D:\test2\拟牛顿法\2\Eigen\Eigen\Eigen"
#include<cmath>
using namespace std;
using index_t = int64_t;
//using index_t = int;
using SpMat = Eigen::SparseMatrix<double, 0, index_t>;



double energy_a(double *faces,int face_number, int point_number, double *a, double *points_to_faces, double *length, double *point_to_points, int max_col, SpMat A_index)
{
	double E = 0;
	for (int i = 0; i < point_number; i++)
	{
		Eigen::Matrix3d M;
		M.setIdentity();

		for (int j = 0; j < points_to_faces[i]; j++)
		{
			double theta = 0;
			for (int k = 0; k < 3; k++)
			{
				if (faces[(int)points_to_faces[(j + 1)*point_number + i] - 1 + k * face_number] == (i + 1))
				{
					int m = points_to_faces[i + (j + 1)*point_number] - 1;
					theta = acos((length[m + k * face_number] * length[m + k * face_number] + length[m + ((k + 2) % 3) * face_number] * length[m + ((k + 2) % 3) * face_number] - length[m + ((k + 1) % 3) * face_number] * length[m + ((k + 1) % 3) * face_number]) / (2 * length[m + k * face_number] * length[m + ((k + 2) % 3) * face_number]));
					break;
				}
			}
		
			Eigen::Matrix3d Z;
			Z<< cos(theta), -sin(theta), 0,
				sin(theta), cos(theta), 0,
				0, 0, 1;
			M = M * Z;

			Eigen::Matrix3d X;
			theta = a[(int)A_index.coeffRef(i, (int)point_to_points[i + (((j + 1) % (int)points_to_faces[i]) + 1)*point_number] - 1)-1];
			X << 1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta);
			M = M * X;
		}
		Eigen::Matrix3d D = M - Eigen::MatrixXd::Identity(3, 3);
		
		E = E + (D.array()*D.array()).sum();
	}
	return E;
}

size_t nnz(const mxArray* m) { return *(mxGetJc(m) + mxGetN(m)); }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *points, *l_target, *value, *col, *row, *faces, *d, *l_temp,*points_to_faces,*edge_length,*point_to_points;
	//int ;
	bool transpose;
	int RowA, ColA, RowB, ColB;
	if (nlhs != 1) {
		mexErrMsgTxt("One output required.");
	}
	else if ( nrhs != 7) {
		mexErrMsgTxt("seven input required.");
	}
    


	//获得矩阵的列数
	int max_col= mxGetN(prhs[3]);
	int face_number= mxGetM(prhs[0]);
	//获取指向输入参数的指针
	faces = mxGetPr(prhs[0]);
	size_t point_number = int(mxGetScalar(prhs[1]));
	double *a = mxGetPr(prhs[2]);
	points_to_faces = mxGetPr(prhs[3]);
	edge_length = mxGetPr(prhs[4]);
	point_to_points = mxGetPr(prhs[5]);
	auto A_index=Eigen::Map<const SpMat>(mxGetM(prhs[6]), mxGetN(prhs[6]), nnz(prhs[6]), (const index_t*)mxGetJc(prhs[6]), (const index_t*)mxGetIr(prhs[6]), mxGetPr(prhs[6]));



		

	Eigen::SparseMatrix < double >  Mk(3 * point_number, 3 * point_number);
	
	
	//生成输出参量的mxArray
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *length;
	length = mxGetPr(plhs[0]);
	*length= energy_a(faces, face_number, point_number, a, points_to_faces, edge_length, point_to_points, max_col, A_index);


}


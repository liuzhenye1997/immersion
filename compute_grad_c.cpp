
# include "mex.h"
#include "D:\test2\拟牛顿法\2\Eigen\Eigen\Dense"
//#include <Eigen/Dense>
#include "D:\test2\拟牛顿法\2\Eigen\Eigen\Eigen"
#include<cmath>
using namespace std;
//非负数返回1，否则返回-1
double sign(double x)
{
	if (x >= 0)
	{
		return 1;
	}
	else
	{
		return -1;
	}
}

//计算各种能量的梯度
void compute_grad(const int face_number, const int point_number, const double * points_temp, const double * faces, const double * l_target, double * l_temp, double *grad, int energy_type = 1, double p = 2,double w=0, double *initial_points=nullptr)
{
	//四次能量
	if (energy_type == 1)
	{
#pragma omp parallel for 
		for (int k = 0; k < 3; k++)
		{
			for (int i = 0; i < face_number; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int ij = faces[i + j * face_number] - 1;
					int ij1 = faces[i + (j + 1) % 3 * face_number] - 1;
					grad[ij + k * point_number] = grad[ij + k * point_number] + 4 * (pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2))*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]);
					grad[ij1 + k * point_number] = grad[ij1 + k * point_number] - 4 * (pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2))*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]);
				}
			}
		}
	}
	else if (energy_type == 2)
	{
#pragma omp parallel for 
		for (int k = 0; k < 3; k++)
		{
			for (int i = 0; i < face_number; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int ij = faces[i + j * face_number] - 1;
					int ij1 = faces[i + (j + 1) % 3 * face_number] - 1;
					grad[ij + k * point_number] = grad[ij + k * point_number] + sign(pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2)) / (sqrt(abs(pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2))))*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]);
					grad[ij1 + k * point_number] = grad[ij1 + k * point_number] - sign(pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2)) / (sqrt(abs(pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2))))*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]);
				}
			}
		}
	}
	else if (energy_type == 3)
	{
#pragma omp parallel for 
		for (int k = 0; k < 3; k++)
		{
			for (int i = 0; i < face_number; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int ij = faces[i + j * face_number] - 1;
					int ij1 = faces[i + (j + 1) % 3 * face_number] - 1;
					grad[ij + k * point_number] = grad[ij + k * point_number] + 2 * (l_temp[i + j * face_number] - l_target[i + j * face_number])*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]) / l_temp[i + j * face_number];
					grad[ij1 + k * point_number] = grad[ij1 + k * point_number] - 2 * (l_temp[i + j * face_number] - l_target[i + j * face_number])*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]) / l_temp[i + j * face_number];
				}
			}
		}
	}
	else if (energy_type == 4)
	{
#pragma omp parallel for 
		for (int k = 0; k < 3; k++)
		{
			for (int i = 0; i < face_number; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int ij = faces[i + j * face_number] - 1;
					int ij1 = faces[i + (j + 1) % 3 * face_number] - 1;
					grad[ij + k * point_number] = grad[ij + k * point_number] + 2 * (l_temp[i + j * face_number] - l_target[i + j * face_number])*exp(pow(l_temp[i + j * face_number] - l_target[i + j * face_number], 2))*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]) / l_temp[i + j * face_number];
					grad[ij1 + k * point_number] = grad[ij1 + k * point_number] - 2 * (l_temp[i + j * face_number] - l_target[i + j * face_number])*exp(pow(l_temp[i + j * face_number] - l_target[i + j * face_number], 2))*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]) / l_temp[i + j * face_number];
				}
			}
		}
	}
	//L2能量
	else if (energy_type == 5)
	{
#pragma omp parallel for 
		for (int k = 0; k < 3; k++)
		{
			for (int i = 0; i < face_number; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int ij = faces[i + j * face_number] - 1;
					int ij1 = faces[i + (j + 1) % 3 * face_number] - 1;
					grad[ij + k * point_number] = grad[ij + k * point_number] + 2 * (l_temp[i + j * face_number] - l_target[i + j * face_number])*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]) / l_temp[i + j * face_number] / l_target[i + j * face_number] / l_target[i + j * face_number];
					grad[ij1 + k * point_number] = grad[ij1 + k * point_number] - 2 * (l_temp[i + j * face_number] - l_target[i + j * face_number])*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]) / l_temp[i + j * face_number] / l_target[i + j * face_number] / l_target[i + j * face_number];
				}
			}
		}
	}
	//LP能量
	else if (energy_type == 6)
	{
#pragma omp parallel for 
		for (int k = 0; k < 3; k++)
		{
			for (int i = 0; i < face_number; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int ij = faces[i + j * face_number] - 1;
					int ij1 = faces[i + (j + 1) % 3 * face_number] - 1;
					grad[ij + k * point_number] = grad[ij + k * point_number] + p * pow(l_temp[i + j * face_number] - l_target[i + j * face_number], p - 1)*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]) / l_temp[i + j * face_number] / pow(l_target[i + j * face_number], p);
					grad[ij1 + k * point_number] = grad[ij1 + k * point_number] - p * pow(l_temp[i + j * face_number] - l_target[i + j * face_number], p - 1)*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]) / l_temp[i + j * face_number] / pow(l_target[i + j * face_number], p);
				}
			}
		}
	}
	else if (energy_type == 7)
	{
	#pragma omp parallel for 
		for (int k = 0; k < 3; k++)
		{
			for (int i = 0; i < face_number; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int ij = faces[i + j * face_number] - 1;
					int ij1 = faces[i + (j + 1) % 3 * face_number] - 1;
					grad[ij + k * point_number] = grad[ij + k * point_number] + 4 * (pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2) )*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number])+w/12.0*(points_temp[ij + k * point_number]-initial_points[ij + k * point_number]);
					grad[ij1 + k * point_number] = grad[ij1 + k * point_number] - 4 * (pow(l_temp[i + j * face_number], 2) - pow(l_target[i + j * face_number], 2) )*(points_temp[ij + k * point_number] - points_temp[ij1 + k * point_number]) + w / 12.0*(points_temp[ij1 + k * point_number] - initial_points[ij1 + k * point_number]);
				}
			}
		}
	}

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *points, *l_target, *value, *col, *row, *faces, *d, *l_temp;
	//int ;
	bool transpose;
	int RowA, ColA, RowB, ColB;
	if (nlhs != 1) {
		mexErrMsgTxt("One output required.");
	}
	else if (nrhs != 4 && nrhs != 5 && nrhs != 6 && nrhs != 8) {
		mexErrMsgTxt("Four or five or six input required.");
	}
	int energy_type;
	if (nrhs != 4)
	{
		energy_type = int(mxGetScalar(prhs[4]));
		if (energy_type != 0 && energy_type != 1 && energy_type != 2 && energy_type != 3 && energy_type != 4 && energy_type != 5 && energy_type != 6 && energy_type != 7 && energy_type != 8)
		{
			mexErrMsgTxt("energ_type必须等于0,1,2,3,4,5,6,7或8");
		}
	}

	//获得矩阵的行数
	size_t point_number = mxGetM(prhs[0]);


	//获得矩阵的列数
	ColA = mxGetN(prhs[0]);
	int face_number = mxGetM(prhs[1]);
	ColB = mxGetN(prhs[1]);

	//获取指向输入参数的指针
	points = mxGetPr(prhs[0]);
	faces = mxGetPr(prhs[1]);
	l_target = mxGetPr(prhs[2]);
	l_temp = mxGetPr(prhs[3]);
	//生成输出参量的mxArray
	plhs[0] = mxCreateDoubleMatrix(point_number, 3, mxREAL);
	double *grad;
	grad = mxGetPr(plhs[0]);

	if (nrhs == 5)
	{
		int energy_type = int(mxGetScalar(prhs[4]));
		compute_grad(face_number, point_number, points, faces, l_target, l_temp, grad, energy_type);
	}
	else if (nrhs == 4)
	{
		compute_grad(face_number, point_number, points, faces, l_target, l_temp, grad);
	}
	else if (nrhs==6)
	{
		int energy_type = int(mxGetScalar(prhs[4]));
		double p = int(mxGetScalar(prhs[5]));
		compute_grad(face_number, point_number, points, faces, l_target, l_temp, grad, energy_type, p);
	}
	else if (nrhs == 8)
	{
		int energy_type = int(mxGetScalar(prhs[4]));
		double p = int(mxGetScalar(prhs[5]));
		double w = double(mxGetScalar(prhs[6]));
		/*mexPrintf("%lf\t%lf\n",p,w);*/
		double *initial_points = mxGetPr(prhs[7]);
		compute_grad(face_number, point_number, points, faces, l_target, l_temp, grad, energy_type, p,w, initial_points);
	}

}


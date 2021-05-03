
#include "mex.h"
#include "D:\test2\拟牛顿法\2\Eigen\Eigen\Dense"
#include "D:\test2\拟牛顿法\2\Eigen\Eigen\Eigen"
#include<cmath>
using namespace std;
using index_t = int64_t;
using SpMat = Eigen::SparseMatrix<double, 0, index_t>;

void compute_Mk(double *faces,int face_number, int point_number, SpMat A, double *points_to_faces, double *length, double *point_to_points, int max_col, Eigen::SparseMatrix < double >&MS)
{
	auto Mk = vector<vector<Eigen::Matrix3d>>(point_number, vector<Eigen::Matrix3d>(max_col));

	for (int i = 0; i < point_number; i++)
	{
		Eigen::Matrix3d M;
		M.setIdentity();
		auto X = vector<Eigen::Matrix3d>(points_to_faces[i]);
		auto Z = vector<Eigen::Matrix3d>(points_to_faces[i]);

		for (int j = 0; j < points_to_faces[i]; j++)
		{
			double theta = 0;
			for (int k = 0; k < 3; k++)
			{
				if (faces[(int)points_to_faces[(j + 1)*point_number + i]-1 + k * face_number] == (i + 1))
				{
					int m = points_to_faces[i + (j + 1)*point_number] - 1;
					theta = acos((length[m + k * face_number] * length[m + k * face_number] + length[m + ((k + 2) % 3) * face_number] * length[m + ((k + 2) % 3) * face_number] - length[m + ((k + 1) % 3) * face_number] * length[m + ((k + 1) % 3) * face_number]) / (2 * length[m + k * face_number] * length[m + ((k + 2) % 3) * face_number]));
					break;
				}
			}

			Z[j] << cos(theta), -sin(theta), 0,
				sin(theta), cos(theta), 0,
				0, 0, 1;
			M = M * Z[j];
			theta=A.coeff(i, point_to_points[i + (((j + 1) % (int)points_to_faces[i]) + 1)*point_number] - 1);

			X[j] << 1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta);
	
			if ((j + 1) != points_to_faces[i])
			{
				M = M * X[j];
			}
		}
		Mk[i][0] = M;
		for (int j = 0; j < points_to_faces[i] - 1; j++)
		{
			int index = (j + (int)points_to_faces[i] - 1) % (int)points_to_faces[i];
			M = M * X[index];		
			M = M * Z[j];
			M = Z[j].transpose()*M;

			M = X[j].transpose() * M;
			index = j + 1;
			Mk[i][index] = M;		
		}	
	}

	std::vector < Eigen::Triplet < double > > triplets;
	for (int i = 0; i < point_number; i++)
	{
		for (int j = 0; j < points_to_faces[i]; j++)
		{
			for (int m = 0; m < 3; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					triplets.emplace_back(3 * i+m, 3 * ((int)point_to_points[i + (j + 1)*point_number]-1)+n, Mk[i][j](m, n));
				}
			}		
		}
	}

	MS.setFromTriplets(triplets.begin(), triplets.end());
}

size_t nnz(const mxArray* m) { return *(mxGetJc(m) + mxGetN(m)); }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *points, *l_target, *value, *col, *row, *faces, *d, *l_temp,*points_to_faces,*edge_length,*point_to_points;
	bool transpose;
	int RowA, ColA, RowB, ColB;
	if (nlhs != 1) {
		mexErrMsgTxt("One output required.");
	}
	else if ( nrhs != 6) {
		mexErrMsgTxt("six input required.");
	}
    


	//获得矩阵的列数
	int max_col= mxGetN(prhs[3]);
	int face_number= mxGetM(prhs[0]);
	//获取指向输入参数的指针
	faces = mxGetPr(prhs[0]);
	size_t point_number = int(mxGetScalar(prhs[1]));
	auto A=Eigen::Map<const SpMat>(mxGetM(prhs[2]), mxGetN(prhs[2]), nnz(prhs[2]), (const index_t*)mxGetJc(prhs[2]), (const index_t*)mxGetIr(prhs[2]), mxGetPr(prhs[2]));


	points_to_faces = mxGetPr(prhs[3]);
	edge_length= mxGetPr(prhs[4]);
	point_to_points= mxGetPr(prhs[5]);

	Eigen::SparseMatrix < double >  Mk(3 * point_number, 3 * point_number);
	compute_Mk(faces,face_number,point_number,A,points_to_faces, edge_length,point_to_points,max_col,Mk);
	
	//生成输出参量的mxArray
	plhs[0] = mxCreateSparse(3 * point_number, 3 * point_number, Mk.nonZeros(), mxREAL);
	double *pr, *pi, *si, *sr;
	mwIndex *irs, *jcs;
	sr = mxGetPr(plhs[0]);
	irs = mxGetIr(plhs[0]);
	jcs = mxGetJc(plhs[0]);
	int *innerIndex = Mk.innerIndexPtr();
	double *valueptr = Mk.valuePtr();
	for (int i = 0; i < Mk.nonZeros(); i++)
	{
		sr[i] = valueptr[i];
	}

	for (int i = 0; i < Mk.nonZeros(); i++)
	{
		irs[i] = innerIndex[i];
	}

	int *OuterStarts = Mk.outerIndexPtr();
	for (int i = 0; i < 3 * point_number + 1; i++)
	{
		jcs[i] = OuterStarts[i];
	}

}


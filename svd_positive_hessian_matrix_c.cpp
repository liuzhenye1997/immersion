#include"math.h"
# include "mex.h"
#include "D:\test2\拟牛顿法\2\Eigen\Eigen\Dense"
//#include <Eigen/Dense>
#include "D:\test2\拟牛顿法\2\Eigen\Eigen\Eigen"

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

void hessian_matrix_triplet(const int face_number, const int point_number, const double * points, const double * faces, const double * l_target, double * col, double *row, double *value, int energy_type,double p=2)
{
	int i, j, k, l;
	if (energy_type == 1)
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
				double ld = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) - 4 * l*l;
				double ld2 = 12 * pow((xi - xj), 2.0) + 12 * pow((yi - yj), 2.0) + 12 * pow((zi - zj), 2.0) - 4 * l*l;
				double norm = sqrt(pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0));
				double x = (xi - xj) / norm;
				double y = (yi - yj) / norm;
				double z = (zi - zj) / norm;

				row[i * 108 + j * 36] = xi_index; col[i * 108 + j * 36] = xi_index; value[i * 108 + j * 36] = ld + 2 * (ld < 0)*abs(ld)*(y*y + z * z) + 2 * pow(2 * xi - 2 * xj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*x *x;
				row[i * 108 + j * 36 + 1] = yi_index; col[i * 108 + j * 36 + 1] = yi_index; value[i * 108 + j * 36 + 1] = ld + 2 * (ld < 0)*abs(ld)*(x*x + z * z) + 2 * pow(2 * yi - 2 * yj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*y*y;
				row[i * 108 + j * 36 + 2] = zi_index; col[i * 108 + j * 36 + 2] = zi_index; value[i * 108 + j * 36 + 2] = ld + 2 * (ld < 0)*abs(ld)*(y*y + x * x) + 2 * pow(2 * zi - 2 * zj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*z *z;
				row[i * 108 + j * 36 + 3] = yi_index; col[i * 108 + j * 36 + 3] = xi_index; value[i * 108 + j * 36 + 3] = -2 * x * y*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 4] = zi_index; col[i * 108 + j * 36 + 4] = xi_index; value[i * 108 + j * 36 + 4] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 5] = zi_index; col[i * 108 + j * 36 + 5] = yi_index; value[i * 108 + j * 36 + 5] = -2 * z * y*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				memcpy(value + i * 108 + j * 36 + 6, value + i * 108 + j * 36 + 3, 3 * sizeof(double));
				row[i * 108 + j * 36 + 6] = xi_index; col[i * 108 + j * 36 + 6] = yi_index; //value[i * 108 + j * 36 + 6] = -2 * x * y*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 7] = xi_index; col[i * 108 + j * 36 + 7] = zi_index; //value[i * 108 + j * 36 + 9] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 8] = yi_index; col[i * 108 + j * 36 + 8] = zi_index;// value[i * 108 + j * 36 + 10] = -2 * z * y*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				memcpy(value + i * 108 + j * 36 + 9, value + i * 108 + j * 36, 9 * sizeof(double));
				row[i * 108 + j * 36 + 9] = xj_index; col[i * 108 + j * 36 + 9] = xj_index;//value[i * 108 + j * 36 + 1] = ld + 2*(ld < 0)*abs(ld)*(y*y + z*z) + 2 * pow(2 * xi - 2 * xj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*x*x;			
				row[i * 108 + j * 36 + 10] = yj_index; col[i * 108 + j * 36 + 10] = yj_index;//value[i * 108 + j * 36 + 3] = ld + 2 * (ld < 0)*abs(ld)*(x*x + z*z) + 2 * pow(2 * yi - 2 * yj, 2.0)+2 * (ld2 < 0)*abs(ld2)*y*y;	
				row[i * 108 + j * 36 + 11] = zj_index; col[i * 108 + j * 36 + 11] = zj_index;// value[i * 108 + j * 36 + 5] = ld + 2 * (ld < 0)*abs(ld)*(y*y + x*x) + 2 * pow(2 * zi - 2 * zj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*z *z;
				row[i * 108 + j * 36 + 12] = yj_index; col[i * 108 + j * 36 + 12] = xj_index; //value[i * 108 + j * 36 + 31] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 13] = zj_index; col[i * 108 + j * 36 + 13] = xj_index; //value[i * 108 + j * 36 + 33] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 14] = zj_index; col[i * 108 + j * 36 + 14] = yj_index; //value[i * 108 + j * 36 + 34] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				row[i * 108 + j * 36 + 15] = xj_index; col[i * 108 + j * 36 + 15] = yj_index; //value[i * 108 + j * 36 + 30] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 16] = xj_index; col[i * 108 + j * 36 + 16] = zj_index; //value[i * 108 + j * 36 + 32] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;			
				row[i * 108 + j * 36 + 17] = yj_index; col[i * 108 + j * 36 + 17] = zj_index; //value[i * 108 + j * 36 + 35] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;

				row[i * 108 + j * 36 + 18] = xi_index; col[i * 108 + j * 36 + 18] = xj_index; value[i * 108 + j * 36 + 18] = -ld - 2 * pow(2 * xi - 2 * xj, 2.0);
				row[i * 108 + j * 36 + 19] = yi_index; col[i * 108 + j * 36 + 19] = yj_index; value[i * 108 + j * 36 + 19] = -ld - 2 * pow(2 * yi - 2 * yj, 2);
				row[i * 108 + j * 36 + 20] = zi_index; col[i * 108 + j * 36 + 20] = zj_index; value[i * 108 + j * 36 + 20] = -ld - 2 * pow(2 * zi - 2 * zj, 2);
				row[i * 108 + j * 36 + 21] = yi_index; col[i * 108 + j * 36 + 21] = xj_index; value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 22] = zi_index; col[i * 108 + j * 36 + 22] = xj_index; value[i * 108 + j * 36 + 22] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 23] = zi_index; col[i * 108 + j * 36 + 23] = yj_index; value[i * 108 + j * 36 + 23] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
				memcpy(value + i * 108 + j * 36 + 24, value + i * 108 + j * 36 + 21, 3 * sizeof(double));
				row[i * 108 + j * 36 + 24] = xi_index; col[i * 108 + j * 36 + 24] = yj_index; //value[i * 108 + j * 36 + 18] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 25] = xi_index; col[i * 108 + j * 36 + 25] = zj_index; //value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 26] = yi_index; col[i * 108 + j * 36 + 26] = zj_index; //value[i * 108 + j * 36 + 24] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
				memcpy(value + i * 108 + j * 36 + 27, value + i * 108 + j * 36 + 18, 9 * sizeof(double));
				row[i * 108 + j * 36 + 27] = xj_index; col[i * 108 + j * 36 + 27] = xi_index; //value[i * 108 + j * 36 + 13] = -ld - 2 * pow(2 * xi - 2 * xj, 2);			
				row[i * 108 + j * 36 + 28] = yj_index; col[i * 108 + j * 36 + 28] = yi_index; //value[i * 108 + j * 36 + 15] = -ld - 2 * pow(2 * yi - 2 * yj, 2);			
				row[i * 108 + j * 36 + 29] = zj_index; col[i * 108 + j * 36 + 29] = zi_index; //value[i * 108 + j * 36 + 17] = -ld - 2 * pow(2 * zi - 2 * zj, 2);
				row[i * 108 + j * 36 + 30] = yj_index; col[i * 108 + j * 36 + 30] = xi_index; //value[i * 108 + j * 36 + 19] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 31] = zj_index; col[i * 108 + j * 36 + 31] = xi_index; //value[i * 108 + j * 36 + 20] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 32] = zj_index; col[i * 108 + j * 36 + 32] = yi_index; //value[i * 108 + j * 36 + 25] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);											
				row[i * 108 + j * 36 + 33] = xj_index; col[i * 108 + j * 36 + 33] = yi_index; //value[i * 108 + j * 36 + 23] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);						
				row[i * 108 + j * 36 + 34] = xj_index; col[i * 108 + j * 36 + 34] = zi_index; //value[i * 108 + j * 36 + 26] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);			
				row[i * 108 + j * 36 + 35] = yj_index; col[i * 108 + j * 36 + 35] = zi_index; //value[i * 108 + j * 36 + 28] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
			}
		}
	}
	else if (energy_type == 2)
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
				double ld = pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0);
				double norm = sqrt(pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0));
				double x = (xi - xj) / norm;
				double y = (yi - yj) / norm;
				double z = (zi - zj) / norm;
				double lambda = (norm - l) / norm * 2;

				Eigen::Matrix3d A;
				A << sign(ld - l * l) / pow(abs(ld - l * l), 0.5) - (1 * 4 * (xi - xj)*(xi - xj)) / (4 * pow(abs(ld - l * l), 1.5)),
					-(1 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (4 * pow(abs(ld - l * l), 1.5)),
					-(1 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (4 * pow(abs(ld - l * l), 1.5)),
					-(1 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (4 * pow(abs(ld - l * l), 1.5)),
					sign(ld - l * l) / pow(abs(ld - l * l), 0.5) - (1 * 4 * (yi - yj)*(yi - yj)) / (4 * pow(abs(ld - l * l), 1.5)),
					-(1 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (4 * pow(abs(ld - l * l), 1.5)),
					-(1 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (4 * pow(abs(ld - l * l), 1.5)),
					-(1 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (4 * pow(abs(ld - l * l), 1.5)),
					sign(ld - l * l) / pow(abs(ld - l * l), 0.5) - (1 * 4 * (zi - zj)*(zi - zj)) / (4 * pow(abs(ld - l * l), 1.5));
				Eigen::JacobiSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
				Eigen::Matrix3d V = svd.matrixV();
				Eigen::Matrix3d U = svd.matrixU();
				Eigen::Vector3d d = svd.singularValues();
				for (int k = 0; k < 3; k++)
				{
					if (V(0, k)*U(0, k) + V(1, k)*U(1, k) + V(2, k)*U(2, k) < 0)
					{
						d(k) = -d(k);
					}
					if (d(k) < 0)
					{
						d(k) = 0;
					}
				}
				Eigen::Matrix3d S = d.asDiagonal();
				Eigen::Matrix3d H = U * S*U.transpose();

				row[i * 108 + j * 36] = xi_index; col[i * 108 + j * 36] = xi_index; value[i * 108 + j * 36] = H(0, 0);
				row[i * 108 + j * 36 + 1] = yi_index; col[i * 108 + j * 36 + 1] = yi_index; value[i * 108 + j * 36 + 1] = H(1, 1);
				row[i * 108 + j * 36 + 2] = zi_index; col[i * 108 + j * 36 + 2] = zi_index; value[i * 108 + j * 36 + 2] = H(2, 2);
				row[i * 108 + j * 36 + 3] = yi_index; col[i * 108 + j * 36 + 3] = xi_index; value[i * 108 + j * 36 + 3] = H(1, 0);
				row[i * 108 + j * 36 + 4] = zi_index; col[i * 108 + j * 36 + 4] = xi_index; value[i * 108 + j * 36 + 4] = H(2, 0);
				row[i * 108 + j * 36 + 5] = zi_index; col[i * 108 + j * 36 + 5] = yi_index; value[i * 108 + j * 36 + 5] = H(2, 1);
				memcpy(value + i * 108 + j * 36 + 6, value + i * 108 + j * 36 + 3, 3 * sizeof(double));
				row[i * 108 + j * 36 + 6] = xi_index; col[i * 108 + j * 36 + 6] = yi_index; //value[i * 108 + j * 36 + 6] = -2 * x * y*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 7] = xi_index; col[i * 108 + j * 36 + 7] = zi_index; //value[i * 108 + j * 36 + 9] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 8] = yi_index; col[i * 108 + j * 36 + 8] = zi_index;// value[i * 108 + j * 36 + 10] = -2 * z * y*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				memcpy(value + i * 108 + j * 36 + 9, value + i * 108 + j * 36, 9 * sizeof(double));
				row[i * 108 + j * 36 + 9] = xj_index; col[i * 108 + j * 36 + 9] = xj_index;//value[i * 108 + j * 36 + 1] = ld + 2*(ld < 0)*abs(ld)*(y*y + z*z) + 2 * pow(2 * xi - 2 * xj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*x*x;			
				row[i * 108 + j * 36 + 10] = yj_index; col[i * 108 + j * 36 + 10] = yj_index;//value[i * 108 + j * 36 + 3] = ld + 2 * (ld < 0)*abs(ld)*(x*x + z*z) + 2 * pow(2 * yi - 2 * yj, 2.0)+2 * (ld2 < 0)*abs(ld2)*y*y;	
				row[i * 108 + j * 36 + 11] = zj_index; col[i * 108 + j * 36 + 11] = zj_index;// value[i * 108 + j * 36 + 5] = ld + 2 * (ld < 0)*abs(ld)*(y*y + x*x) + 2 * pow(2 * zi - 2 * zj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*z *z;
				row[i * 108 + j * 36 + 12] = yj_index; col[i * 108 + j * 36 + 12] = xj_index; //value[i * 108 + j * 36 + 31] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 13] = zj_index; col[i * 108 + j * 36 + 13] = xj_index; //value[i * 108 + j * 36 + 33] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 14] = zj_index; col[i * 108 + j * 36 + 14] = yj_index; //value[i * 108 + j * 36 + 34] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				row[i * 108 + j * 36 + 15] = xj_index; col[i * 108 + j * 36 + 15] = yj_index; //value[i * 108 + j * 36 + 30] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 16] = xj_index; col[i * 108 + j * 36 + 16] = zj_index; //value[i * 108 + j * 36 + 32] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;			
				row[i * 108 + j * 36 + 17] = yj_index; col[i * 108 + j * 36 + 17] = zj_index; //value[i * 108 + j * 36 + 35] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;

				row[i * 108 + j * 36 + 18] = xi_index; col[i * 108 + j * 36 + 18] = xj_index; value[i * 108 + j * 36 + 18] = -H(0, 0);
				row[i * 108 + j * 36 + 19] = yi_index; col[i * 108 + j * 36 + 19] = yj_index; value[i * 108 + j * 36 + 19] = -H(1, 1);
				row[i * 108 + j * 36 + 20] = zi_index; col[i * 108 + j * 36 + 20] = zj_index; value[i * 108 + j * 36 + 20] = -H(2, 2);
				row[i * 108 + j * 36 + 21] = yi_index; col[i * 108 + j * 36 + 21] = xj_index; value[i * 108 + j * 36 + 21] = -H(0, 1);
				row[i * 108 + j * 36 + 22] = zi_index; col[i * 108 + j * 36 + 22] = xj_index; value[i * 108 + j * 36 + 22] = -H(0, 2);
				row[i * 108 + j * 36 + 23] = zi_index; col[i * 108 + j * 36 + 23] = yj_index; value[i * 108 + j * 36 + 23] = -H(1, 2);
				memcpy(value + i * 108 + j * 36 + 24, value + i * 108 + j * 36 + 21, 3 * sizeof(double));
				row[i * 108 + j * 36 + 24] = xi_index; col[i * 108 + j * 36 + 24] = yj_index; //value[i * 108 + j * 36 + 18] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 25] = xi_index; col[i * 108 + j * 36 + 25] = zj_index; //value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 26] = yi_index; col[i * 108 + j * 36 + 26] = zj_index; //value[i * 108 + j * 36 + 24] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
				memcpy(value + i * 108 + j * 36 + 27, value + i * 108 + j * 36 + 18, 9 * sizeof(double));
				row[i * 108 + j * 36 + 27] = xj_index; col[i * 108 + j * 36 + 27] = xi_index; //value[i * 108 + j * 36 + 13] = -ld - 2 * pow(2 * xi - 2 * xj, 2);			
				row[i * 108 + j * 36 + 28] = yj_index; col[i * 108 + j * 36 + 28] = yi_index; //value[i * 108 + j * 36 + 15] = -ld - 2 * pow(2 * yi - 2 * yj, 2);			
				row[i * 108 + j * 36 + 29] = zj_index; col[i * 108 + j * 36 + 29] = zi_index; //value[i * 108 + j * 36 + 17] = -ld - 2 * pow(2 * zi - 2 * zj, 2);
				row[i * 108 + j * 36 + 30] = yj_index; col[i * 108 + j * 36 + 30] = xi_index; //value[i * 108 + j * 36 + 19] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 31] = zj_index; col[i * 108 + j * 36 + 31] = xi_index; //value[i * 108 + j * 36 + 20] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 32] = zj_index; col[i * 108 + j * 36 + 32] = yi_index; //value[i * 108 + j * 36 + 25] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);											
				row[i * 108 + j * 36 + 33] = xj_index; col[i * 108 + j * 36 + 33] = yi_index; //value[i * 108 + j * 36 + 23] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);						
				row[i * 108 + j * 36 + 34] = xj_index; col[i * 108 + j * 36 + 34] = zi_index; //value[i * 108 + j * 36 + 26] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);			
				row[i * 108 + j * 36 + 35] = yj_index; col[i * 108 + j * 36 + 35] = zi_index; //value[i * 108 + j * 36 + 28] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);			
			}
		}
	}
	else if (energy_type == 3)
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
				double ld = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) - 4 * l*l;
				double ld2 = 12 * pow((xi - xj), 2.0) + 12 * pow((yi - yj), 2.0) + 12 * pow((zi - zj), 2.0) - 4 * l*l;
				double norm = sqrt(pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0));
				double x = (xi - xj) / norm;
				double y = (yi - yj) / norm;
				double z = (zi - zj) / norm;
				double lambda = (norm - l) / norm * 2;

				row[i * 108 + j * 36] = xi_index; col[i * 108 + j * 36] = xi_index; value[i * 108 + j * 36] = 4 * (xi - xj)*(xi - xj) / (2 * (norm*norm)) - (2 * (l - norm)) / norm + (4 * (xi - xj)*(xi - xj) * (l - norm)) / (2 * norm*norm*norm) + 2 * (lambda < 0)*abs(lambda)*(1 - x * x);
				row[i * 108 + j * 36 + 1] = yi_index; col[i * 108 + j * 36 + 1] = yi_index; value[i * 108 + j * 36 + 1] = 4 * (yi - yj)*(yi - yj) / (2 * (norm*norm)) - (2 * (l - norm)) / norm + (4 * (yi - yj)*(yi - yj) * (l - norm)) / (2 * norm*norm*norm) + 2 * (lambda < 0)*abs(lambda)*(1 - y * y);
				row[i * 108 + j * 36 + 2] = zi_index; col[i * 108 + j * 36 + 2] = zi_index; value[i * 108 + j * 36 + 2] = 4 * (zi - zj)*(zi - zj) / (2 * (norm*norm)) - (2 * (l - norm)) / norm + (4 * (zi - zj)*(zi - zj) * (l - norm)) / (2 * norm*norm*norm) + 2 * (lambda < 0)*abs(lambda)*(1 - z * z);
				row[i * 108 + j * 36 + 3] = yi_index; col[i * 108 + j * 36 + 3] = xi_index; value[i * 108 + j * 36 + 3] = ((2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (2 * (norm*norm)) + ((2 * xi - 2 * xj)*(2 * yi - 2 * yj)*(l - norm)) / (2 * norm*norm*norm) - 2 * (lambda < 0)*abs(lambda)*x*y;
				row[i * 108 + j * 36 + 4] = zi_index; col[i * 108 + j * 36 + 4] = xi_index; value[i * 108 + j * 36 + 4] = ((2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) + ((2 * xi - 2 * xj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) - 2 * (lambda < 0)*abs(lambda)*x*z;
				row[i * 108 + j * 36 + 5] = zi_index; col[i * 108 + j * 36 + 5] = yi_index; value[i * 108 + j * 36 + 5] = ((2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) + ((2 * yi - 2 * yj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) - 2 * (lambda < 0)*abs(lambda)*z*y;
				memcpy(value + i * 108 + j * 36 + 6, value + i * 108 + j * 36 + 3, 3 * sizeof(double));
				row[i * 108 + j * 36 + 6] = xi_index; col[i * 108 + j * 36 + 6] = yi_index; //value[i * 108 + j * 36 + 6] = -2 * x * y*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 7] = xi_index; col[i * 108 + j * 36 + 7] = zi_index; //value[i * 108 + j * 36 + 9] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 8] = yi_index; col[i * 108 + j * 36 + 8] = zi_index;// value[i * 108 + j * 36 + 10] = -2 * z * y*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				memcpy(value + i * 108 + j * 36 + 9, value + i * 108 + j * 36, 9 * sizeof(double));
				row[i * 108 + j * 36 + 9] = xj_index; col[i * 108 + j * 36 + 9] = xj_index;//value[i * 108 + j * 36 + 1] = ld + 2*(ld < 0)*abs(ld)*(y*y + z*z) + 2 * pow(2 * xi - 2 * xj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*x*x;			
				row[i * 108 + j * 36 + 10] = yj_index; col[i * 108 + j * 36 + 10] = yj_index;//value[i * 108 + j * 36 + 3] = ld + 2 * (ld < 0)*abs(ld)*(x*x + z*z) + 2 * pow(2 * yi - 2 * yj, 2.0)+2 * (ld2 < 0)*abs(ld2)*y*y;	
				row[i * 108 + j * 36 + 11] = zj_index; col[i * 108 + j * 36 + 11] = zj_index;// value[i * 108 + j * 36 + 5] = ld + 2 * (ld < 0)*abs(ld)*(y*y + x*x) + 2 * pow(2 * zi - 2 * zj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*z *z;
				row[i * 108 + j * 36 + 12] = yj_index; col[i * 108 + j * 36 + 12] = xj_index; //value[i * 108 + j * 36 + 31] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 13] = zj_index; col[i * 108 + j * 36 + 13] = xj_index; //value[i * 108 + j * 36 + 33] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 14] = zj_index; col[i * 108 + j * 36 + 14] = yj_index; //value[i * 108 + j * 36 + 34] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				row[i * 108 + j * 36 + 15] = xj_index; col[i * 108 + j * 36 + 15] = yj_index; //value[i * 108 + j * 36 + 30] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 16] = xj_index; col[i * 108 + j * 36 + 16] = zj_index; //value[i * 108 + j * 36 + 32] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;			
				row[i * 108 + j * 36 + 17] = yj_index; col[i * 108 + j * 36 + 17] = zj_index; //value[i * 108 + j * 36 + 35] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;

				row[i * 108 + j * 36 + 18] = xi_index; col[i * 108 + j * 36 + 18] = xj_index; value[i * 108 + j * 36 + 18] = (2 * (l - norm)) / norm - 4 * (xi - xj)*(xi - xj) / (2 * (norm*norm)) - (4 * (xi - xj)*(xi - xj) * (l - norm)) / (2 * norm*norm*norm);
				row[i * 108 + j * 36 + 19] = yi_index; col[i * 108 + j * 36 + 19] = yj_index; value[i * 108 + j * 36 + 19] = (2 * (l - norm)) / norm - 4 * (yi - yj)*(yi - yj) / (2 * (norm*norm)) - (4 * (yi - yj)*(yi - yj) * (l - norm)) / (2 * norm*norm*norm);
				row[i * 108 + j * 36 + 20] = zi_index; col[i * 108 + j * 36 + 20] = zj_index; value[i * 108 + j * 36 + 20] = (2 * (l - norm)) / norm - 4 * (zi - zj)*(zi - zj) / (2 * (norm*norm)) - (4 * (zi - zj)*(zi - zj) * (l - norm)) / (2 * norm*norm*norm);
				row[i * 108 + j * 36 + 21] = yi_index; col[i * 108 + j * 36 + 21] = xj_index; value[i * 108 + j * 36 + 21] = -((2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (2 * (norm*norm)) - ((2 * xi - 2 * xj)*(2 * yi - 2 * yj)*(l - norm)) / (2 * norm*norm*norm);
				row[i * 108 + j * 36 + 22] = zi_index; col[i * 108 + j * 36 + 22] = xj_index; value[i * 108 + j * 36 + 22] = -((2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) - ((2 * xi - 2 * xj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm);
				row[i * 108 + j * 36 + 23] = zi_index; col[i * 108 + j * 36 + 23] = yj_index; value[i * 108 + j * 36 + 23] = -((2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) - ((2 * yi - 2 * yj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm);
				memcpy(value + i * 108 + j * 36 + 24, value + i * 108 + j * 36 + 21, 3 * sizeof(double));
				row[i * 108 + j * 36 + 24] = xi_index; col[i * 108 + j * 36 + 24] = yj_index; //value[i * 108 + j * 36 + 18] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 25] = xi_index; col[i * 108 + j * 36 + 25] = zj_index; //value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 26] = yi_index; col[i * 108 + j * 36 + 26] = zj_index; //value[i * 108 + j * 36 + 24] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
				memcpy(value + i * 108 + j * 36 + 27, value + i * 108 + j * 36 + 18, 9 * sizeof(double));
				row[i * 108 + j * 36 + 27] = xj_index; col[i * 108 + j * 36 + 27] = xi_index; //value[i * 108 + j * 36 + 13] = -ld - 2 * pow(2 * xi - 2 * xj, 2);			
				row[i * 108 + j * 36 + 28] = yj_index; col[i * 108 + j * 36 + 28] = yi_index; //value[i * 108 + j * 36 + 15] = -ld - 2 * pow(2 * yi - 2 * yj, 2);			
				row[i * 108 + j * 36 + 29] = zj_index; col[i * 108 + j * 36 + 29] = zi_index; //value[i * 108 + j * 36 + 17] = -ld - 2 * pow(2 * zi - 2 * zj, 2);
				row[i * 108 + j * 36 + 30] = yj_index; col[i * 108 + j * 36 + 30] = xi_index; //value[i * 108 + j * 36 + 19] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 31] = zj_index; col[i * 108 + j * 36 + 31] = xi_index; //value[i * 108 + j * 36 + 20] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 32] = zj_index; col[i * 108 + j * 36 + 32] = yi_index; //value[i * 108 + j * 36 + 25] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);											
				row[i * 108 + j * 36 + 33] = xj_index; col[i * 108 + j * 36 + 33] = yi_index; //value[i * 108 + j * 36 + 23] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);						
				row[i * 108 + j * 36 + 34] = xj_index; col[i * 108 + j * 36 + 34] = zi_index; //value[i * 108 + j * 36 + 26] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);			
				row[i * 108 + j * 36 + 35] = yj_index; col[i * 108 + j * 36 + 35] = zi_index; //value[i * 108 + j * 36 + 28] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);

			}
		}
	}
	else if (energy_type == 4)
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
				double ld = pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0);
				double norm = sqrt(pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0));
				double x = (xi - xj) / norm;
				double y = (yi - yj) / norm;
				double z = (zi - zj) / norm;
				double lambda = (norm - l) / norm * 2;
				double d = l - sqrt(pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0));

				Eigen::Matrix3d A;
				A << (exp(d*d) * 4 * (xi - xj)*(xi - xj)) / (2 * (ld)) - (2 * exp(d*d)*d) / norm + (exp(d*d) * 4 * (xi - xj)*(xi - xj) * d*d) / (ld)+(exp(d*d) * 4 * (xi - xj)*(xi - xj) * d) / (2 * pow(norm, 3)),
					(exp(d*d)*(2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (2 * (ld)) + (exp(d*d)*(2 * xi - 2 * xj)*(2 * yi - 2 * yj)*d) / (2 * pow(norm, 3)) + (exp(d*d)*(2 * xi - 2 * xj)*(2 * yi - 2 * yj)*d*d) / (ld),
					(exp(d*d)*(2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (2 * (ld)) + (exp(d*d)*(2 * xi - 2 * xj)*(2 * zi - 2 * zj)*d) / (2 * pow(norm, 3)) + (exp(d*d)*(2 * xi - 2 * xj)*(2 * zi - 2 * zj)*d*d) / (ld),
					(exp(d*d)*(2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (2 * (ld)) + (exp(d*d)*(2 * xi - 2 * xj)*(2 * yi - 2 * yj)*d) / (2 * pow(norm, 3)) + (exp(d*d)*(2 * xi - 2 * xj)*(2 * yi - 2 * yj)*d*d) / (ld),
					(exp(d*d) * 4 * (yi - yj)*(yi - yj)) / (2 * (ld)) - (2 * exp(d*d)*d) / norm + (exp(d*d) * 4 * (yi - yj)*(yi - yj) * d*d) / (ld)+(exp(d*d) * 4 * (yi - yj)*(yi - yj) * d) / (2 * pow(norm, 3)),
					(exp(d*d)*(2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (2 * (ld)) + (exp(d*d)*(2 * yi - 2 * yj)*(2 * zi - 2 * zj)*d) / (2 * pow(norm, 3)) + (exp(d*d)*(2 * yi - 2 * yj)*(2 * zi - 2 * zj)*d*d) / (ld),
					(exp(d*d)*(2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (2 * (ld)) + (exp(d*d)*(2 * xi - 2 * xj)*(2 * zi - 2 * zj)*d) / (2 * pow(norm, 3)) + (exp(d*d)*(2 * xi - 2 * xj)*(2 * zi - 2 * zj)*d*d) / (ld),
					(exp(d*d)*(2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (2 * (ld)) + (exp(d*d)*(2 * yi - 2 * yj)*(2 * zi - 2 * zj)*d) / (2 * pow(norm, 3)) + (exp(d*d)*(2 * yi - 2 * yj)*(2 * zi - 2 * zj)*d*d) / (ld),
					(exp(d*d) * 4 * (zi - zj)*(zi - zj)) / (2 * (ld)) - (2 * exp(d*d)*d) / norm + (exp(d*d) * 4 * (zi - zj)*(zi - zj) * d*d) / (ld)+(exp(d*d) * 4 * (zi - zj)*(zi - zj) * d) / (2 * pow(norm, 3));
				Eigen::JacobiSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
				Eigen::Matrix3d V = svd.matrixV();
				Eigen::Matrix3d U = svd.matrixU();
				Eigen::Vector3d singularValues = svd.singularValues();
				for (int k = 0; k < 3; k++)
				{
					if (V(0, k)*U(0, k) + V(1, k)*U(1, k) + V(2, k)*U(2, k) < 0)
					{
						singularValues(k) = -singularValues(k);
					}
					if (singularValues(k) < 0)
					{
						singularValues(k) = 0;
					}
				}
				Eigen::Matrix3d S = singularValues.asDiagonal();
				Eigen::Matrix3d H = U * S*U.transpose();

				row[i * 108 + j * 36] = xi_index; col[i * 108 + j * 36] = xi_index; value[i * 108 + j * 36] = 2 * H(0, 0) - A(0, 0);
				row[i * 108 + j * 36 + 1] = yi_index; col[i * 108 + j * 36 + 1] = yi_index; value[i * 108 + j * 36 + 1] = 2 * H(1, 1) - A(1, 1);
				row[i * 108 + j * 36 + 2] = zi_index; col[i * 108 + j * 36 + 2] = zi_index; value[i * 108 + j * 36 + 2] = 2 * H(2, 2) - A(2, 2);
				row[i * 108 + j * 36 + 3] = yi_index; col[i * 108 + j * 36 + 3] = xi_index; value[i * 108 + j * 36 + 3] = 2 * H(1, 0) - A(1, 0);
				row[i * 108 + j * 36 + 4] = zi_index; col[i * 108 + j * 36 + 4] = xi_index; value[i * 108 + j * 36 + 4] = 2 * H(2, 0) - A(0, 2);
				row[i * 108 + j * 36 + 5] = zi_index; col[i * 108 + j * 36 + 5] = yi_index; value[i * 108 + j * 36 + 5] = 2 * H(2, 1) - A(1, 2);
				memcpy(value + i * 108 + j * 36 + 6, value + i * 108 + j * 36 + 3, 3 * sizeof(double));
				row[i * 108 + j * 36 + 6] = xi_index; col[i * 108 + j * 36 + 6] = yi_index; //value[i * 108 + j * 36 + 6] = -2 * x * y*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 7] = xi_index; col[i * 108 + j * 36 + 7] = zi_index; //value[i * 108 + j * 36 + 9] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 8] = yi_index; col[i * 108 + j * 36 + 8] = zi_index;// value[i * 108 + j * 36 + 10] = -2 * z * y*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				memcpy(value + i * 108 + j * 36 + 9, value + i * 108 + j * 36, 9 * sizeof(double));
				row[i * 108 + j * 36 + 9] = xj_index; col[i * 108 + j * 36 + 9] = xj_index;//value[i * 108 + j * 36 + 1] = ld + 2*(ld < 0)*abs(ld)*(y*y + z*z) + 2 * pow(2 * xi - 2 * xj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*x*x;			
				row[i * 108 + j * 36 + 10] = yj_index; col[i * 108 + j * 36 + 10] = yj_index;//value[i * 108 + j * 36 + 3] = ld + 2 * (ld < 0)*abs(ld)*(x*x + z*z) + 2 * pow(2 * yi - 2 * yj, 2.0)+2 * (ld2 < 0)*abs(ld2)*y*y;	
				row[i * 108 + j * 36 + 11] = zj_index; col[i * 108 + j * 36 + 11] = zj_index;// value[i * 108 + j * 36 + 5] = ld + 2 * (ld < 0)*abs(ld)*(y*y + x*x) + 2 * pow(2 * zi - 2 * zj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*z *z;
				row[i * 108 + j * 36 + 12] = yj_index; col[i * 108 + j * 36 + 12] = xj_index; //value[i * 108 + j * 36 + 31] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 13] = zj_index; col[i * 108 + j * 36 + 13] = xj_index; //value[i * 108 + j * 36 + 33] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 14] = zj_index; col[i * 108 + j * 36 + 14] = yj_index; //value[i * 108 + j * 36 + 34] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				row[i * 108 + j * 36 + 15] = xj_index; col[i * 108 + j * 36 + 15] = yj_index; //value[i * 108 + j * 36 + 30] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 16] = xj_index; col[i * 108 + j * 36 + 16] = zj_index; //value[i * 108 + j * 36 + 32] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;			
				row[i * 108 + j * 36 + 17] = yj_index; col[i * 108 + j * 36 + 17] = zj_index; //value[i * 108 + j * 36 + 35] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;

				row[i * 108 + j * 36 + 18] = xi_index; col[i * 108 + j * 36 + 18] = xj_index; value[i * 108 + j * 36 + 18] = -A(0, 0);
				row[i * 108 + j * 36 + 19] = yi_index; col[i * 108 + j * 36 + 19] = yj_index; value[i * 108 + j * 36 + 19] = -A(1, 1);
				row[i * 108 + j * 36 + 20] = zi_index; col[i * 108 + j * 36 + 20] = zj_index; value[i * 108 + j * 36 + 20] = -A(2, 2);
				row[i * 108 + j * 36 + 21] = yi_index; col[i * 108 + j * 36 + 21] = xj_index; value[i * 108 + j * 36 + 21] = -A(0, 1);
				row[i * 108 + j * 36 + 22] = zi_index; col[i * 108 + j * 36 + 22] = xj_index; value[i * 108 + j * 36 + 22] = -A(0, 2);
				row[i * 108 + j * 36 + 23] = zi_index; col[i * 108 + j * 36 + 23] = yj_index; value[i * 108 + j * 36 + 23] = -A(1, 2);
				memcpy(value + i * 108 + j * 36 + 24, value + i * 108 + j * 36 + 21, 3 * sizeof(double));
				row[i * 108 + j * 36 + 24] = xi_index; col[i * 108 + j * 36 + 24] = yj_index; //value[i * 108 + j * 36 + 18] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 25] = xi_index; col[i * 108 + j * 36 + 25] = zj_index; //value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 26] = yi_index; col[i * 108 + j * 36 + 26] = zj_index; //value[i * 108 + j * 36 + 24] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
				memcpy(value + i * 108 + j * 36 + 27, value + i * 108 + j * 36 + 18, 9 * sizeof(double));
				row[i * 108 + j * 36 + 27] = xj_index; col[i * 108 + j * 36 + 27] = xi_index; //value[i * 108 + j * 36 + 13] = -ld - 2 * pow(2 * xi - 2 * xj, 2);			
				row[i * 108 + j * 36 + 28] = yj_index; col[i * 108 + j * 36 + 28] = yi_index; //value[i * 108 + j * 36 + 15] = -ld - 2 * pow(2 * yi - 2 * yj, 2);			
				row[i * 108 + j * 36 + 29] = zj_index; col[i * 108 + j * 36 + 29] = zi_index; //value[i * 108 + j * 36 + 17] = -ld - 2 * pow(2 * zi - 2 * zj, 2);
				row[i * 108 + j * 36 + 30] = yj_index; col[i * 108 + j * 36 + 30] = xi_index; //value[i * 108 + j * 36 + 19] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 31] = zj_index; col[i * 108 + j * 36 + 31] = xi_index; //value[i * 108 + j * 36 + 20] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 32] = zj_index; col[i * 108 + j * 36 + 32] = yi_index; //value[i * 108 + j * 36 + 25] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);											
				row[i * 108 + j * 36 + 33] = xj_index; col[i * 108 + j * 36 + 33] = yi_index; //value[i * 108 + j * 36 + 23] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);						
				row[i * 108 + j * 36 + 34] = xj_index; col[i * 108 + j * 36 + 34] = zi_index; //value[i * 108 + j * 36 + 26] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);			
				row[i * 108 + j * 36 + 35] = yj_index; col[i * 108 + j * 36 + 35] = zi_index; //value[i * 108 + j * 36 + 28] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);			
			}
		}
	}
	else if (energy_type == 5)
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
				double ld = 4 * pow((xi - xj), 2.0) + 4 * pow((yi - yj), 2.0) + 4 * pow((zi - zj), 2.0) - 4 * l*l;
				double ld2 = 12 * pow((xi - xj), 2.0) + 12 * pow((yi - yj), 2.0) + 12 * pow((zi - zj), 2.0) - 4 * l*l;
				double norm = sqrt(pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0));
				double x = (xi - xj) / norm;
				double y = (yi - yj) / norm;
				double z = (zi - zj) / norm;
				double lambda = (norm - l) / norm * 2;
				double weight = 1.0 / (l*l);

				row[i * 108 + j * 36] = xi_index; col[i * 108 + j * 36] = xi_index; value[i * 108 + j * 36] = 4 * (xi - xj)*(xi - xj) / (2 * (norm*norm)) - (2 * (l - norm)) / norm + (4 * (xi - xj)*(xi - xj) * (l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*(1 - x * x);
				row[i * 108 + j * 36 + 1] = yi_index; col[i * 108 + j * 36 + 1] = yi_index; value[i * 108 + j * 36 + 1] = 4 * (yi - yj)*(yi - yj) / (2 * (norm*norm)) - (2 * (l - norm)) / norm + (4 * (yi - yj)*(yi - yj) * (l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*(1 - y * y);
				row[i * 108 + j * 36 + 2] = zi_index; col[i * 108 + j * 36 + 2] = zi_index; value[i * 108 + j * 36 + 2] = 4 * (zi - zj)*(zi - zj) / (2 * (norm*norm)) - (2 * (l - norm)) / norm + (4 * (zi - zj)*(zi - zj) * (l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*(1 - z * z);
				row[i * 108 + j * 36 + 3] = yi_index; col[i * 108 + j * 36 + 3] = xi_index; value[i * 108 + j * 36 + 3] = ((2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (2 * (norm*norm)) + ((2 * xi - 2 * xj)*(2 * yi - 2 * yj)*(l - norm)) / (2 * norm*norm*norm) - (lambda < 0)*abs(lambda)*x*y;
				row[i * 108 + j * 36 + 4] = zi_index; col[i * 108 + j * 36 + 4] = xi_index; value[i * 108 + j * 36 + 4] = ((2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) + ((2 * xi - 2 * xj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) - (lambda < 0)*abs(lambda)*x*z;
				row[i * 108 + j * 36 + 5] = zi_index; col[i * 108 + j * 36 + 5] = yi_index; value[i * 108 + j * 36 + 5] = ((2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) + ((2 * yi - 2 * yj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) - (lambda < 0)*abs(lambda)*z*y;
				memcpy(value + i * 108 + j * 36 + 6, value + i * 108 + j * 36 + 3, 3 * sizeof(double));
				row[i * 108 + j * 36 + 6] = xi_index; col[i * 108 + j * 36 + 6] = yi_index; //value[i * 108 + j * 36 + 6] = -2 * x * y*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 7] = xi_index; col[i * 108 + j * 36 + 7] = zi_index; //value[i * 108 + j * 36 + 9] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 8] = yi_index; col[i * 108 + j * 36 + 8] = zi_index;// value[i * 108 + j * 36 + 10] = -2 * z * y*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				memcpy(value + i * 108 + j * 36 + 9, value + i * 108 + j * 36, 9 * sizeof(double));
				row[i * 108 + j * 36 + 9] = xj_index; col[i * 108 + j * 36 + 9] = xj_index;//value[i * 108 + j * 36 + 1] = ld + 2*(ld < 0)*abs(ld)*(y*y + z*z) + 2 * pow(2 * xi - 2 * xj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*x*x;			
				row[i * 108 + j * 36 + 10] = yj_index; col[i * 108 + j * 36 + 10] = yj_index;//value[i * 108 + j * 36 + 3] = ld + 2 * (ld < 0)*abs(ld)*(x*x + z*z) + 2 * pow(2 * yi - 2 * yj, 2.0)+2 * (ld2 < 0)*abs(ld2)*y*y;	
				row[i * 108 + j * 36 + 11] = zj_index; col[i * 108 + j * 36 + 11] = zj_index;// value[i * 108 + j * 36 + 5] = ld + 2 * (ld < 0)*abs(ld)*(y*y + x*x) + 2 * pow(2 * zi - 2 * zj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*z *z;
				row[i * 108 + j * 36 + 12] = yj_index; col[i * 108 + j * 36 + 12] = xj_index; //value[i * 108 + j * 36 + 31] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 13] = zj_index; col[i * 108 + j * 36 + 13] = xj_index; //value[i * 108 + j * 36 + 33] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 14] = zj_index; col[i * 108 + j * 36 + 14] = yj_index; //value[i * 108 + j * 36 + 34] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				row[i * 108 + j * 36 + 15] = xj_index; col[i * 108 + j * 36 + 15] = yj_index; //value[i * 108 + j * 36 + 30] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 16] = xj_index; col[i * 108 + j * 36 + 16] = zj_index; //value[i * 108 + j * 36 + 32] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;			
				row[i * 108 + j * 36 + 17] = yj_index; col[i * 108 + j * 36 + 17] = zj_index; //value[i * 108 + j * 36 + 35] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;

				row[i * 108 + j * 36 + 18] = xi_index; col[i * 108 + j * 36 + 18] = xj_index; value[i * 108 + j * 36 + 18] = (2 * (l - norm)) / norm - 4 * (xi - xj)*(xi - xj) / (2 * (norm*norm)) - (4 * (xi - xj)*(xi - xj) * (l - norm)) / (2 * norm*norm*norm) - (lambda < 0)*abs(lambda)*(1 - x * x);
				row[i * 108 + j * 36 + 19] = yi_index; col[i * 108 + j * 36 + 19] = yj_index; value[i * 108 + j * 36 + 19] = (2 * (l - norm)) / norm - 4 * (yi - yj)*(yi - yj) / (2 * (norm*norm)) - (4 * (yi - yj)*(yi - yj) * (l - norm)) / (2 * norm*norm*norm) - (lambda < 0)*abs(lambda)*(1 - y * y);
				row[i * 108 + j * 36 + 20] = zi_index; col[i * 108 + j * 36 + 20] = zj_index; value[i * 108 + j * 36 + 20] = (2 * (l - norm)) / norm - 4 * (zi - zj)*(zi - zj) / (2 * (norm*norm)) - (4 * (zi - zj)*(zi - zj) * (l - norm)) / (2 * norm*norm*norm) - (lambda < 0)*abs(lambda)*(1 - z * z);
				row[i * 108 + j * 36 + 21] = yi_index; col[i * 108 + j * 36 + 21] = xj_index; value[i * 108 + j * 36 + 21] = -((2 * xi - 2 * xj)*(2 * yi - 2 * yj)) / (2 * (norm*norm)) - ((2 * xi - 2 * xj)*(2 * yi - 2 * yj)*(l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*x*y;
				row[i * 108 + j * 36 + 22] = zi_index; col[i * 108 + j * 36 + 22] = xj_index; value[i * 108 + j * 36 + 22] = -((2 * xi - 2 * xj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) - ((2 * xi - 2 * xj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*x*z;
				row[i * 108 + j * 36 + 23] = zi_index; col[i * 108 + j * 36 + 23] = yj_index; value[i * 108 + j * 36 + 23] = -((2 * yi - 2 * yj)*(2 * zi - 2 * zj)) / (2 * (norm*norm)) - ((2 * yi - 2 * yj)*(2 * zi - 2 * zj)*(l - norm)) / (2 * norm*norm*norm) + (lambda < 0)*abs(lambda)*z*y;
				memcpy(value + i * 108 + j * 36 + 24, value + i * 108 + j * 36 + 21, 3 * sizeof(double));
				row[i * 108 + j * 36 + 24] = xi_index; col[i * 108 + j * 36 + 24] = yj_index; //value[i * 108 + j * 36 + 18] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 25] = xi_index; col[i * 108 + j * 36 + 25] = zj_index; //value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 26] = yi_index; col[i * 108 + j * 36 + 26] = zj_index; //value[i * 108 + j * 36 + 24] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
				memcpy(value + i * 108 + j * 36 + 27, value + i * 108 + j * 36 + 18, 9 * sizeof(double));
				row[i * 108 + j * 36 + 27] = xj_index; col[i * 108 + j * 36 + 27] = xi_index; //value[i * 108 + j * 36 + 13] = -ld - 2 * pow(2 * xi - 2 * xj, 2);			
				row[i * 108 + j * 36 + 28] = yj_index; col[i * 108 + j * 36 + 28] = yi_index; //value[i * 108 + j * 36 + 15] = -ld - 2 * pow(2 * yi - 2 * yj, 2);			
				row[i * 108 + j * 36 + 29] = zj_index; col[i * 108 + j * 36 + 29] = zi_index; //value[i * 108 + j * 36 + 17] = -ld - 2 * pow(2 * zi - 2 * zj, 2);
				row[i * 108 + j * 36 + 30] = yj_index; col[i * 108 + j * 36 + 30] = xi_index; //value[i * 108 + j * 36 + 19] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 31] = zj_index; col[i * 108 + j * 36 + 31] = xi_index; //value[i * 108 + j * 36 + 20] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 32] = zj_index; col[i * 108 + j * 36 + 32] = yi_index; //value[i * 108 + j * 36 + 25] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);											
				row[i * 108 + j * 36 + 33] = xj_index; col[i * 108 + j * 36 + 33] = yi_index; //value[i * 108 + j * 36 + 23] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);						
				row[i * 108 + j * 36 + 34] = xj_index; col[i * 108 + j * 36 + 34] = zi_index; //value[i * 108 + j * 36 + 26] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);			
				row[i * 108 + j * 36 + 35] = yj_index; col[i * 108 + j * 36 + 35] = zi_index; //value[i * 108 + j * 36 + 28] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);

				for (int k = 0; k < 36; k++)
				{
					value[i * 108 + j * 36 + k] *= weight;
				}
			}
		}
	}
	else if (energy_type == 6)
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


				double norm = sqrt(pow((xi - xj), 2.0) + pow((yi - yj), 2.0) + pow((zi - zj), 2.0));
				double x = (xi - xj) / norm;
				double y = (yi - yj) / norm;
				double z = (zi - zj) / norm;
				double lambda = pow(norm - l, p - 1) / norm * p;
				double weight = 1.0 / (pow(l, p));
				double ld = pow(norm - l, p - 2);

				row[i * 108 + j * 36] = xi_index; col[i * 108 + j * 36] = xi_index; value[i * 108 + j * 36] = p * ld / norm * (norm - l + (xi - xj)*(xi - xj) / norm * (p - 1) - (norm - l) / (norm*norm)*(xi - xj)*(xi - xj)) + (lambda < 0)*abs(lambda)*(1 - x * x);
				row[i * 108 + j * 36 + 1] = yi_index; col[i * 108 + j * 36 + 1] = yi_index; value[i * 108 + j * 36 + 1] = p * ld / norm * (norm - l + (yi - yj)*(yi - yj) / norm * (p - 1) - (norm - l) / (norm*norm)*(yi - yj)*(yi - yj)) + +(lambda < 0)*abs(lambda)*(1 - y * y);
				row[i * 108 + j * 36 + 2] = zi_index; col[i * 108 + j * 36 + 2] = zi_index; value[i * 108 + j * 36 + 2] = p * ld / norm * (norm - l + (zi - zj)*(zi - zj) / norm * (p - 1) - (norm - l) / (norm*norm)*(zi - zj)*(zi - zj)) + +(lambda < 0)*abs(lambda)*(1 - z * z);
				row[i * 108 + j * 36 + 3] = yi_index; col[i * 108 + j * 36 + 3] = xi_index; value[i * 108 + j * 36 + 3] = p * (p - 1)*(xi - xj)*(yi - yj)*ld / (norm*norm) - p * (xi - xj)*(yi - yj)*ld*(norm - l) / (norm*norm*norm) - (lambda < 0)*abs(lambda)*x*y;
				row[i * 108 + j * 36 + 4] = zi_index; col[i * 108 + j * 36 + 4] = xi_index; value[i * 108 + j * 36 + 4] = p * (p - 1)*(xi - xj)*(zi - zj)*ld / (norm*norm) - p * (xi - xj)*(zi - zj)*ld*(norm - l) / (norm*norm*norm) - (lambda < 0)*abs(lambda)*x*z;
				row[i * 108 + j * 36 + 5] = zi_index; col[i * 108 + j * 36 + 5] = yi_index; value[i * 108 + j * 36 + 5] = p * (p - 1)*(zi - zj)*(yi - yj)*ld / (norm*norm) - p * (zi - zj)*(yi - yj)*ld*(norm - l) / (norm*norm*norm) - (lambda < 0)*abs(lambda)*z*y;
				memcpy(value + i * 108 + j * 36 + 6, value + i * 108 + j * 36 + 3, 3 * sizeof(double));
				row[i * 108 + j * 36 + 6] = xi_index; col[i * 108 + j * 36 + 6] = yi_index; //value[i * 108 + j * 36 + 6] = -2 * x * y*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 7] = xi_index; col[i * 108 + j * 36 + 7] = zi_index; //value[i * 108 + j * 36 + 9] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 8] = yi_index; col[i * 108 + j * 36 + 8] = zi_index;// value[i * 108 + j * 36 + 10] = -2 * z * y*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				memcpy(value + i * 108 + j * 36 + 9, value + i * 108 + j * 36, 9 * sizeof(double));
				row[i * 108 + j * 36 + 9] = xj_index; col[i * 108 + j * 36 + 9] = xj_index;//value[i * 108 + j * 36 + 1] = ld + 2*(ld < 0)*abs(ld)*(y*y + z*z) + 2 * pow(2 * xi - 2 * xj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*x*x;			
				row[i * 108 + j * 36 + 10] = yj_index; col[i * 108 + j * 36 + 10] = yj_index;//value[i * 108 + j * 36 + 3] = ld + 2 * (ld < 0)*abs(ld)*(x*x + z*z) + 2 * pow(2 * yi - 2 * yj, 2.0)+2 * (ld2 < 0)*abs(ld2)*y*y;	
				row[i * 108 + j * 36 + 11] = zj_index; col[i * 108 + j * 36 + 11] = zj_index;// value[i * 108 + j * 36 + 5] = ld + 2 * (ld < 0)*abs(ld)*(y*y + x*x) + 2 * pow(2 * zi - 2 * zj, 2.0) + 2 * (ld2 < 0)*abs(ld2)*z *z;
				row[i * 108 + j * 36 + 12] = yj_index; col[i * 108 + j * 36 + 12] = xj_index; //value[i * 108 + j * 36 + 31] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 13] = zj_index; col[i * 108 + j * 36 + 13] = xj_index; //value[i * 108 + j * 36 + 33] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;
				row[i * 108 + j * 36 + 14] = zj_index; col[i * 108 + j * 36 + 14] = yj_index; //value[i * 108 + j * 36 + 34] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;
				row[i * 108 + j * 36 + 15] = xj_index; col[i * 108 + j * 36 + 15] = yj_index; //value[i * 108 + j * 36 + 30] = -2 * y * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj) + 2 * (ld2 < 0)*abs(ld2)*x*y;
				row[i * 108 + j * 36 + 16] = xj_index; col[i * 108 + j * 36 + 16] = zj_index; //value[i * 108 + j * 36 + 32] = -2 * z * x*(ld < 0)*abs(ld) + 2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*x*z;			
				row[i * 108 + j * 36 + 17] = yj_index; col[i * 108 + j * 36 + 17] = zj_index; //value[i * 108 + j * 36 + 35] = -2 * y * z*(ld < 0)*abs(ld) + 2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj) + 2 * (ld2 < 0)*abs(ld2)*z*y;

				row[i * 108 + j * 36 + 18] = xi_index; col[i * 108 + j * 36 + 18] = xj_index; value[i * 108 + j * 36 + 18] = -value[i * 108 + j * 36];
				row[i * 108 + j * 36 + 19] = yi_index; col[i * 108 + j * 36 + 19] = yj_index; value[i * 108 + j * 36 + 19] = -value[i * 108 + j * 36 + 1];
				row[i * 108 + j * 36 + 20] = zi_index; col[i * 108 + j * 36 + 20] = zj_index; value[i * 108 + j * 36 + 20] = -value[i * 108 + j * 36 + 2];
				row[i * 108 + j * 36 + 21] = yi_index; col[i * 108 + j * 36 + 21] = xj_index; value[i * 108 + j * 36 + 21] = -value[i * 108 + j * 36 + 3];
				row[i * 108 + j * 36 + 22] = zi_index; col[i * 108 + j * 36 + 22] = xj_index; value[i * 108 + j * 36 + 22] = -value[i * 108 + j * 36 + 4];
				row[i * 108 + j * 36 + 23] = zi_index; col[i * 108 + j * 36 + 23] = yj_index; value[i * 108 + j * 36 + 23] = -value[i * 108 + j * 36 + 5];
				memcpy(value + i * 108 + j * 36 + 24, value + i * 108 + j * 36 + 21, 3 * sizeof(double));
				row[i * 108 + j * 36 + 24] = xi_index; col[i * 108 + j * 36 + 24] = yj_index; //value[i * 108 + j * 36 + 18] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 25] = xi_index; col[i * 108 + j * 36 + 25] = zj_index; //value[i * 108 + j * 36 + 21] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 26] = yi_index; col[i * 108 + j * 36 + 26] = zj_index; //value[i * 108 + j * 36 + 24] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);
				memcpy(value + i * 108 + j * 36 + 27, value + i * 108 + j * 36 + 18, 9 * sizeof(double));
				row[i * 108 + j * 36 + 27] = xj_index; col[i * 108 + j * 36 + 27] = xi_index; //value[i * 108 + j * 36 + 13] = -ld - 2 * pow(2 * xi - 2 * xj, 2);			
				row[i * 108 + j * 36 + 28] = yj_index; col[i * 108 + j * 36 + 28] = yi_index; //value[i * 108 + j * 36 + 15] = -ld - 2 * pow(2 * yi - 2 * yj, 2);			
				row[i * 108 + j * 36 + 29] = zj_index; col[i * 108 + j * 36 + 29] = zi_index; //value[i * 108 + j * 36 + 17] = -ld - 2 * pow(2 * zi - 2 * zj, 2);
				row[i * 108 + j * 36 + 30] = yj_index; col[i * 108 + j * 36 + 30] = xi_index; //value[i * 108 + j * 36 + 19] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);
				row[i * 108 + j * 36 + 31] = zj_index; col[i * 108 + j * 36 + 31] = xi_index; //value[i * 108 + j * 36 + 20] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);
				row[i * 108 + j * 36 + 32] = zj_index; col[i * 108 + j * 36 + 32] = yi_index; //value[i * 108 + j * 36 + 25] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);											
				row[i * 108 + j * 36 + 33] = xj_index; col[i * 108 + j * 36 + 33] = yi_index; //value[i * 108 + j * 36 + 23] = -2 * (2 * xi - 2 * xj)*(2 * yi - 2 * yj);						
				row[i * 108 + j * 36 + 34] = xj_index; col[i * 108 + j * 36 + 34] = zi_index; //value[i * 108 + j * 36 + 26] = -2 * (2 * xi - 2 * xj)*(2 * zi - 2 * zj);			
				row[i * 108 + j * 36 + 35] = yj_index; col[i * 108 + j * 36 + 35] = zi_index; //value[i * 108 + j * 36 + 28] = -2 * (2 * yi - 2 * yj)*(2 * zi - 2 * zj);

				for (int k = 0; k < 36; k++)
				{
					value[i * 108 + j * 36 + k] *= weight;
				}
			}
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
	else if (nrhs != 4 && nrhs != 5) {
	mexErrMsgTxt("Four input required.");
	}
	int energy_type = int(mxGetScalar(prhs[3]));
	if (energy_type != 1 && energy_type != 2 && energy_type != 3 && energy_type != 4 && energy_type != 5 && energy_type != 6)
	{
		mexErrMsgTxt("energ_type必须等于1,2,3,4,5或6");
	}
	if (energy_type != 6 && nrhs == 5)
	{
		mexErrMsgTxt("输入幂次p");
	}

	//获得矩阵的行数
	size_t point_number = mxGetM(prhs[0]);


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


	//获取指向输出参数的指针
	col = mxGetPr(plhs[0]);
	row = mxGetPr(plhs[1]);
	value = mxGetPr(plhs[2]);
	if (nrhs == 6)
	{
		double p = int(mxGetScalar(prhs[5]));
		hessian_matrix_triplet(face_number, point_number, points, faces, l_target, col, row, value, energy_type, p);
	}
	else
	{
		hessian_matrix_triplet(face_number, point_number, points, faces, l_target, col, row, value, energy_type);
	}
	

	Eigen::SparseMatrix < double >  A(3 * point_number, 3 * point_number);
	std::vector < Eigen::Triplet < double > > triplets;
	for (size_t i = 0; i < 108 * face_number; i++)
	{
		triplets.emplace_back(col[i] - 1, row[i] - 1, value[i]);
	}

	A.setFromTriplets(triplets.begin(), triplets.end());

	int *innerIndex = A.innerIndexPtr();
	double *valueptr = A.valuePtr();

	plhs[3] = mxCreateSparse(3 * point_number, 3 * point_number, A.nonZeros(), mxREAL);
	double *pr, *pi, *si, *sr;
	mwIndex *irs, *jcs;
	sr = mxGetPr(plhs[3]);

	irs = mxGetIr(plhs[3]);
	jcs = mxGetJc(plhs[3]);
	for (int i = 0; i < A.nonZeros(); i++)
	{
		sr[i] = valueptr[i];
	}

	for (int i = 0; i < A.nonZeros(); i++)
	{
		irs[i] = innerIndex[i];
	}

	int *OuterStarts = A.outerIndexPtr();
	for (int i = 0; i < 3 * point_number + 1; i++)
	{
		jcs[i] = OuterStarts[i];
	}
}


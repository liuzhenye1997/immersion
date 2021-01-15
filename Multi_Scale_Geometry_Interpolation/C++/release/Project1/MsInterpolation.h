#ifndef MS_INTERPOLATION__H
#define MS_INTERPOLATION__H

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
using namespace OpenMesh;

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Pair.h"

class Pair;
class PairVertex;

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef OpenMesh::VertexHandle MyVertexHandle;
typedef OpenMesh::FaceHandle MyFaceHandle;
typedef Eigen::Vector3d MyVector3f;
typedef std::vector<PairVertex> PointList;
typedef Eigen::MatrixXd MyMatrixXf;

/*** ϡ����� ***/
typedef Eigen::SparseMatrix<double> MySMatrixXf;

const unsigned int NR_MIN_PATCH_PER_LEVEL = 6;
const unsigned int LEAF_NODE_MIN_NR = 6;

/*** blending �Ƿ������� ***/
const bool BLENDING = true;

#define CHKROPT( Option ) \
	std::cout << "  provides " << #Option \
	<< (ropt.check(IO::Options:: Option)?": yes\n":": no\n")
#define CHKWOPT( Option ) \
	std::cout << "  write " << #Option \
	<< (wopt.check(IO::Options:: Option)?": yes\n":": no\n")
#define MESHOPT( msg, tf ) \
	std::cout << "  " << msg << ": " << ((tf)?"yes\n":"no\n")


class FaceNode
{
public:
	/*** ���index ***/
	std::vector<int> idx;

	/*** ��¼�ýڵ�ĵ���Ϣ��index������ֵ ***/
	std::vector<MyVector3f> pts;
	std::vector<int> pts_index;
	std::map<int, int> r_idx;

	/*** multiregistration ��һЩ��Ϣ ***/
	std::vector<Pair> P;
	std::vector<std::vector<int> > pl;// ��¼ÿ��Pair���еĵ��index
	

	/*** ������νṹ ***/
	FaceNode *next[ NR_MIN_PATCH_PER_LEVEL ];
	std::vector<int> boundray; // �͸���patch֮�乲�е�������index

	/*** blending������ Mv = e ***/
	MySMatrixXf M;

	/*** �����������index ***/
	std::vector<Pair> edges_vector[NR_MIN_PATCH_PER_LEVEL];


	FaceNode ()
	{
		for(int i = 0; i < NR_MIN_PATCH_PER_LEVEL; i++)
		{
			next[i] = NULL;
		}
	}
};


class MsInterpolation
{

public:
	MsInterpolation(void);
	~MsInterpolation(void);

	void read_mesh_data(std::string filename, std::string target_mesh);
	void write_mesh_data(std::string output);

	void test();
	void test(int n, char *output_path);
	void build_hierarchy_on_face();

private:
	void build_interpolation(FaceNode *subroot, MyMesh &src_mesh, MyMesh &target_mesh, double t);

	void build_registration_pair(FaceNode *subroot, int level = 0);
	void build_registration_pair();
	void build_hierarchy_based_on_face(FaceNode *subroot, int level = 0);
	void random_seed_point_helper(int, int *);

	void pre_blending_process(FaceNode * subroot);

	void release_facenode_help(FaceNode * &subroot);

private:
	//һЩ��������
	void max_min_height(FaceNode *subroot, int &min, int &max);

private:
	MyMesh m_mesh;     // source mesh
	MyMesh target_mesh;// target mesh

	/***���������νṹ***/
	FaceNode *face_root;
	int min_size_node;  // ״̬�����������鿴���ڵ��а��������������ж���
};


#endif //endif 
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MsInterpolation.h"

#include "LeafNodeInterpolation.h"

#include <iostream>
using namespace std;


int main(int argc, char *argv[])
{
	MsInterpolation ms;
	std::string src;
	std::string target;
	if (argc == 1)
	{
		src = "./bar_twisted_more_s.obj";
		target = "./bar_twisted_more_t.obj";
		ms.read_mesh_data(src, target);
		ms.test();
	}
	else
	{
		src = argv[1];
		target = argv[2];
		ms.read_mesh_data(src, target);
		int n = atoi(argv[3]);
		ms.test(n, argv[4]);
	}
	std::cout << src << std::endl;

	
	return 0;
}

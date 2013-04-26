// MixEnergParam.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include "Mesh.h"
using namespace std;
int _tmain(int argc, _TCHAR* argv[])
{
	Mesh mesh;
	string filename;
	cout<<"------------------------------------------"<<endl;
	cout<<"参数化 程序: MinEnergyParam.Cpp"<<endl;
	cout<<"By iTair , @ Mar 2013"<<endl;
	cout<<"------------------------------------------"<<endl<<endl;
	cout<<"输入网格文件名..." <<endl;
//	cin>>Mesh.filename;
	filename="nefertiti.off";
	cout<<"文件名: " <<filename<<endl<<endl;
	cout<<"------------------------------------------"<<endl;
	mesh.ReadMesh(filename);
	cout<<"网格读取完毕"<<endl;
	cout<<"------------------------------------------"<<endl;
	mesh.InitMesh();
	cout<<"网格初始化完毕-"<<endl;
	cout<<"------------------------------------------"<<endl;
	cout<<"开始参数化过程-"<<endl;
	mesh.RunFlatPara();
	cout<<"网格计算完毕"<<endl;
	cout<<"------------------------------------------"<<endl;
	mesh.MeshesOutput(filename);
	mesh.ResultOutput(filename);
	system("pause");
	return 0;
}


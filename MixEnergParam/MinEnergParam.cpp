// MixEnergParam.cpp : �������̨Ӧ�ó������ڵ㡣
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
	cout<<"������ ����: MinEnergyParam.Cpp"<<endl;
	cout<<"By iTair , @ Mar 2013"<<endl;
	cout<<"------------------------------------------"<<endl<<endl;
	cout<<"���������ļ���..." <<endl;
//	cin>>Mesh.filename;
	filename="nefertiti.off";
	cout<<"�ļ���: " <<filename<<endl<<endl;
	cout<<"------------------------------------------"<<endl;
	mesh.ReadMesh(filename);
	cout<<"�����ȡ���"<<endl;
		cout<<"------------------------------------------"<<endl;
	// TODO:
	//���������㷨 
	mesh.InitMesh();
	cout<<"�����ʼ�����-"<<endl;
		cout<<"------------------------------------------"<<endl;
	cout<<"��ʼ����������-"<<endl;
	mesh.RunFlatPara();

	//Mesh.RunFlatPara();

	cout<<"����������"<<endl;
	cout<<"------------------------------------------"<<endl;
	mesh.MeshesOutput(filename);
	mesh.ResultOutput(filename);

	system("pause");
	//getchar();
	return 0;
}


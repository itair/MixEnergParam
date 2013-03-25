//���������� ���� ��ȡ ������������
//
//
//
//
//
#pragma once
#include <string>
#include <vector>
#include <iostream>

using namespace std;

class Normal{
	double x,y,z;
};


class TriAngle{
public:
	unsigned long  v1,v2,v3;
	Normal normal;
	unsigned long e1,e2,e3;   //�����
	vector <unsigned long> adjFace;// �ڽ���Ƭ����. 
public:	
	void SetVextexIndex(int a,int b,int c){
		v1=a;		v2=b;		v3=c;
	}
	bool isEqual(TriAngle& tri){
// 		return((v1==tri.v1 && v2==tri.v2 && v3==tri.v3)||
// 		(v1==tri.v2 && v2==tri.v3 && v3==tri.v1)||
// 		(v1==tri.v3 && v2==tri.v1 && v3==tri.v3)||
// 		(v1==tri.v1 && v2==tri.v3 && v3==tri.v2)||
// 		(v1==tri.v2 && v2==tri.v1 && v3==tri.v3)||
// 		(v1==tri.v3 && v2==tri.v2 && v3==tri.v1));  //̫2��
		return ( (v1+v2+v3)==(tri.v1+tri.v2+tri.v3) && ((v1*v2*v3)==(tri.v1*tri.v2*tri.v3)) );
	}
};

class Edge{
public:
	unsigned long v1,v2;
	vector <unsigned long> adjFace;
public:
	void SetVextexIndex(int a,int b){
		v1=a;		v2=b;
	}
	bool isEqual(Edge& ed){
		if(v1==ed.v1 && v2==ed.v2)  return true;
		if (v1==ed.v2 && v2==ed.v1) return true;
		return false;
	}
};

class Vetex{
public:
	double x,y,z;
	vector <unsigned long> adjFace;
public:
	void SetVextxPos(double dx,double dy,double dz){
	x=dx;	y=dy;	z=dz;
	}
};

// class BorderEdge{
// public:
// 	unsigned long v1;   // ��1 ������ 
// 	unsigned long v2;   //��2 ������ 
// 	unsigned long faceIndex;    // ������Ƭ����  
// };

class PlanePara{
public:
	double u;
	double v;
public:
    void setParaPos(PlanePara &pl){
	   u=pl.u;	   v=pl.v;
   }
};

class TriMesh{
	//����
public:
	vector<TriAngle>	m_Mesh_faces;
	vector<Vetex>	m_Mesh_vetexs;
	vector<Edge> m_Mesh_edges;

	vector<TriAngle> m_Mesh_T1;//�Ѵ���
	vector<TriAngle> m_Mesh_T2;//δ����
	vector<TriAngle> m_Mesh_T0;//��ջ

	vector<Edge> m_Edge_T1;//�߳�
	vector<Edge> m_Edge_T2;//��ջ
	vector<Edge> m_Edge_Border; //�߽�

	int m_epsilon;

	Vetex	m_CurrentVex;
	TriAngle m_CurrentTri;
	Edge	m_CurrentEdge;

	vector<PlanePara>  m_Plane_Vertex;

	char  m_Mesh_format[5];
	string filename;

//����
public:
	TriMesh();
	~TriMesh();

	 void ReadMesh(string filename);
	 void ResultOutput(string filename);
	 void MeshesOutput(string filename);

	 void InitMesh(void);   //��ʼ������,�����������
	 void RunFlatPara(void);
		void GetNextVetex(void);
			void Flatten1stTir(void);
			void GetNextTri(void);
			void GetT1Boundary(void);
//				bool CompareEdgesinT2();
		void ComputeCurrEnery(void);	 
		void SloveMinEnery(void);
		void MapToSquare(void);
	 //void ShowMesh();
	 //void showParam();
	 

};
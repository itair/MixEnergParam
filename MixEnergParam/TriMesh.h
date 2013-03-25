//三角网格类 用于 读取 处理三角网格
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
	unsigned long e1,e2,e3;   //边序号
	vector <unsigned long> adjFace;// 邻接面片索引. 
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
// 		(v1==tri.v3 && v2==tri.v2 && v3==tri.v1));  //太2了
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
// 	unsigned long v1;   // 第1 点索引 
// 	unsigned long v2;   //第2 点索引 
// 	unsigned long faceIndex;    // 所属面片索引  
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
	//变量
public:
	vector<TriAngle>	m_Mesh_faces;
	vector<Vetex>	m_Mesh_vetexs;
	vector<Edge> m_Mesh_edges;

	vector<TriAngle> m_Mesh_T1;//已处理
	vector<TriAngle> m_Mesh_T2;//未处理
	vector<TriAngle> m_Mesh_T0;//面栈

	vector<Edge> m_Edge_T1;//边池
	vector<Edge> m_Edge_T2;//边栈
	vector<Edge> m_Edge_Border; //边界

	int m_epsilon;

	Vetex	m_CurrentVex;
	TriAngle m_CurrentTri;
	Edge	m_CurrentEdge;

	vector<PlanePara>  m_Plane_Vertex;

	char  m_Mesh_format[5];
	string filename;

//方法
public:
	TriMesh();
	~TriMesh();

	 void ReadMesh(string filename);
	 void ResultOutput(string filename);
	 void MeshesOutput(string filename);

	 void InitMesh(void);   //初始化网格,求得填满参数
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
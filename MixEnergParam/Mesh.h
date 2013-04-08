/*************************************************************
//  网格类 包含 参数化方法实现
//  
//
//
/
/
/
/ 修改自 TriMesh 类
***********************************************************/

#pragma once
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>

using namespace std;

typedef  unsigned long Index;

class sVector{
public:
	double x,y,z;
	sVector& operator- (const sVector &sv) {
		this->x=x-sv.x;
		this->y=y-sv.y;
		this->z=z-sv.z;
		return *this;
	}
	sVector& operator/ (const double real){
		this->x = x / real;
		this->y = y / real;
		this->z = z / real;
		return *this;
	}
	double operator* (sVector &sv){
		this->x = x * sv.x;
		this->y = y * sv.y;
		this->z = z * sv.z;
		return  (x*x + y*y+ z*z);
	}
};

class Vetex {
public:
	Index index;
	sVector pos;
	bool mark; //多用途标记....
	vector <Index> adjFace;
	vector <Index> adjEgde;
public:
	void SetPos(double dx,double dy,double dz){
		pos.x=dx; pos.y=dy;pos.z=dz;
	}
};

class Edge{
public:
	Index index;
	Index v1,v2;
	vector <Index> adjFace;
	bool mark; //多用途标记....
public:
 Edge(){	 mark=false;	 }
 void SetVexIndex(Index a,Index b){
		v1= ( a < b ) ? a : b; 
		v2= ( a >= b )? a : b;		
	}
	bool operator == (Edge edge){
		return (v1==edge.v1 && v2==edge.v2);
	}

};

class Face{
public:
	Index index;
	Index v1,v2,v3;
	Index e1,e2,e3;
	bool mark; //多用途标记....
	sVector normal;
	vector <Index> adjFace;
public:	
	void SetVextexIndex(Index a,Index b,Index c){
		v1=a;		v2=b;		v3=c;
	}
	void SetEdgeIndex(Index a,Index b,Index c){
		e1=a;		e2=b;		e3=c;
	}
	void SwapEdgeIndex(Index old,  Index newEdge ) {
		if (old == e1) { e1 = newEdge ; return; }
		if (old == e2) { e2 = newEdge ; return;}
		if (old == e3) { e3 = newEdge ; return; }
		
	}
	bool operator == (Face fa){
		return index==fa.index;
	}
	
};

class Border{
	Index v1;   // 第1 点索引 
	Index v2;   //第2 点索引 
	Index faceIndex;    // 所属面片索引  
};

class PlanePara{
public:
	Index index; //对应顶点的 索引
	double u;
	double v;
public:
	void setParaPos(PlanePara &pl){
		u=pl.u;	   v=pl.v;
	}
};

class Mesh{
	//变量
public:
	vector<Face>	m_Mesh_faces;
	vector<Vetex>	m_Mesh_vetexs;
	vector<Edge>    m_Mesh_edges;

	list<Index> m_Face_T1;	//已处理 //序号数组
	list<Index> m_Face_T2;	//未处理
	list<Index> m_Face_T0;	//面栈

	vector<Edge> m_Edge_T1;	//边池
	list<Index> m_Edge_T2;	//边栈
	list<Index> m_Edge_Border; //边界

	list<Index> m_Vertex_T1;	//点池
	list<Index> m_Vertex_T2;	//点栈
	list<Index> m_Vertex_Border; //边界点
	list<Index> m_Vertex_Free;	//自由点

	int m_epsilon;
	sVector m_CurrentPos;
	Vetex	m_CurrentVex;
	Face  m_CurrentTri;
	Edge	m_CurrentEdge;
	PlanePara m_CurrentPlane;

	sVector m_PlaneOrigin; //正交坐标系原点 对应的顶点坐标,firstTri设定v1
	sVector m_PlaneU;
	sVector m_PlaneV;
	vector<PlanePara>  m_Plane_Vertex;

	char  m_Mesh_format[5];
	string filename;


	//方法
public:
	Mesh();
	~Mesh();

	sVector NormalCross(sVector &v1,sVector &v2);
		sVector Cross (sVector &v1,sVector &v2);
		double More(sVector &v1);

	void ReadMesh(string filename);
	void ResultOutput(string filename);
	void MeshesOutput(string filename);

	void InitMesh(void);   // 初始化 网格,求得填满参数;
		void ProcessVetex(vector<Face>::iterator iter,Index index);
		void ProcessEdge(vector<Face>::iterator iter, Index index);
		void SimplyEdge(void);
		void ScanEdge(void);
		void ProcessFace();

	bool RunFlatPara(void);
		void GetFirstTri(void);
		bool GetNextVetex(void);
		void Flatten1stTir(void);
	
	void ComputeCurrEnergy(void);	 
	void SloveMinEnery(void);
	void MapToSquare(void);
	//void ShowMesh();
	//void showParam();


};
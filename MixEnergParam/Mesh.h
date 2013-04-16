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
#include <complex>
#include "BasicDS.h"
#include "TriMesh.h"
#include "PDerivative.h"
//#include "Complex.cpp"
using namespace std;


class Mesh{
	//变量
private:
	vector<Face>	m_Mesh_faces;
	vector<Vertex>	m_Mesh_vetexs;
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

	double  m_epsilon ;
	sVector m_CurrentPos;
	Vertex	m_CurrentVex;
	Face  m_CurrentTri;
	Face  m_OldTri;
	Edge	m_CurrentEdge;
	PlanePara m_CurrentPlane;

	sVector m_PlaneOrigin; //正交坐标系原点 对应的顶点坐标,firstTri设定v1
 	sVector m_PlaneU;  //正交基 X
 	 sVector m_PlaneV;	//正交基 Y
	 sVector m_PlaneNormal; //全局平面法向量
	vector<PlanePara>  m_Plane_T0; //T0中未优化的平面点 uv参数 集
	vector<PlanePara>  m_Plane_Vertex; //输出的自由边界 平面uv参数 点集 

	Energy E;
	Energy E1;
	Energy E2;

	char  m_Mesh_format[5];
	string filename;
	Index v123[3],v12[2],e123[3];
	double uv[2];
	PlanePara p123[3];

	//方法
public:
	Mesh();
	~Mesh();

	sVector NormalCross(const sVector &v1, const sVector &v2);
	//	double More(sVector &v1);

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
	    PDerivative PartialDerivative (Face Old, Face New);
			void FreeVertexProjection(void);
 
	void SloveMinEnery(void);
	void MapToSquare(void);
	// 自由定点在初始平面上的投影点


	
};
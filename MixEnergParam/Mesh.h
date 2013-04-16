/*************************************************************
//  ������ ���� ����������ʵ��
//  
//
//
/
/
/
/ �޸��� TriMesh ��
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
	//����
private:
	vector<Face>	m_Mesh_faces;
	vector<Vertex>	m_Mesh_vetexs;
	vector<Edge>    m_Mesh_edges;

	list<Index> m_Face_T1;	//�Ѵ��� //�������
	list<Index> m_Face_T2;	//δ����
	list<Index> m_Face_T0;	//��ջ

	vector<Edge> m_Edge_T1;	//�߳�
	list<Index> m_Edge_T2;	//��ջ
	list<Index> m_Edge_Border; //�߽�

	list<Index> m_Vertex_T1;	//���
	list<Index> m_Vertex_T2;	//��ջ
	list<Index> m_Vertex_Border; //�߽��
	list<Index> m_Vertex_Free;	//���ɵ�

	double  m_epsilon ;
	sVector m_CurrentPos;
	Vertex	m_CurrentVex;
	Face  m_CurrentTri;
	Face  m_OldTri;
	Edge	m_CurrentEdge;
	PlanePara m_CurrentPlane;

	sVector m_PlaneOrigin; //��������ϵԭ�� ��Ӧ�Ķ�������,firstTri�趨v1
 	sVector m_PlaneU;  //������ X
 	 sVector m_PlaneV;	//������ Y
	 sVector m_PlaneNormal; //ȫ��ƽ�淨����
	vector<PlanePara>  m_Plane_T0; //T0��δ�Ż���ƽ��� uv���� ��
	vector<PlanePara>  m_Plane_Vertex; //��������ɱ߽� ƽ��uv���� �㼯 

	Energy E;
	Energy E1;
	Energy E2;

	char  m_Mesh_format[5];
	string filename;
	Index v123[3],v12[2],e123[3];
	double uv[2];
	PlanePara p123[3];

	//����
public:
	Mesh();
	~Mesh();

	sVector NormalCross(const sVector &v1, const sVector &v2);
	//	double More(sVector &v1);

	void ReadMesh(string filename);
	void ResultOutput(string filename);
	void MeshesOutput(string filename);

	void InitMesh(void);   // ��ʼ�� ����,�����������;
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
	// ���ɶ����ڳ�ʼƽ���ϵ�ͶӰ��


	
};
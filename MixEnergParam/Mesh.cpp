// Mesh ����ʵ��

#include "StdAfx.h"
#include "Mesh.h"
#include <cmath>
using namespace std;

Mesh::Mesh(){
	filename.clear();
	m_Mesh_faces.clear();
	m_Mesh_edges.clear();
	m_Mesh_vetexs.clear();
	m_Edge_Border.clear();
	m_Edge_T1.clear();	
	m_epsilon = 0.5;
	
}

Mesh::~Mesh(){

}

void Mesh::ReadMesh(string filename){
	FILE *in=NULL;
	const char* file;
	file=filename.c_str();
	errno_t err;
	err  = fopen_s( &in, file, "r" );
	if( err == 0 ) cout<<"The file '"<<filename<<"' was opened"<<endl; 
	else cout<<"The file '"<<filename<<"' was not opened"<<endl; 

	char ext[5];
	int dV=0;
	int dF=0;
	int dE=0;
	int i=0;
	int temp(0); //for assert
	fscanf_s(in,"%s",ext,_countof(ext));
	fscanf_s(in,"%d",&dV);
	fscanf_s(in,"%d",&dF);
	fscanf_s(in,"%d",&dE);
	//cin>>ext>>dV>>dF>>dE;
	strcpy_s(m_Mesh_format,ext);
	cout<<"MeshFile Format: " <<ext<<endl;
	//memoryallocate(dV,dF);
	//read vertex
	Vertex dvetex;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;
	for(i=0;i<dV;i++){
		fscanf_s(in,"%lf %lf %lf",&dx,&dy,&dz);
		dvetex.SetPos(dx,dy,dz);
		m_Mesh_vetexs.push_back(dvetex);
	}
	temp=m_Mesh_vetexs.size();
	assert( temp == dV );
	cout<<"Vetex Num: "<<dV<<endl;
	//read faces
	int val=3;
	int di=0;
	int dj=0;
	int dk=0;
	Face dface;
	for(i=0;i<dF;i++){
		fscanf_s(in,"%d %d %d %d",&val,&di,&dj,&dk);
		dface.SetVextexIndex(di,dj,dk);
		m_Mesh_faces.push_back(dface);
	}
	assert(m_Mesh_faces.size()==dF);
	cout<<"Face Num: "<<dF<<endl;
	//Read Edgs:
	//cout<<endl<<"Edge Num: "<<dE;
	fclose(in);
	//cout<<"Mesh Read.\n";
	//cout<<"------------------------------------------"<<endl;
}

void Mesh::ResultOutput(string filename){
	FILE *in=NULL;
	const char* file;
	string file_result("result.txt");
	file=file_result.c_str();
	const char* sourcefile;
	sourcefile=filename.c_str();

	errno_t err;
	err  = fopen_s( &in, file, "w+" );
	if( err == 0 ) cout<<"The file '"<<file<<"' was opened"<<endl; 
	else cout<<"The file '"<<file<<"' was not opened"<<endl; 
//todo save���� ������������: ��������Ķ�����ƽ�����ϵ�����(2ά����) 
//  ��Ӧ�����ε����� ( ������ )
// ƽ��������꼯            vector<vetex>  m_Plane_Vertex
//������ ������             vector<Triangle> m_Mesh_faces
	//write �ļ�ͷ
	fseek(in,0L,SEEK_SET);	
	fprintf_s(in,"Result file : %s\n",file);
	fprintf_s(in,"Source file : %s\n",sourcefile);
	fprintf_s(in,"MeshFile format: ,%s\n",m_Mesh_format);
	//write ����
	fprintf_s(in,"Vetex Number: %d \n",m_Plane_Vertex.size());
	double uv[2];
	for (int i=0; i!=m_Plane_Vertex.size(); i++){
				m_Plane_Vertex.at(i).Getuv(uv);
			fprintf_s(in, "%lf\t%lf\n",uv);
	}
	//write faces 
	fprintf_s(in,"Faces Number: %d \n",m_Mesh_faces.size());
	Index v123[3];
	for (int i=0; i!=m_Mesh_faces.size(); i++){
		m_Mesh_faces.at(i).Getv123(v123)	;
		fprintf_s(in,"%d\t%d\t%d\n", v123);
	}	
	fclose(in);
	cout<<"Mesh Saved.\n";
	cout<<"------------------------------------------"<<endl;
	cout<<endl;
}

void Mesh::InitMesh(){
	//������߱�  ������ջ �ͱ�ջ
	m_Mesh_edges.clear();
	Index faceIndex=0;
	//Index vetexIndex=0;	
	for (vector<Face>::iterator iter = m_Mesh_faces.begin();
															iter != m_Mesh_faces.end();)	{
			//�����ڽӵ�ǰ��;
			(*iter).ResetMark();
			m_Face_T2.push_back(faceIndex);
			ProcessVetex(iter,faceIndex);
			ProcessEdge(iter,faceIndex);
			(*iter).SetIndex(faceIndex++) ;
			++iter;
	}		
	SimplyEdge();// �ϲ� �ظ�ͳ�Ƶı�/. ������ð��;
	ScanEdge(); // ��ȫ ����ٽ��� ����ڽӱ�;
	ProcessFace();		
}

void Mesh::ProcessVetex(vector<Face>::iterator iter, Index faceIndex){
		(*iter).Getv123(v123);
	for (int i =0 ; i != 3; i++)	{
		// �������� ���������� == ��Ƭ�е�������;
		m_Mesh_vetexs.at(v123[i]).SetIndex(v123[i]);
		// �ѵ�ǰ��������ż��� ����������ڽ����;
		m_Mesh_vetexs.at(v123[i]).AddadjFace(faceIndex);
	}
}

void Mesh::ProcessEdge(vector<Face>::iterator iter,Index faceIndex){
	// ���ι��������� ѹ��߱� ��������;
	//Index *v123,*v12,*e123;
	(*iter).Getv123(v123);
	for (int i= 0; i != 3; i++){
		v12[0]=v123[i];
		v12[1]=v123[((i+1)%3)];
		m_CurrentEdge.SetVexIndex(v12[0],v12[1]); // �ں��������;
		Index edgeIndex = faceIndex + i;
		m_CurrentEdge.SetIndex(edgeIndex);
		m_CurrentEdge.ClearAdjFace();
		m_CurrentEdge.AddadjFace(faceIndex);
		m_Edge_T1.push_back(m_CurrentEdge);
		e123[i] = edgeIndex;
	}
	(*iter).Sete123(e123);
}

void Mesh::SimplyEdge(){
	//ȥ���ظ��ı� ���±��;
	Index face1, face2,oldEdge1,oldEdge2;
	Index edgeIndex=0;
	list<Edge>::iterator iter = m_Edge_T1.begin();		
	list<Edge>::iterator iter2; 
	for (list<Edge>::iterator iter= m_Edge_T1.begin();iter!= m_Edge_T1.end();){		
		for (list<Edge>::iterator iter2 = iter ; iter2 != m_Edge_T1.end(); ++iter2){
			if (*iter2 == *iter){
				face2 = iter2->GetLastAdjFace();
				oldEdge2 = iter2->GetIndex();
				m_Mesh_faces.at(face2).SwapEdgeIndex(oldEdge2, edgeIndex);
				iter->AddadjFace(face2);
				m_Edge_T1.erase(iter2);			
			}
		}
		face1 = iter->GetLastAdjFace();
		oldEdge1 = iter->GetIndex(); 
		m_Mesh_faces.at(face1).SwapEdgeIndex(oldEdge1, edgeIndex);
		if ((iter->GetAdjFaceSize()) > 2 ){
			cout<<"Edge Num: "<<iter->GetIndex()<<endl;
			throw	runtime_error("adjFace.size()>2!");
		}	
		m_CurrentEdge = *iter;
		m_Mesh_edges.push_back(m_CurrentEdge);		
		edgeIndex++;
		iter++;
	}	
}
void Mesh::ScanEdge(){
	//ɨ��߱�, ���䶥����ڽӱ� �� �ߵ��ڽ���;
	//Index *v12;
	Index edgeIndex, faceIndex1, faceIndex2 ;
	for (vector<Edge>::iterator iter=m_Mesh_edges.begin();iter != m_Mesh_edges.end();){
		iter->Getv12(v12);
		edgeIndex = iter->GetIndex();
		m_Mesh_vetexs.at(v12[0]).AddadjEdge(edgeIndex);
		m_Mesh_vetexs.at(v12[0]).AddadjEdge(edgeIndex);
		faceIndex1 = iter->GetFirstAdjFace();
		faceIndex2 = iter->GetLastAdjFace();
		m_Mesh_faces.at(faceIndex1).AddadjFace(faceIndex2);
		m_Mesh_faces.at(faceIndex2).AddadjFace(faceIndex1);
		++iter;
	}		
}
 
 sVector Mesh::NormalCross( const sVector &v1, const sVector &v2){   //�淶�����		
	sVector norm;
	double more; 	
	norm. Cross( v1 , v2 );
	more = norm.More();
	norm = norm / more;
	return  norm;	
}
void Mesh::ProcessFace(){
	//���� ����� ÿ�������ε� ������ �� �ֲ������� ƽ������ ���; �������.
	sVector a,b,c,A,B,X,Y,N,XB;
	double more,Area,p2u,p2v,p3u,p3v;
	//Index *v123;	
	PlanePara p123[3];
	for (vector<Face>::iterator iter=m_Mesh_faces.begin(); 
															iter !=m_Mesh_faces.end();)	{
		iter->Getv123(v123); 		
		a = m_Mesh_vetexs.at(v123[0]).GetPos();
		b = m_Mesh_vetexs.at(v123[1]).GetPos();
		c = m_Mesh_vetexs.at(v123[2]).GetPos();
 		//��һ��������		
		A =  b - a ;
		B =  c - a ;
		N = NormalCross( A , B);
		(*iter).SetNormal( N );
		X = NormalCross( N, B ); //u����		
		XB.Cross( X, B );
		more = XB.More();
		Y = XB / more;   //v����
		p2u = A*X;
		p2v = A*Y;
		p3u = B*X ;
		p3v =  B*Y;

		m_CurrentPlane.Setuv(0.0 , 0.0);
		p123[0]=m_CurrentPlane;
		m_CurrentPlane.Setuv( p2u , p2v);
		p123[1]=m_CurrentPlane;
		m_CurrentPlane.Setuv( p3u , p3v);
		p123[2]=m_CurrentPlane;
		iter->Getp123(p123);

		 Area = 0.5 * ( (p2u * p3v) - (p2v * p3u));
		 (*iter).SetArea (fabs(Area));
		//assert(  (*iter).Area<=0 );
		++iter;
	}
}

void Mesh::MeshesOutput(string filename) {
	FILE *in=NULL;
	const char* file;
	string file_result("MeshesOutput.txt");
	file=file_result.c_str();
	const char* sourcefile;
	sourcefile=filename.c_str();

	errno_t err;
	err  = fopen_s( &in, file, "w+" );
	if( err == 0 ) cout<<"The file '"<<file<<"' was opened"<<endl; 
	else cout<<"The file '"<<file<<"' was not opened"<<endl; 
	//TODO ���� ��ʼ������������ݼ�;
	// ƽ��������꼯              vector<vetex>  m_Plane_Vertex
	//������ ������              vector<Triangle> m_Mesh_faces
	//write �ļ�ͷ
	fseek(in,0L,SEEK_SET);	
	fprintf_s(in,"Meshes file : %s\n",file);
	fprintf_s(in,"Source file : %s\n",sourcefile);
	fprintf_s(in,"MeshFile format: %s\n",m_Mesh_format);
	//write ����
	fprintf_s(in,"Vetex Number: %d \n",m_Mesh_vetexs.size());
	for (int i=0; i!=m_Mesh_vetexs.size(); i++){
		m_CurrentVex = m_Mesh_vetexs.at(i);
		double xyz[3];
		m_CurrentVex.GetXYZ(xyz);
		fprintf_s(in, "%lf\t%lf\t%lf\n", xyz[0], xyz[1], xyz[2] );
	}
	//write faces 
	fprintf_s(in,"Faces Number: %d \n",m_Mesh_faces.size());
  for (int i=0; i!=m_Mesh_faces.size(); i++){
		m_CurrentTri = m_Mesh_faces.at(i);
		Index v123[3];
		m_CurrentTri.Gete123(v123);
		fprintf_s(in, "%ld\t%ld\t%ld\n", v123[0],v123[1],v123[2]);
	}
	//write edges
	fprintf_s(in,"Edges Number: %d \n",m_Mesh_edges.size());
	for (int i=0; i!=m_Mesh_edges.size(); i++)	{
		m_CurrentEdge=m_Mesh_edges.at(i);
		Index v12[2];
		m_CurrentEdge.Getv12(v12);
		fprintf_s(in, "%ld\t%ld\n",v12[0],v12[1]);
	}
	fclose(in);
	cout<<"Mesh Saved.\n";
	cout<<"------------------------------------------"<<endl;
	cout<<endl;
}

bool Mesh::RunFlatPara(){
	m_Face_T1.clear();
	m_Face_T0.clear();
	m_Plane_Vertex.clear();
	bool check;
	//ѡȡ��һ��������
	GetFirstTri();
	Flatten1stTir();
	//�㷨��ѭ��
	while (m_Face_T2.empty()!=true)	{		
		check = GetNextVetex();
		if (check==false && m_Face_T2.empty()!=true) {
			cout<<"�����ж�, �޷��ҵ���һ��T0... " <<endl;
			return false;
		  }	
		FreeVertexProjection();   //�������ɵ� �� ��ʼƽ����ϵ�ͶӰ ��ת��Ϊƽ������ uv ���� m_Plane_T0 
		ComputeCurrEnergy()  ;   // ���� T0 ������ʾ;
//TODO ��С������ ��� m_Vertex_Free ��Ӧ 2D���� ���� m_Plane_Vertex
//		SloveMinEnery();		
		for (list<Index>::iterator iter = m_Face_T0.begin();
				iter != m_Face_T0.end();	iter ++){  //���� T1 T2;
			m_Face_T2.remove(*iter);
			m_Face_T1.push_back(*iter);
		}
}
	return true;
}
void Mesh::GetFirstTri(){
	//���õ�һ��������
	Index firstface=m_Mesh_faces.size()/2;
	m_Face_T0.push_back(firstface);
	m_Face_T1.push_back(firstface);
	m_Face_T2.remove(firstface);
	m_Mesh_faces.at(firstface).IsMarked();
}

void Mesh::Flatten1stTir(){
	// Դ������ ������2Dƽ��;
	Index index_face,index_vertex;	
	Index v123[3];	
	PlanePara p123[3];
	index_face=m_Face_T0.back();
	m_CurrentTri = m_Mesh_faces.at(index_face);
	m_CurrentTri.Getv123(v123);
	m_CurrentTri.Getp123(p123);	
 	m_PlaneNormal = m_CurrentTri.GetNormal(); //Ƕ��ƽ���ͳһ������
	for (int i = 0; i< 2 ; i++)	{
		index_vertex = v123[i];
		m_CurrentPos = m_Mesh_vetexs.at(index_vertex).GetPos();
		if (i=0) { m_PlaneOrigin = m_CurrentPos; }    //2Dԭ��
		m_CurrentPlane = p123[i];
		m_CurrentPlane.SetIndex( index_vertex );
		m_Plane_Vertex.push_back(m_CurrentPlane);	
	}	
}

bool Mesh::GetNextVetex(){
	//T0�����б��� ����δmark�ı� ��Ӧ�Ķ��㼯 m_Vertex_Border;
	//m_Vertex_Borderi��ȫ������δmark�� Ϊ��T0;
	//��T0�� ����m_Vertex_Border�еĶ���Ϊ���ɶ���;
	Index index_face;
	Index e123[3];
	for (list<Index>::iterator iter = m_Face_T0.begin();
				iter != m_Face_T0.end();	iter ++){
		m_CurrentTri = m_Mesh_faces.at(*iter); 
		m_CurrentTri.Gete123(e123);
		for (int i=0; i<2; i++) {
			m_Edge_Border.push_back(e123[i]);
		}
	}

	m_Edge_Border.sort();
	m_Edge_Border.unique();
	m_Vertex_Border.clear();// 

	for (list<Index>::iterator iter = m_Edge_Border.begin();
				iter != m_Edge_Border.end();	iter ++){
		m_CurrentEdge = m_Mesh_edges.at(*iter);
		while (!(m_CurrentEdge.IsAdjFaceEmpty())) {
			index_face = m_CurrentEdge.GetLastAdjFace();
			m_CurrentEdge.DeladjFace();
			if (!(m_Mesh_faces.at(index_face).IsMarked())){
				Index v12[2];
				m_CurrentEdge.Getv12(v12);
				m_Vertex_Border.push_back(v12[0]);
				m_Vertex_Border.push_back(v12[1]);
			}
		}
	}
	m_Vertex_Border.sort();
	m_Vertex_Border.unique();
	m_Face_T0.clear();
	m_Vertex_Free.clear();//��ջ
	for (list<Index>::iterator iter = m_Vertex_Border.begin();
				iter != m_Vertex_Border.end();	iter ++)	{
		m_CurrentVex = m_Mesh_vetexs.at(*iter);
		while (!(m_CurrentVex.IsAdjFaceEmpty()) ){
			index_face = m_CurrentVex.GetLastAdjFace();;
			m_CurrentVex.DeladjFace();
			if (!(m_Mesh_faces.at(index_face).IsMarked())){
				m_Mesh_faces.at(index_face).GetMarked();
				m_Face_T0.push_back(index_face);//��T0�ҵ�
				Index v123[3];
				m_Mesh_faces.at(index_face).Getv123(v123);
				for (int i=0; i <2 ; i++){
					m_Vertex_Free.push_back(v123[i]);
				}				
			}//if
		}		
	}//for
	m_Vertex_Free.sort();
	m_Vertex_Free.unique();
	return !m_Face_T0.empty();
}

void Mesh::ComputeCurrEnergy(){
	//����T0��������   m_Face_T0 ,m_Plane_T0	
	E=0;
	for ( list<Index>::iterator iter = m_Face_T0.begin();
															iter != m_Face_T0.end(); iter++){	
		m_OldTri = m_Mesh_faces.at(*iter);
		m_CurrentTri=m_OldTri;
		Index v123[3];
		PlanePara p123[3];
		m_CurrentTri.Getv123(v123);
		m_CurrentTri.Getp123(p123);
		Index AimIndex;
		for (vector<PlanePara>::iterator iter2 = m_Plane_T0.begin();
							iter2 != m_Plane_T0.end(); iter2++){
								AimIndex = iter2->GetIndex();
								for (int i=0; i < 2 ; i++){
									if (AimIndex == v123[i] ){
										p123[i] = (*iter2);
									}
								}
								m_CurrentTri.Setp123(p123);// 
				
		}
		PDerivative pl; //ԭʼ�������������������
		Complex fu,fv;
		Complex f,I;
		Energy detI;
		pl = PartialDerivative ( m_OldTri, m_CurrentTri);
		Energy EE=0.0;		
		E1=0.0;
		E2=0.0;
		fu = pl.Getfu();
		fv = pl.Getfv();
			I = Complex (0,1);
			f = fu + I*fv ;
			E1 = abs(f);
			E1 = E1 *E1;
			detI = (pl.DotProduct(fu, fu) * pl.DotProduct(fv, fv)) - 
						 (pl.DotProduct(fu, fv) * pl.DotProduct(fu, fv));
			E2 = detI -1 ;
			E2 = E2 * E2;
			EE = m_epsilon * E1 + (1 -m_epsilon ) *E2;
			E = E + EE;
	}
	//ȡ�� E  ��ʼֵ ��Ӧ��ǰT0
}
// ���ɶ����ڳ�ʼƽ���ϵ�ͶӰ��
void Mesh::FreeVertexProjection(void){
	//m_Vertex_Free ; m_Vertex_Border; m_PlaneNormal; m_Face_T0;
	sVector footpoint; //���� 
	double t;//parameter_t
 m_Plane_T0.clear();
	for (list<Index>::iterator iter = m_Vertex_Free.begin();
											iter != m_Vertex_Free.end(); iter++)	{
				m_CurrentPos = (m_Mesh_vetexs.at(*iter )).GetPos();
				footpoint = m_CurrentPos - m_PlaneOrigin ;
					t = footpoint * m_PlaneNormal;
				footpoint = m_CurrentPos - ( m_PlaneNormal * t);
				m_CurrentPlane.SetIndex( *iter);
				double uu = footpoint * m_PlaneU ;
				double vv = footpoint * m_PlaneV ;
				m_CurrentPlane.Setuv(uu , vv );
				m_Plane_T0.push_back(m_CurrentPlane)	;
	}		
}
PDerivative Mesh::PartialDerivative(Face OldTri, Face NewTri){
	//��ƫ��
	PDerivative pl;
	Complex fu,fv;
	Complex p[3];
	double u[3],v[3];
	double Area = OldTri.GetArea();
	PlanePara p123[3];
	double uv[2];
	
	NewTri.Getp123(p123);
	for (int i=0; i<2; i++)	{
		p123[i].Getuv(uv);
		p[i] = Complex (uv[0], uv[1]);
	}

	OldTri.Getp123(p123);
	for (int i=0; i<2; i++)	{
		p123[i].Getuv(uv);
		u[i] = uv[0];
		v[i] = uv[1] ;
	}
	fu = ((v[1]-v[2])*p[0]) + ((v[2]-v[0])*p[1]) + ((v[0]-v[1])*p[2]);
	fu = fu / (2*Area);
	pl.Setfu(fu);
	fv = ((u[1]-u[2])*p[0]) + ((u[2]-u[0])*p[1]) + ((u[0]-u[1])*p[2]);
	fv = fv / (2*Area);
	pl.Setfv(fv);
	return pl;
}

void Mesh::SloveMinEnery(){
	//�ݶ��½��� ����С������;
	//����: ��ʼ���� E , Ŀ�꺯�� E(Tx), �������� a , �������� e ;
	//��Ȼ����֪�� �ݶȱ��ʽ ;��� �ݶȹ�С���µ�������̫�� ����ʹ�ýض��ݶȼ����������ݶȼ���


}
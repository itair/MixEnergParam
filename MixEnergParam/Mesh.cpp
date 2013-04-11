// Mesh ����ʵ��

#include "StdAfx.h"
#include "Mesh.h"
using namespace std;

Mesh::Mesh(){
	filename.clear();
	m_Mesh_faces.clear();
	m_Mesh_edges.clear();
	m_Mesh_vetexs.clear();
	m_Edge_Border.clear();
	m_Edge_T1.clear();	
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
	Vetex dvetex;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;
	for(i=0;i<dV;i++){
		fscanf_s(in,"%lf %lf %lf",&dx,&dy,&dz);
		dvetex.SetPos(dx,dy,dz);
		m_Mesh_vetexs.push_back(dvetex);
	}
	temp=m_Mesh_vetexs.size();
	assert(temp==dV);
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
	PlanePara vetex2d;
	for (int i=0; i!=m_Plane_Vertex.size(); i++){
		vetex2d.setParaPos(m_Plane_Vertex[i]);
		fprintf_s(in,"%lf\t%lf\n",vetex2d.u,vetex2d.v);
	}
	//write faces 
	fprintf_s(in,"Faces Number: %d \n",m_Mesh_faces.size());
	Face face2d;
	for (int i=0; i!=m_Mesh_faces.size(); i++){
		face2d=m_Mesh_faces[i];
		fprintf_s(in,"%d\t%d\t%d\n",face2d.v1,face2d.v2,face2d.v3);
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
	for (vector<Face>::iterator iter=m_Mesh_faces.begin();	iter!=m_Mesh_faces.end();)	{
			//�����ڽӵ�ǰ��;
			iter->mark=false;
			m_Face_T2.push_back(faceIndex);
			ProcessVetex(iter,faceIndex);
			ProcessEdge(iter,faceIndex);
			iter->index=faceIndex++ ;
			++iter;
	}		
	SimplyEdge();//�ϲ� �ظ�ͳ�Ƶı�/. ������ð��;
	ScanEdge(); //��ȫ ����ٽ��� ����ڽӱ�;
	ProcessFace();
		
}

void Mesh::ProcessVetex(vector<Face>::iterator iter, Index faceIndex){
	//�������� ���������� == ��Ƭ�е�������
	m_Mesh_vetexs.at(iter->v1).index=iter->v1;
	m_Mesh_vetexs.at(iter->v2).index=iter->v2;
	m_Mesh_vetexs.at(iter->v3).index=iter->v3;
	//�ѵ�ǰ��������� ���� ����������ڽ����;	
	m_Mesh_vetexs.at(iter->v1).adjFace.push_back(faceIndex);
	m_Mesh_vetexs.at(iter->v2).adjFace.push_back(faceIndex);
	m_Mesh_vetexs.at(iter->v3).adjFace.push_back(faceIndex);
}

void Mesh::ProcessEdge(vector<Face>::iterator iter,Index faceIndex){
		m_CurrentEdge.adjFace.clear();
	Index edgeIndex=0;
	//���ι��������� ѹ��߱� ��������;
	m_CurrentEdge.SetVexIndex(iter->v1,iter->v2); //�ں��������;
	m_CurrentEdge.index=edgeIndex;
	iter->e1=edgeIndex++;
	m_CurrentEdge.adjFace.clear();
	m_CurrentEdge.adjFace.push_back(faceIndex);
	m_Edge_T1.push_back(m_CurrentEdge);

	m_CurrentEdge.SetVexIndex(iter->v2,iter->v3); //�ں��������;
	m_CurrentEdge.index=edgeIndex;
	iter->e2=edgeIndex++;
	m_CurrentEdge.adjFace.clear();
	m_CurrentEdge.adjFace.push_back(faceIndex);
	m_Edge_T1.push_back(m_CurrentEdge);

	m_CurrentEdge.SetVexIndex(iter->v1,iter->v3); //�ں��������;
	m_CurrentEdge.index=edgeIndex;
	iter->e3=edgeIndex++;
	m_CurrentEdge.adjFace.clear();
	m_CurrentEdge.adjFace.push_back(faceIndex);
	m_Edge_T1.push_back(m_CurrentEdge);
}

void Mesh::SimplyEdge(){
	//ȥ���ظ��ı� ���±��;
	Index face1, face2,oldedgeIndex1,oldedgeIndex2;
	Index edgeIndex=0;
	vector<Edge>::iterator iter = m_Edge_T1.begin();		
	vector<Edge>::iterator iter2 = m_Edge_T1.begin();
	++iter2;	
	while(iter != m_Edge_T1.end()){
		if (iter->mark==true)  {
			++iter;
			iter2=iter;
			continue;//�����ظ���;
		}
		++iter2;
		while( iter2!=m_Edge_T1.end()){
			if ( *iter2 == *iter ){	//�ظ���	;	
				face2=iter2->adjFace.back();		
				oldedgeIndex2=iter2->index;//ԭ�߱����� ����	;	
				m_Mesh_faces.at(face2).SwapEdgeIndex(oldedgeIndex2,edgeIndex);//���¶�Ӧ�ڽ��� �ı߱�;
				iter->adjFace.push_back(face2);
				iter2->mark=true;
				break;
			}else{
				++iter2;
			}
		}// iter2
		face1=iter->adjFace.front();
		oldedgeIndex1=iter->index;
		m_Mesh_faces.at(face1).SwapEdgeIndex(oldedgeIndex1,edgeIndex);
		m_CurrentEdge.adjFace.clear();
		m_CurrentEdge=*iter;	

		if (iter->adjFace.size()>2){
			cout<<"Edge Num: "<<iter->index<<endl;
			throw	runtime_error("adjFace.size()>2!");
		}	

		m_CurrentEdge.index=edgeIndex++;
		m_Mesh_edges.push_back(m_CurrentEdge);
		++iter;
		iter2=iter;		    
	}//iter
}
void Mesh::ScanEdge(){
	//ɨ��߱�, ���䶥����ڽӱ� �� �ߵ��ڽ���;
	for (vector<Edge>::iterator iter=m_Mesh_edges.begin();iter != m_Mesh_edges.end();){
		m_Mesh_vetexs.at(iter->v1).adjEgde.push_back(iter->index);
		m_Mesh_vetexs.at(iter->v2).adjEgde.push_back(iter->index);

		m_Mesh_faces.at(iter->adjFace.front()).adjFace.push_back(iter->adjFace.back());
		m_Mesh_faces.at(iter->adjFace.back()).adjFace.push_back(iter->adjFace.front());

		++iter;
	}
		
}
sVector Mesh::NormalCross( sVector &v1, sVector&v2 ){
	//�淶�����
	sVector norm;
	double more; 	
	norm = Cross( v1 , v2 );
	more = More( norm );
	norm =norm / more;
	return  norm;	
}

sVector Mesh::Cross(sVector &v1,sVector &v2){
	sVector ab,bc,norm;
	ab = v1;
	bc = v2;
	norm.x = (ab.y * bc.z) - (ab.z * bc.y);
	norm.y = -((ab.x * bc.z) - (ab.z * bc.x));
	norm.z = (ab.x * bc.y) - (ab.y * bc.x);
	return norm;
}

double Mesh::More(sVector &norm){
	return sqrt(norm.x *norm.x + norm.y*norm.y + norm.z*norm.z);
}

void Mesh::ProcessFace(){
	//���� ����� ÿ�������ε� ������ �� �ֲ������� ƽ������ ���;
	sVector a,b,c,A,B,X,Y,N,XB;
	double more;
	for (vector<Face>::iterator iter=m_Mesh_faces.begin(); iter !=m_Mesh_faces.end();)	{
		a = m_Mesh_vetexs.at(iter->v1).pos;
		b = m_Mesh_vetexs.at(iter->v2).pos;
		c = m_Mesh_vetexs.at(iter->v3).pos;
 		//��һ��������		
		A =  b - a ;
		B =  c - a ;
		N = NormalCross( A , B);
		(*iter).normal = N ;
		X = NormalCross( N, B ); //u����
		Y = Cross( N, B );
		XB= Cross( X, B );
		more=More(XB);
		Y = Y / more;   //v����
		//p1 = (0,0);
		(*iter).p1.u =0.0;
		(*iter). p1.v=0.0;
		// p2= ( A.X , A.Y);
		(*iter).p2.u = A*X ;	
		(*iter).p2.v = A*Y;
		//p3= ( B.X , B.Y);
		(*iter).p3.u = B * X ;	
		(*iter).p3.v = B * Y ;
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
	//PlanePara vetex2d;
	for (int i=0; i!=m_Mesh_vetexs.size(); i++){
		m_CurrentVex=m_Mesh_vetexs[i];
		fprintf_s(in, "%lf\t%lf\t%lf\n",
			m_CurrentVex.pos.x, m_CurrentVex.pos.y, m_CurrentVex.pos.z, m_CurrentVex.adjFace);
	}
	//write faces 
	fprintf_s(in,"Faces Number: %d \n",m_Mesh_faces.size());
	//TriAngle face2d;
	for (int i=0; i!=m_Mesh_faces.size(); i++){
		m_CurrentTri=m_Mesh_faces[i];
		fprintf_s(in, "%ld\t%ld\t%ld\n", m_CurrentTri.v1, m_CurrentTri.v2, m_CurrentTri.v3);
	}
	//write edges
	fprintf_s(in,"Edges Number: %d \n",m_Mesh_edges.size());
	for (int i=0; i!=m_Mesh_edges.size(); i++)	{
		m_CurrentEdge=m_Mesh_edges[i];
		fprintf_s(in, "%ld\t%ld\n", m_CurrentEdge.v1, m_CurrentEdge.v2, m_CurrentEdge.adjFace);
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
		//list<Index> m_Face_T0;	// ��ǰ ��ջ
		//  list<Index> m_Vertex_Free;	//��ǰ�� ���ɵ�		
		FreeVertexProjection();//�������ɵ� �� ��ʼƽ����ϵ�ͶӰ ��ת��Ϊƽ������ uv ���� m_Plane_T0 
//TODO ���� T0 ������ʾ;
//		ComputeCurrEnery();
//TODO ��С������ ��� m_Vertex_Free ��Ӧ 2D���� ���� m_Plane_Vertex
//		SloveMinEnery();

		//���� T1 T2;
		for (list<Index>::iterator iter = m_Face_T0.begin();
				iter != m_Face_T0.end();	iter ++){
			m_Face_T2.remove(*iter);
			m_Face_T1.push_back(*iter);
		}
		//m_Face_T1.merge(m_Face_T0);
	}
	return true;
}
void Mesh::GetFirstTri(){
	//���õ�һ��������
	Index firstface=m_Mesh_faces.size()/2;
	m_Face_T0.push_back(firstface);
	m_Face_T1.push_back(firstface);
	m_Face_T2.remove(firstface);
	m_Mesh_faces.at(firstface).mark=true;
}

void Mesh::Flatten1stTir(){
	// Դ������ ������2Dƽ��;
	Index index_face,index_vertex;
	sVector a,b,c,A,B,X,Y,N,XB;
	double more;
	//��һ��������
	index_face=m_Face_T0.back();
	m_CurrentTri=m_Mesh_faces.at(index_face);
	a = m_Mesh_vetexs.at(m_CurrentTri.v1).pos;
	b = m_Mesh_vetexs.at(m_CurrentTri.v2).pos;
	c = m_Mesh_vetexs.at(m_CurrentTri.v3).pos;
	A =  b - a ;
	B =  c - a ;
	N = m_CurrentTri.normal;
	m_PlaneNormal=N; //Ƕ��ƽ���ͳһ������
	X = NormalCross( N, B ); //u����
	m_PlaneU=X;
	Y = Cross( N, B );
	XB= Cross( X, B );
	more=More(XB);
	Y = Y / more;   //v����
	m_PlaneV=Y;
	// Ƕ�� 2d ƽ��;
	//v1 = (0,0);
	index_vertex = m_CurrentTri.v1;
	m_CurrentPos = m_Mesh_vetexs.at(index_vertex).pos;
	m_PlaneOrigin = m_CurrentPos;    //2Dԭ��

	m_CurrentPlane.index=index_vertex;
	m_CurrentPlane.u = 0.0 ;
	m_CurrentPlane.v = 0.0 ;
	m_Plane_Vertex.push_back(m_CurrentPlane);
	
	// v2= ( A.X , A.Y);
	index_vertex = m_CurrentTri.v2;
	m_CurrentPlane.index=index_vertex;
	m_CurrentPlane.u = A*X ;
	m_CurrentPlane.v = A*Y;
	m_Plane_Vertex.push_back(m_CurrentPlane);
	//v3= ( B.X , B.Y);
	index_vertex = m_CurrentTri.v3;
	m_CurrentPlane.index=index_vertex;
	m_CurrentPlane.u = B * X ;
	m_CurrentPlane.v = B * Y;
	m_Plane_Vertex.push_back(m_CurrentPlane);

}

bool Mesh::GetNextVetex(){
	//T0�����б��� ����δmark�ı� ��Ӧ�Ķ��㼯 m_Vertex_Border;
	//m_Vertex_Borderi��ȫ������δmark�� Ϊ��T0;
	//��T0�� ����m_Vertex_Border�еĶ���Ϊ���ɶ���;
	Index index_face;

	for (list<Index>::iterator iter = m_Face_T0.begin();
				iter != m_Face_T0.end();	iter ++){
		m_CurrentTri = m_Mesh_faces.at(*iter); 
		m_Edge_Border.push_back(m_CurrentTri.e1);
		m_Edge_Border.push_back(m_CurrentTri.e2);
		m_Edge_Border.push_back(m_CurrentTri.e3);
	}

	m_Edge_Border.sort();
	m_Edge_Border.unique();
	m_Vertex_Border.clear();// 
	for (list<Index>::iterator iter = m_Edge_Border.begin();
				iter != m_Edge_Border.end();	iter ++){
		m_CurrentEdge = m_Mesh_edges.at(*iter);
		while (m_CurrentEdge.adjFace.empty() != true){
			index_face = m_CurrentEdge.adjFace.back();
			m_CurrentEdge.adjFace.pop_back();
			if (m_Mesh_faces.at(index_face).mark==false){
				m_Vertex_Border.push_back(m_CurrentEdge.v1);
				m_Vertex_Border.push_back(m_CurrentEdge.v2);
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
		while (m_CurrentVex.adjFace.empty() != true ){
			index_face = m_CurrentVex.adjFace.back();
			m_CurrentVex.adjFace.pop_back();
			if (m_Mesh_faces.at(index_face).mark==false){
				m_Mesh_faces.at(index_face).mark=true;
				m_Face_T0.push_back(index_face);//��T0�ҵ�
				m_Vertex_Free.push_back((m_Mesh_faces.at(index_face).v1));
				m_Vertex_Free.push_back((m_Mesh_faces.at(index_face).v2));
				m_Vertex_Free.push_back((m_Mesh_faces.at(index_face).v3));
			}//if
		}
		//ȥ�� �ϱ߽��, һ��������ֻͳ��һ��,�����ظ����(*iter) �� Free,ɾһ�ξ͹���
		m_Vertex_Free.remove((*iter));
	}//for
	m_Vertex_Free.sort();
	m_Vertex_Free.unique();
	return !m_Face_T0.empty();
}

void Mesh::ComputeCurrEnergy(){
	//����T0��������   m_Face_T0 ,m_Plane_T0
	E=0.0;
	E1=0.0;
	E2=0.0;
	PlanePara Pl; //ԭʼ�������������������
	Energy fu,fv;
	for (list<Index>iterator iter = m_Face_T0.begin();
													iter != m_Face_T0.end(); iter++){	
		m_CurrentTri = m_Mesh_faces.at(*iteir);
		m_CurrentPlane=ComputePlaneBasePos(m_CurrentTri);
		pl = PartialDerivative (m_CurrentPlane);
		fu = Pl.u;
		fv = Pl.v;
		E1 = fu*fu + fv*fv ;
		E2 = (fu*fu)*(fv*fv)-(fu*fv)*(fu*fv)-1;
		E2= E2*E2;
		E= m_epsilon * E1 + (1-m_epsilon) * e2;
	}

}

void Mesh::SloveMinEnery(){
	//����С������
}

// ���ɶ����ڳ�ʼƽ���ϵ�ͶӰ��
void Mesh::FreeVertexProjection(void){
	//m_Vertex_Free ; m_Vertex_Border; m_PlaneNormal; m_Face_T0;
	sVector footpoint; //���� 
	double t;//parameter_t
 m_Plane_T0.clear();
	for (list<Index>::iterator iter = m_Vertex_Free.begin();
											iter != m_Vertex_Free.end(); iter++)	{
				m_CurrentPos = (m_Mesh_vetexs.at(*iter )).pos;
				t = (m_CurrentPos - m_PlaneOrigin) * m_PlaneNormal;
				t = t / (m_PlaneNormal * m_PlaneNormal);
				footpoint = ( m_PlaneNormal * t ) + m_CurrentPos ;
				m_CurrentPlane.index = *iter ;
				m_CurrentPlane.u = footpoint * m_PlaneU ;
				m_CurrentPlane.v = footpoint * m_PlaneV ;
				m_Plane_T0.push_back(m_CurrentPlane);
	}		
}

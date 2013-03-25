// TriMesh 方法实现

#include "StdAfx.h"
#include "TriMesh.h"
//#include <fstream>
using namespace std;

TriMesh::TriMesh(){
	filename.clear();
	m_Mesh_faces.clear();
	m_Mesh_edges.clear();
	m_Mesh_vetexs.clear();
	 m_Edge_Border.clear();
}

TriMesh::~TriMesh(){
	
}

void TriMesh::ReadMesh(string filename){
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
	//read vetex
	Vetex dvetex;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;
	for(i=0;i<dV;i++){
		fscanf_s(in,"%lf %lf %lf",&dx,&dy,&dz);
		dvetex.SetVextxPos(dx,dy,dz);
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
	TriAngle dface;
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
	cout<<"Mesh Read.\n";
		cout<<"------------------------------------------"<<endl;
}

void TriMesh::ResultOutput(string filename){
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
//todo save代码 保存两个矩阵: 参数化后的顶点在平面域上的坐标(2维参量) 
		//                   对应三角形的拓扑 ( 顶点编号 )
		// 平面参数坐标集              vector<vetex>  m_Plane_Vertex
		//三角面 顶点编号              vector<Triangle> m_Mesh_faces
//write 文件头
	 fseek(in,0L,SEEK_SET);	
	 fprintf_s(in,"Result file : %s\n",file);
	 fprintf_s(in,"Source file : %s\n",sourcefile);
	 fprintf_s(in,"MeshFile format: ,%s\n",m_Mesh_format);
//write 顶点
	 fprintf_s(in,"Vetex Number: %d \n",m_Plane_Vertex.size());
	PlanePara vetex2d;
	for (int i=0; i!=m_Plane_Vertex.size(); i++){
		vetex2d.setParaPos(m_Plane_Vertex[i]);
		fprintf_s(in,"%lf\t%lf\n",vetex2d.u,vetex2d.v);
	}
	//write faces 
	fprintf_s(in,"Faces Number: %d \n",m_Mesh_faces.size());
	TriAngle face2d;
	for (int i=0; i!=m_Mesh_faces.size(); i++){
		face2d=m_Mesh_faces[i];
		fprintf_s(in,"%d\t%d\t%d\n",face2d.v1,face2d.v2,face2d.v3);
	}	
	fclose(in);
	cout<<"Mesh Saved.\n";
	cout<<"------------------------------------------"<<endl;
	cout<<endl;
}

void TriMesh::InitMesh(){
		//求各个边表  利用面栈 和边栈
	m_Mesh_T0.clear();	
	//m_Mesh_T0=m_Mesh_faces;
	m_Mesh_edges.clear();
	unsigned long vetexIndex=0;
	unsigned long edgeIndex=0;
	//unsigned long blakvetex=m_Mesh_vetexs.size()+1;

	for (vector<TriAngle>::iterator iter=m_Mesh_faces.begin();
		iter!=m_Mesh_faces.end(); iter++,vetexIndex++)	{
			//顶点邻接当前面
		m_CurrentEdge.adjFace.clear();
 		m_Mesh_vetexs.at(iter->v1).adjFace.push_back(vetexIndex);
		m_CurrentEdge.SetVextexIndex(iter->v1,iter->v2);
		m_CurrentEdge.adjFace.push_back(vetexIndex);
		m_Mesh_edges.push_back(m_CurrentEdge);
		iter->e1=edgeIndex++;

		m_Mesh_vetexs.at(iter->v2).adjFace.push_back(vetexIndex);
		m_CurrentEdge.SetVextexIndex(iter->v2,iter->v3);
		m_CurrentEdge.adjFace.at(0)=vetexIndex;
		m_Mesh_edges.push_back(m_CurrentEdge);
		iter->e2=edgeIndex++;

		m_Mesh_vetexs.at(iter->v3).adjFace.push_back(vetexIndex);
		m_CurrentEdge.SetVextexIndex(iter->v3,iter->v1);
		m_CurrentEdge.adjFace.at(0)=vetexIndex;
		m_Mesh_edges.push_back(m_CurrentEdge);
		iter->e3=edgeIndex++;//压边	
	}
	//交互不同方向的同一条边的邻接三角形编号
	for (vector<Edge>::iterator iter=m_Mesh_edges.begin(); iter !=m_Mesh_edges.end();iter++){

	}
}


void TriMesh::MeshesOutput(string filename){
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
	//todo 保存 初始化后的网格数据集;
	// 平面参数坐标集              vector<vetex>  m_Plane_Vertex
	//三角面 顶点编号              vector<Triangle> m_Mesh_faces
	//write 文件头
	fseek(in,0L,SEEK_SET);	
	fprintf_s(in,"Meshes file : %s\n",file);
	fprintf_s(in,"Source file : %s\n",sourcefile);
	fprintf_s(in,"MeshFile format: ,%s\n",m_Mesh_format);
	//write 顶点
	fprintf_s(in,"Vetex Number: %d \n",m_Mesh_vetexs.size());
	//PlanePara vetex2d;
	for (int i=0; i!=m_Mesh_vetexs.size(); i++){
		m_CurrentVex=m_Mesh_vetexs[i];
		fprintf_s(in,"%lf\t%lf\t%lf\n",m_CurrentVex.x,m_CurrentVex.y,m_CurrentVex.z,m_CurrentVex.adjFace);
			}
	//write faces 
	fprintf_s(in,"Faces Number: %d \n",m_Mesh_faces.size());
	//TriAngle face2d;
	for (int i=0; i!=m_Mesh_faces.size(); i++){
		m_CurrentTri=m_Mesh_faces[i];
		fprintf_s(in,"%ld\t%ld\t%ld\n",m_CurrentTri.v1,m_CurrentTri.v2,m_CurrentTri.v3);
	}
	//write edges
	fprintf_s(in,"Edges Number: %d \n",m_Mesh_edges.size());
	for (int i=0; i!=m_Mesh_edges.size(); i++)	{
		m_CurrentEdge=m_Mesh_edges[i];
		fprintf_s(in,"%ld\t%ld\n",m_CurrentEdge.v1,m_CurrentEdge.v2,m_CurrentEdge.adjFace);
	}
	fclose(in);
	cout<<"Mesh Saved.\n";
	cout<<"------------------------------------------"<<endl;
	cout<<endl;
}


void TriMesh::RunFlatPara(){
	m_Mesh_T2=m_Mesh_faces;
	m_Mesh_T1.clear();
	m_Mesh_T0.clear();
	while (m_Mesh_T2.empty()!=true)	{
		GetNextVetex();
//		ComputeCurrEnery();
//		SloveMinEnery();
	}
		
}

void TriMesh::GetNextVetex(){
	int nextTri;
	if (m_Mesh_T1.empty()==true)	{
		//放置第一个三角形
		nextTri=m_Mesh_faces.size()/2;
		m_CurrentTri=m_Mesh_faces[nextTri];
		m_Mesh_T0.push_back(m_CurrentTri);
		Flatten1stTir();
	}
	else{		
		GetNextTri();//查找出 t0
	}
	while (m_Mesh_T0.empty()!=true)	{
		m_CurrentTri=m_Mesh_T0.back();
		m_Mesh_T0.pop_back();
		m_Mesh_T1.push_back(m_CurrentTri);
		//swap(m_CurrentTri,m_Mesh_T2.back());
		//m_Mesh_T2.pop_back();
		for (vector<TriAngle>::iterator iter=m_Mesh_T2.begin();iter!=m_Mesh_T2.end(); iter++){
			if (m_CurrentTri.isEqual(*iter)) { m_Mesh_T2.erase(iter);break;}
		}
	}
	 m_Edge_Border.clear();
	GetT1Boundary();//更新T1边缘;存于m_Edge_Border

}

void TriMesh::Flatten1stTir(){
	//源三角形 放置于2D平面
}

void TriMesh::GetNextTri(){
	//找到与已处理的网格片的边界相邻的三角形..一次性找全找出来..
	//以顶点的邻接三角形?还是以边的邻接三角形?(以边的邻接三角形)存入 m_Mesh_T0 
	m_Mesh_T0.clear();
	for (vector<Edge>::iterator iter1=m_Edge_Border.begin(); iter1 !=m_Edge_Border.end();iter1++)	{
		for (vector<Edge>::iterator iter2=m_Mesh_edges.begin(); iter2 !=m_Mesh_edges.end();iter2++)		{
			if (iter1->v1==iter2->v2 && iter1->v2==iter2->v1){
				m_CurrentTri=m_Mesh_faces.at(iter2->adjFace.back());
				m_Mesh_T0.push_back(m_CurrentTri);
			}
		}
	}


}

void TriMesh::GetT1Boundary(){  
	for (int i=0; i<m_Mesh_T1.size(); i++)	{
		//int e=(int)m_Mesh_T1.at(i).v1;
		m_CurrentEdge=m_Mesh_edges.at(m_Mesh_T1.at(i).e1);//e1编号
		m_Edge_T1.push_back(m_CurrentEdge);
		m_CurrentEdge=m_Mesh_edges.at(m_Mesh_T1.at(i).e2);
		m_Edge_T1.push_back(m_CurrentEdge);
		m_CurrentEdge=m_Mesh_edges.at(m_Mesh_T1.at(i).e3);
		m_Edge_T1.push_back(m_CurrentEdge);
	}//压入边栈
	
	while (m_Edge_T1.empty()!=true){   //栈空否
		m_CurrentEdge=m_Edge_T1.back();
		m_Edge_T1.pop_back();			//弹出一个边
		if (m_Edge_T2.empty()!=true){			//池空否
			for (vector<Edge>::iterator iter=m_Edge_T2.begin();iter!=m_Edge_T2.end(); iter++)	{
				if (m_CurrentEdge.isEqual(*iter)) {
					m_Edge_T2.erase(iter); //删除池中相同边
					break;
				}
			}
			m_Edge_T2.push_back(m_CurrentEdge);// 边 压入边池
		}else{
			m_Edge_T2.push_back(m_CurrentEdge);// 边 压入边池
		}
	}
	//边界保存在 m_Edge_T2 中
	
	 m_Edge_Border=m_Edge_T2;
}

// bool TriMesh::CompareEdgesinT2(){	
// 	for (int i=0; i<m_Edge_T2.size(); i++)	{
// 		bool is=m_CurrentEdge.isEqual(m_Edge_T2.at(i));
// 		if (is==true) return true; 
// 	}
// 	return false;
// }
void TriMesh::ComputeCurrEnery(){
	//计算T1变形能量
}

void TriMesh::SloveMinEnery(){
	//解最小化问题
}
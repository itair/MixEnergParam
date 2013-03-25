// Mesh 方法实现

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
//todo save代码 保存两个矩阵: 参数化后的顶点在平面域上的坐标(2维参量) 
//  对应三角形的拓扑 ( 顶点编号 )
// 平面参数坐标集            vector<vetex>  m_Plane_Vertex
//三角面 顶点编号             vector<Triangle> m_Mesh_faces
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
	//求各个边表  利用面栈 和边栈
	//m_Mesh_T0.clear();	
	//m_Mesh_T0=m_Mesh_faces;
	m_Mesh_edges.clear();
	Index faceIndex=0;
	//Index vetexIndex=0;	
	for (vector<Face>::iterator iter=m_Mesh_faces.begin();
								iter!=m_Mesh_faces.end();)	{
			//顶点邻接当前面;
			iter->mark=false;
			m_Face_T2.push_back(faceIndex);
			//m_CurrentVex.adjEgde.clear();
			ProcessVetex(iter,faceIndex);
			ProcessEdge(iter,faceIndex);
			iter->index=faceIndex++ ;
			++iter;
	}
		
	SimplyEdge();//合并 重复统计的边/. 迭代器冒泡;
	ScanEdge(); //补全 面的临街面 点的邻接边;
	ProcessFace();
		
}

void Mesh::ProcessVetex(vector<Face>::iterator iter,Index faceIndex){
	//三个顶点 保存索引号 == 面片中的索引号
	m_Mesh_vetexs.at(iter->v1).index=iter->v1;
	m_Mesh_vetexs.at(iter->v2).index=iter->v2;
	m_Mesh_vetexs.at(iter->v3).index=iter->v3;

	//把当前面的索引号 加入 三个顶点的邻接面表;	
	m_Mesh_vetexs.at(iter->v1).adjFace.push_back(faceIndex);
	m_Mesh_vetexs.at(iter->v2).adjFace.push_back(faceIndex);
	m_Mesh_vetexs.at(iter->v3).adjFace.push_back(faceIndex);


}

void Mesh::ProcessEdge(vector<Face>::iterator iter,Index faceIndex){
		m_CurrentEdge.adjFace.clear();
	Index edgeIndex=0;
	//依次构造三条边 压入边表 加入索引;
	m_CurrentEdge.SetVexIndex(iter->v1,iter->v2); //内含序号升序;
	m_CurrentEdge.index=edgeIndex;
	iter->e1=edgeIndex++;
	m_CurrentEdge.adjFace.clear();
	m_CurrentEdge.adjFace.push_back(faceIndex);
	m_Edge_T1.push_back(m_CurrentEdge);

	m_CurrentEdge.SetVexIndex(iter->v2,iter->v3); //内含序号升序;
	m_CurrentEdge.index=edgeIndex;
	iter->e2=edgeIndex++;
	m_CurrentEdge.adjFace.clear();
	m_CurrentEdge.adjFace.push_back(faceIndex);
	m_Edge_T1.push_back(m_CurrentEdge);

	m_CurrentEdge.SetVexIndex(iter->v1,iter->v3); //内含序号升序;
	m_CurrentEdge.index=edgeIndex;
	iter->e3=edgeIndex++;
	m_CurrentEdge.adjFace.clear();
	m_CurrentEdge.adjFace.push_back(faceIndex);
	m_Edge_T1.push_back(m_CurrentEdge);
}

void Mesh::SimplyEdge(){
	//去除重复的边 重新编号;
	Index face1, face2,oldedgeIndex1,oldedgeIndex2;
	Index edgeIndex=0;
	vector<Edge>::iterator iter = m_Edge_T1.begin();		
	vector<Edge>::iterator iter2 = m_Edge_T1.begin();
	++iter2;	
	while(iter != m_Edge_T1.end()){
		if (iter->mark==true)  {
			++iter;
			iter2=iter;
			continue;//跳过重复的;
		}
		++iter2;
		while( iter2!=m_Edge_T1.end()){
			if ( *iter2 == *iter ){	//重复边	;	
				face2=iter2->adjFace.back();		
				oldedgeIndex2=iter2->index;//原边表数据 待改	;	
				m_Mesh_faces.at(face2).SwapEdgeIndex(oldedgeIndex2,edgeIndex);//更新对应邻接面 的边表;
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
	//扫描边表, 补充顶点的邻接边 和 边的邻接面;
	for (vector<Edge>::iterator iter=m_Mesh_edges.begin();iter != m_Mesh_edges.end();){
		m_Mesh_vetexs.at(iter->v1).adjEgde.push_back(iter->index);
		m_Mesh_vetexs.at(iter->v2).adjEgde.push_back(iter->index);

		m_Mesh_faces.at(iter->adjFace.front()).adjFace.push_back(iter->adjFace.back());
		m_Mesh_faces.at(iter->adjFace.back()).adjFace.push_back(iter->adjFace.front());

		++iter;
	}
		
}
sVector Mesh::NormalCross( sVector &v1, sVector&v2 ){
	//规范化叉积
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
	//计算 面表中 每个三角形的 法向量 填表;
	sVector a,b,c,ab,bc,norm;
	for (vector<Face>::iterator iter=m_Mesh_faces.begin();iter !=m_Mesh_faces.end();)	{
		a = m_Mesh_vetexs.at(iter->v1).pos;
		b = m_Mesh_vetexs.at(iter->v2).pos;
		c = m_Mesh_vetexs.at(iter->v3).pos;
 		ab = a - b;
 		bc = b - c;
		norm=NormalCross(ab,bc);
		iter->normal = norm;
		++iter;
	}

}

void Mesh::MeshesOutput(string filename){
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
		fprintf_s(in,"%lf\t%lf\t%lf\n",m_CurrentVex.pos.x,m_CurrentVex.pos.y,m_CurrentVex.pos.z,m_CurrentVex.adjFace);
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

bool Mesh::RunFlatPara(){
	m_Face_T1.clear();
	m_Face_T0.clear();
	m_Plane_Vertex.clear();
	bool check;
	//选取第一个三角形
	GetFirstTri();
	Flatten1stTir();
	//算法主循环
	while (m_Face_T2.empty()!=true)	{		
		check = GetNextVetex();
		if (check==false && m_Face_T2.empty()!=true) {
			cout<<"网格有洞 ,无法找到下一层T0" <<endl;
			return false;
		}	
		//list<Index> m_Face_T0;	// 当前 面栈
		//  list<Index> m_Vertex_Free;	//当前的 自由点

//TODO 计算 T0 能量表示;
//		ComputeCurrEnery();
//TODO 最小化问题 求出 m_Vertex_Free 对应 2D坐标 加入 m_Plane_Vertex
//		SloveMinEnery();

		//更新 T1 T2;
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
	//放置第一个三角形
	Index firstface=m_Mesh_faces.size()/2;
	m_Face_T0.push_back(firstface);
	m_Face_T1.push_back(firstface);
	m_Face_T2.remove(firstface);
	m_Mesh_faces.at(firstface).mark=true;
}

void Mesh::Flatten1stTir(){
	// 源三角形 放置于2D平面;
	Index index_face,index_vertex;
	sVector a,b,c,A,B,X,Y,N,XB;
	double more;
	//求一组正交基
	index_face=m_Face_T0.back();
	m_CurrentTri=m_Mesh_faces.at(index_face);
	a = m_Mesh_vetexs.at(m_CurrentTri.v1).pos;
	b = m_Mesh_vetexs.at(m_CurrentTri.v2).pos;
	c = m_Mesh_vetexs.at(m_CurrentTri.v3).pos;
	A =  b - a ;
	B =  c - a ;
	N = m_CurrentTri.normal;
	X = NormalCross( N, B ); //u方向
	Y = Cross( N, B );
	XB= Cross( X, B );
	more=More(XB);
	Y = Y / more;   //v方向
	// 嵌入 2d 平面;
	//v1 = (0,0);
	index_vertex = m_CurrentTri.v1;
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
	//T0中所有边中 邻面未mark的边 对应的顶点集 m_Vertex_Border;
	//m_Vertex_Borderi中全部邻面未mark的 为新T0;
	//新T0中 不在m_Vertex_Border中的顶点为自由顶点;
	Index index_face;
	//m_Edge_Border;
	for (list<Index>::iterator iter = m_Face_T0.begin();
				iter != m_Face_T0.end();	iter ++){
		m_CurrentTri = m_Mesh_faces.at(*iter); 
		m_Edge_Border.push_back(m_CurrentTri.e1);
		m_Edge_Border.push_back(m_CurrentTri.e2);
		m_Edge_Border.push_back(m_CurrentTri.e3);
	}
	m_Edge_Border.sort();
	m_Edge_Border.unique();
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
	for (list<Index>::iterator iter = m_Vertex_Border.begin();
				iter != m_Vertex_Border.end();	iter ++)	{
		m_CurrentVex = m_Mesh_vetexs.at(*iter);
		while (m_CurrentVex.adjFace.empty() != true ){
			index_face = m_CurrentVex.adjFace.back();
			m_CurrentVex.adjFace.pop_back();
			if (m_Mesh_faces.at(index_face).mark==false){
				m_Mesh_faces.at(index_face).mark=true;
				m_Face_T0.push_back(index_face);//新T0找到
				m_Vertex_Free.push_back((m_Mesh_faces.at(index_face).v1));
				m_Vertex_Free.push_back((m_Mesh_faces.at(index_face).v2));
				m_Vertex_Free.push_back((m_Mesh_faces.at(index_face).v3));
			}//if
		}
		//去除 老边界点, 一个三角形只统计一次,不会重复添加(*iter) 进 Free,删一次就够了
		m_Vertex_Free.remove((*iter));
	}//for
// 	if (m_Face_T0.empty()==true) return false; //下一层找不到了
// 	else return true;
	return !m_Face_T0.empty();
 }



void Mesh::GetNextTri(){
	//找到与已处理的网格片的边界相邻的三角形..一次性找全找出来..
	//以顶点的邻接三角形?还是以边的邻接三角形?(以边的邻接三角形)存入 m_Mesh_T0 
// 	m_Face_T0.clear();
// 	for (vector<Edge>::iterator iter1=m_Edge_Border.begin(); iter1 !=m_Edge_Border.end();iter1++)	{
// 		for (vector<Edge>::iterator iter2=m_Mesh_edges.begin(); iter2 !=m_Mesh_edges.end();iter2++)		{
// 			if (iter1->v1==iter2->v2 && iter1->v2==iter2->v1){
// 				m_CurrentTri=m_Mesh_faces.at(iter2->adjFace.back());
// 				m_Mesh_T0.push_back(m_CurrentTri);
// 			}
// 		}
// 	}
// 

}

void Mesh::GetT1Boundary(){  
// 	for (vector<Mesh>::size_type i=0; i< m_Mesh_T1.size(); i++)	{
// 		//int e=(int)m_Mesh_T1.at(i).v1;
// 		m_CurrentEdge=m_Mesh_edges.at(m_Mesh_T1.at(i).e1);//e1编号
// 		m_Edge_T1.push_back(m_CurrentEdge);
// 		m_CurrentEdge=m_Mesh_edges.at(m_Mesh_T1.at(i).e2);
// 		m_Edge_T1.push_back(m_CurrentEdge);
// 		m_CurrentEdge=m_Mesh_edges.at(m_Mesh_T1.at(i).e3);
// 		m_Edge_T1.push_back(m_CurrentEdge);
// 	}//压入边栈
// 
// 	while (m_Edge_T1.empty()!=true){   //栈空否
// 		m_CurrentEdge=m_Edge_T1.back();
// 		m_Edge_T1.pop_back();			//弹出一个边
//		if (m_Edge_T2.empty()!=true){			//池空否
// 			for (vector<Edge>::iterator iter=m_Edge_T2.begin();iter!=m_Edge_T2.end(); iter++)	{
// 				if (m_CurrentEdge.isEqual(*iter)) {
// 					m_Edge_T2.erase(iter); //删除池中相同边
// 					break;
// 				}
//			}
//			m_Edge_T2.push_back(m_CurrentEdge);// 边 压入边池
//		}else{
//			m_Edge_T2.push_back(m_CurrentEdge);// 边 压入边池
//		}
//	}
	//边界保存在 m_Edge_T2 中

//	m_Edge_Border=m_Edge_T2;
//}
}

// bool Mesh::CompareEdgesinT2(){	
// 	for (int i=0; i<m_Edge_T2.size(); i++)	{
// 		bool is=m_CurrentEdge.isEqual(m_Edge_T2.at(i));
// 		if (is==true) return true; 
// 	}
// 	return false;
// }
void Mesh::ComputeCurrEnery(){
	//计算T1变形能量
}

void Mesh::SloveMinEnery(){
	//解最小化问题
}
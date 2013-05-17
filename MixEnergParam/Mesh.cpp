// Mesh 方法实现

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
	strcpy_s(m_Mesh_format,ext);
	cout<<"MeshFile Format: " <<ext<<endl;
	//read vertex
	Vertex dvetex;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;
	for(i=0;i<dV;i++){
		fscanf_s(in,"%lf%lf%lf",&dx,&dy,&dz);
		dvetex.SetPos(dx,dy,dz);
		m_Mesh_vetexs.push_back(dvetex);
	}
	temp=m_Mesh_vetexs.size();
	assert( temp == dV );
	cout<<"Vetex Num: "<<dV<<endl;
	//read faces
	int val;
	int di;
	int dj;
	int dk;
	Face dface;
	for(i=0;i<dF;i++){
		fscanf_s(in,"%d%d%d%d",&val,&di,&dj,&dk);
		dface.SetVextexIndex(di,dj,dk);
		m_Mesh_faces.push_back(dface);
	}
	assert(m_Mesh_faces.size()==dF);
	cout<<"Face Num: "<<dF<<endl;
	//Read Edgs:
		fclose(in);
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
	double uv[2];
	for (int i=0; i!=m_Plane_Vertex.size(); i++){
				m_Plane_Vertex.at(i).Getuv(uv);
			fprintf_s(in, "%lf\t%lf\n",uv[0],uv[1]);
	}
	//write faces 
	fprintf_s(in,"Faces Number: %d \n",m_Mesh_faces.size());
	Index v123[3];
	for (int i=0; i!=m_Mesh_faces.size(); i++){
		m_Mesh_faces.at(i).Getv123(v123)	;
		fprintf_s(in,"%d\t%d\t%d\n", v123[0], v123[1], v123[2]);
	}	
	fclose(in);
	cout<<"Mesh Saved.\n";
	cout<<"------------------------------------------"<<endl;
	cout<<endl;
}

void Mesh::InitMesh(){
	//求各个边表  利用面栈 和边栈
	m_Mesh_edges.clear();
	Index faceIndex=0;
	//Index vetexIndex=0;	
	for (vector<Face>::iterator iter = m_Mesh_faces.begin();
															iter != m_Mesh_faces.end();)	{
			//顶点邻接当前面;
			(*iter).ResetMark();
			m_Face_T2.push_back(faceIndex);
			ProcessVetex(iter,faceIndex);
			ProcessEdge(iter,faceIndex);
			(*iter).SetIndex(faceIndex++) ;
			++iter;
	}		
	SimplyEdge();// 合并 重复统计的边/. 迭代器冒泡;
	ScanEdge(); // 补全 面的临街面 点的邻接边;
	ProcessFace();		
}

void Mesh::ProcessVetex(vector<Face>::iterator iter, Index faceIndex){
		(*iter).Getv123(v123);
	for (int i =0 ; i != 3; i++)	{
		// 三个顶点 保存索引号 == 面片中的索引号;
		m_Mesh_vetexs.at(v123[i]).SetIndex(faceIndex*3+i);
		// 把当前面的索引号加入 三个顶点的邻接面表;
		m_Mesh_vetexs.at(v123[i]).AddadjFace(faceIndex);
		m_Mesh_vetexs.at(v123[i]).ResetMark();
	}
}

void Mesh::ProcessEdge(vector<Face>::iterator iter,Index faceIndex){
	// 依次构造三条边 压入边表 加入索引;
	//Index *v123,*v12,*e123;
	Index e123[3];
	(*iter).Getv123(v123);
	Index edgeIndex=0;
	for (int i= 0; i != 3; i++){
		v12[0]=v123[i];
		v12[1]=v123[((i+1)%3)];
	
		m_CurrentEdge.SetVexIndex(v12[0],v12[1]); // 内含序号升序;
		edgeIndex = 3*faceIndex + i;
		m_CurrentEdge.SetIndex(edgeIndex);
		m_CurrentEdge.ClearAdjFace();
		m_CurrentEdge.ResetMark();
		m_CurrentEdge.AddadjFace(faceIndex);
		m_Edge_T1.push_back(m_CurrentEdge);
		e123[i] = edgeIndex;
	}
	(*iter).Sete123(e123);
}

void Mesh::SimplyEdge(){
	//去除重复的边 重新编号;
	Index face1, face2,oldEdge1,oldEdge2;
	Index edgeIndex=0;
	list<Edge>::iterator iter = m_Edge_T1.begin();
	list<Edge>::iterator iter2 = m_Edge_T1.begin();
	while(iter != m_Edge_T1.end()){
		if ( (*iter).IsMarked()	)		{
			iter++;
			continue;
		}
		iter2 = iter;
		iter2++;
		while (iter2 != m_Edge_T1.end()){
			if ((*iter2).IsEqual(*iter)){
				face2 = iter2->GetLastAdjFace();
				oldEdge2 = iter2->GetIndex();
				m_Mesh_faces.at(face2).SwapEdgeIndex(oldEdge2, edgeIndex);
				iter->AddadjFace(face2);
				(*iter2).GetMarked();
				break;
			}else{
			iter2++;
			}
		}
	face1 = iter->GetFirstAdjFace();
	oldEdge1 = iter->GetIndex(); 
	m_Mesh_faces.at(face1).SwapEdgeIndex(oldEdge1, edgeIndex);
	m_CurrentEdge=(*iter);
	m_CurrentEdge.SetIndex(edgeIndex);
	m_Mesh_edges.push_back(m_CurrentEdge);
	edgeIndex++;
	iter++;
	}	
}
void Mesh::ScanEdge(){
	//扫描边表, 补充顶点的邻接边 和 边的邻接面;
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
 
 sVector Mesh::NormalCross( const sVector &v1, const sVector &v2){   //规范化叉积		
	sVector norm;
	double more; 	
	norm= norm.Cross( v1 , v2 );//无返回值!!!
	more = norm.More();	
	return  norm/ more;	
}
void Mesh::ProcessFace(){
	//计算 面表中 每个三角形的 法向量 和 局部正交基 平面坐标 填表; 还有面积.
	sVector a,b,c,A,B,X,Y,N;
	double Area,p2u,p2v,p3u,p3v;
	//Index *v123;	
	PlanePara p123[3];
	for (vector<Face>::iterator iter=m_Mesh_faces.begin(); 
															iter !=m_Mesh_faces.end();)	{
		iter->Getv123(v123); 		
		a = m_Mesh_vetexs.at(v123[0]).GetPos();
		b = m_Mesh_vetexs.at(v123[1]).GetPos();
		c = m_Mesh_vetexs.at(v123[2]).GetPos();
 		//求一组正交基		
		A =  b - a ;
		B =  c - a ;
		N = NormalCross( A, B); // 法向量问题
		iter->SetNormal( N );
		X = NormalCross( N, B); //u方向		
		Y = NormalCross( X, B); //v方向
		p2u = A*X;
		p2v = A*Y;
		p3u = B*X;
		p3v = B*Y;

		m_CurrentPlane.Setuv(0.0 , 0.0);
		p123[0]=m_CurrentPlane;
		m_CurrentPlane.Setuv( p2u , p2v);
		p123[1]=m_CurrentPlane;
		m_CurrentPlane.Setuv( p3u , p3v);
		p123[2]=m_CurrentPlane;
		iter->Setp123(p123);

		 Area = 0.5 * ( (p2u * p3v) - (p2v * p3u));
		 (*iter).SetArea (fabs(Area));
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
	//TODO 保存 初始化后的网格数据集;
	// 平面参数坐标集              vector<vetex>  m_Plane_Vertex
	//三角面 顶点编号              vector<Triangle> m_Mesh_faces
	//write 文件头
	fseek(in,0L,SEEK_SET);	
	fprintf_s(in,"Meshes file : %s\n",file);
	fprintf_s(in,"Source file : %s\n",sourcefile);
	fprintf_s(in,"MeshFile format: %s\n",m_Mesh_format);
	//write 顶点
	
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
		m_CurrentTri.Getv123(v123);
		Index e123[3];
		m_CurrentTri.Gete123(e123);
		PlanePara p123[3];
		m_CurrentTri.Getp123(p123);
		sVector norm = m_CurrentTri.GetNormal();
		double xyz[3] ;
		norm.GetXYZ(xyz);
		double Area = m_CurrentTri.GetArea();
		fprintf_s(in, "%lu\t%lu\t%lu\t", v123[0],v123[1],v123[2]);
		fprintf_s(in, "%lf\t%lf\t%lf\n", xyz[0],xyz[1],xyz[2]);
		fprintf_s(in, "%lu\t%lu\t%lu\n", e123[0],e123[1],e123[2]);
		fprintf_s(in , "%lf\n", Area);
		for (int i=0; i<3 ; i++){
			double uv[2];
			p123[i].Getuv(uv);
			fprintf_s(in, "%lf\t%lf\n", uv[0],uv[1]);
		}
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
	//选取第一个三角形
	GetFirstTri();
	Flatten1stTir();
	//算法主循环
	while (m_Face_T2.empty()!=true)	{		
		check = GetNextVetex();
		if (check==false && m_Face_T2.empty()!=true) {
			cout<<"网格有洞, 无法找到下一层T0... " <<endl;
			return false;
		  }	
		FreeVertexProjection();   //计算自由点 在 初始平面的上的投影 并转换为平面坐标 uv 存入 m_Plane_T0 
		//ComputeCurrEnergy()  ;   // 计算 T0 能量表示;
//TODO 最小化问题 求出 m_Vertex_Free 对应 2D坐标 加入 m_Plane_Vertex
//		SloveMinEnery();	

		m_Plane_Vertex.insert(m_Plane_Vertex.end(),m_Plane_T0.begin(),m_Plane_T0.end());

		for (list<Index>::iterator iter = m_Face_T0.begin();
				iter != m_Face_T0.end();	iter ++){  //更新 T1 T2;
			m_Face_T2.remove(*iter);
			m_Face_T1.push_back(*iter);
		}
}
	return true;
}
void Mesh::GetFirstTri(){
	//放置第一个三角形
	Index firstpoint,firstface;
	double xyz[3];
	double temp,temp2;
	m_Mesh_vetexs.at(0).GetXYZ(xyz);
	temp=fabs(xyz[0])+fabs(xyz[1])+fabs(xyz[2]);
	firstpoint=0;
	for (vector<Vertex>::iterator iter = m_Mesh_vetexs.begin();
		iter != m_Mesh_vetexs.end();	iter ++){
		(*iter).GetXYZ(xyz);
		temp2=fabs(xyz[0])+fabs(xyz[1])+fabs(xyz[2]);
		if (temp>temp2) {	
			temp=temp2;
			firstpoint=iter->GetIndex();
			firstface=iter->GetFirstAdjFace();
		}
		
 	}
 	cout<<"firstpoint= "<<firstpoint<<endl;
// 	Index v123[3];
// 	firstface=0;
// 	for (vector<Face>::iterator iter = m_Mesh_faces.begin();
// 		iter != m_Mesh_faces.end();	iter ++){
// 			iter->Getv123(v123);
// 			if (v123[0]==firstpoint || v123[1]==firstpoint || v123[2]==firstpoint){
// 				break;
// 			}
// 			firstface++;
// 	}
	if (firstface >=m_Mesh_faces.size() ){
		cout<<"无法找到firstface = "<<firstface<<endl;
	}else{
		cout<<"firstface= "<<firstface<<endl;
	}


	//firstface=m_Mesh_faces.size()/2;
	m_Face_T0.push_back(firstface);
	m_Face_T1.push_back(firstface);
	m_Face_T2.remove(firstface);
	m_Mesh_faces.at(firstface).GetMarked();
}

void Mesh::Flatten1stTir(){
	// 源三角形 放置于2D平面;
	Index index_face,index_vertex;	
	Index v123[3];	
	PlanePara p123[3];
	index_face=m_Face_T0.back();
	m_CurrentTri = m_Mesh_faces.at(index_face);
	m_CurrentTri.Getv123(v123); 		
	sVector a = m_Mesh_vetexs.at(v123[0]).GetPos();
	sVector b = m_Mesh_vetexs.at(v123[1]).GetPos();
	sVector c = m_Mesh_vetexs.at(v123[2]).GetPos();
	//求一组正交基		
	sVector A =  b - a ;
	sVector B =  c - a ;
	sVector N = NormalCross( A, B); // 法向量问题
	m_PlaneU = NormalCross( N, B); //u方向		
	m_PlaneV = NormalCross( m_PlaneU, B); //v方向

	m_CurrentTri.Getv123(v123);
	m_CurrentTri.Getp123(p123);	
 	m_PlaneNormal = m_CurrentTri.GetNormal(); //嵌入平面的统一法向量
	for (int i = 0; i< 3 ; i++)	{
		index_vertex = v123[i];
		m_CurrentPos = m_Mesh_vetexs.at(index_vertex).GetPos();
		if (i == 0) { m_PlaneOrigin = m_CurrentPos; }    //2D原点
		m_CurrentPlane.Clear();
		m_CurrentPlane = p123[i];
		m_CurrentPlane.SetIndex( index_vertex );
		m_Plane_Vertex.push_back(m_CurrentPlane);	
	}	
}

bool Mesh::GetNextVetex(){
	
	Index index_face,index_edge;
	Index e123[3];
	Index v123[3];
	vector<Index> adjface_index;
	bool marked;
	//t0 所有未被mark的点 为边界点,mark它们; 边界点所有邻接面中未被mark的为新T0
	//新T0 中所有未被mark的点为 自由点,mark新t0
	for (list<Index>::iterator iter = m_Face_T0.begin();
		iter != m_Face_T0.end();	iter ++){
			index_face = *iter;
			m_Mesh_faces.at(index_face).Getv123(v123);
			for (int i=0; i!=3; i++){
				m_CurrentVex = m_Mesh_vetexs.at(v123[i]);
				if ( !m_CurrentVex.IsMarked()){
					m_Mesh_vetexs.at(v123[i]).GetMarked();
					m_Vertex_T1.push_back(v123[i]);					
					}
				}
			}	
	if (m_Vertex_T1.empty()){return false;}
	m_Vertex_T1.sort();
	m_Vertex_T1.unique();
	m_Face_T0.clear();
	for (list<Index>::iterator iter = m_Vertex_T1.begin();
		iter != m_Vertex_T1.end();	iter ++){
		adjface_index = m_Mesh_vetexs.at(*iter).GetadjFace();
		for (vector<Index>::iterator iter2 = adjface_index.begin();
			iter2 != adjface_index.end(); iter2 ++){
				index_face = *iter2;
				m_CurrentTri = m_Mesh_faces.at(index_face);
				if ( !m_CurrentTri.IsMarked())		{//新T0
					m_Face_T0.push_back(index_face);
					m_Mesh_faces.at(index_face).GetMarked();
				}
			}
		}
	m_Face_T0.sort();
	m_Face_T0.unique();
	m_Vertex_Free.clear();
	for (list<Index>::iterator iter = m_Face_T0.begin();
		iter != m_Face_T0.end();	iter ++){
			index_face = *iter;
			m_Mesh_faces.at(index_face).Getv123(v123);
			for (int i=0; i!=3; i++){
				m_CurrentVex = m_Mesh_vetexs.at(v123[i]);
				if ( !m_CurrentVex.IsMarked()){
					//m_Mesh_vetexs.at(v123[i]).GetMarked();
					m_Vertex_Free.push_back(v123[i]);					
				}
			}
	}	
	m_Vertex_Free.sort();
	m_Vertex_Free.unique();
	return !m_Face_T0.empty();
}

void Mesh::ComputeCurrEnergy(){
	//计算T0变形能量   m_Face_T0 ,m_Plane_T0	
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
		PDerivative pl; //原始正交基下三角面的坐标
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
	//取得 E  初始值 对应当前T0
}
// 自由定点在初始平面上的投影点
void Mesh::FreeVertexProjection(void){
	//m_Vertex_Free ; m_Vertex_Border; m_PlaneNormal; m_Face_T0;
	sVector footpoint; //垂足 
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
	//求偏导
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



// 检查初始化网格数据
bool Mesh::CheckMesh(void)
{
	Index vetNum = m_Mesh_vetexs.size();
	Index faceNum = m_Mesh_faces.size();
	Index edgeNum = m_Mesh_edges.size();

	Index index_check;
	Index e123[3];
	Index v123[3];
	//PlanePara p123[3];
	double area_check;
	const char* errLoc; 
	vector<Index> adjFace_check;
	
		for (vector<Face>::iterator iter = m_Mesh_faces.begin();
				iter != m_Mesh_faces.end(); iter++ ){
		try{ //面检测
					index_check=(*iter).GetIndex();
			if (index_check >faceNum) { throw "面片序号index 溢出"; }
			(*iter).Getv123(v123);
			if (v123[0]>vetNum ||v123[1]>vetNum ||v123[2]>vetNum ){
				throw "顶点编号v123 溢出";
			}
			(*iter).Gete123(e123);
			if (e123[0]>edgeNum ||e123[1]>edgeNum ||e123[2]>edgeNum ){
				throw "顶点编号e123 溢出";
			}
			area_check=(*iter).GetArea();
			if (area_check <= 0) {
				throw "面积数据错误";
				}
			adjFace_check=(*iter).GetAdjFaces();
			for (vector<Index>::iterator iter2=adjFace_check.begin();
				iter2 !=adjFace_check.end(); iter2++ ){
					if ((*iter2)>faceNum)	{
						throw "面表邻接面序号溢出";
					}
				}
			}//try
	catch(const char* errLoc ){
			cout<<"m_Mesh_faces 数据错误"<<endl;
			cout<<errLoc<<endl;
			abort();	
			}
		}	 
	Index v12[2];	
		for (vector<Edge>::iterator iter = m_Mesh_edges.begin();
			iter != m_Mesh_edges.end(); iter ++){
				try{//边表检查
					index_check=(*iter).GetIndex();
					if (index_check >edgeNum) { throw "边表序号index 溢出"; }
					(*iter).Getv12(v12);
					if (v12[0]>vetNum || v12[1] >vetNum){throw "边表顶点序号溢出";}
					adjFace_check = (*iter).GetAdjFaces();
					for (vector<Index>::iterator iter2=adjFace_check.begin();
						iter2 !=adjFace_check.end(); iter2++ ){
							if ((*iter2)>faceNum)	{
								throw "边表邻接面序号溢出";
							}
					}
				}
				catch (const char* errLoc){
					cout<<"m_Mesh_edges 数据错误"<<endl;
						cout<<errLoc<<endl;
						abort();	
				}
		}
	return true;
}
void Mesh::SloveMinEnery(){
	//梯度下降法 解最小化问题;
	//输入: 初始能量 E , 目标函数 E(Tx), 迭代步长 a , 收敛精度 e ;
	//当然还得知道 梯度表达式 ;如果 梯度过小导致迭代步数太对 可以使用截断梯度计算或者随机梯度计算
	//输入 m_Face_T0 ; m_Vertex_Free;

}
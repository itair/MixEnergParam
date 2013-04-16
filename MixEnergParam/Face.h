/*************************************************************
	三角网格的基本元素 :
	三维三角面片 : Face
	声明: Face.h
	@by Itair
***********************************************************/
#pragma once
#include "BasicDS.h"
#include "sVector.h"
#include "PlanePara.h"
#include <vector>

class PlanePara;

class Face{
private:
	Index index;
	Index v1,v2,v3;
	Index e1,e2,e3;
	bool mark; //多用途标记....
	PlanePara p1,p2,p3;
	sVector normal;
	double Area;
	vector <Index> adjFace;
public:	
	void SetVextexIndex(const Index a, const Index b, const Index c){
		v1 = a;		v2= b;		v3 = c;
	}
	void SetEdgeIndex(const Index a, const Index b, const Index c){
		e1 = a;		e2 = b;		e3 = c;
	}
	void SwapEdgeIndex(const Index old,  const Index newEdge ) {
		if (old == e1) { e1 = newEdge ; return; }
		if (old == e2) { e2 = newEdge ; return; }
		if (old == e3) { e3 = newEdge ; return; }		
	}
	bool operator == (const Face& fa){		return ( index == fa.index );	}	

	const Index GetIndex(){		return index;	}
	void SetIndex(const Index i){		this->index = i;	}

	const sVector GetNormal(void){  return normal; }
	void SetNormal(const sVector norm){ normal = norm; }

	void  Getv123(Index v123[]){
		v123[0] = v1;
		v123[1] = v2;
		v123[2] = v3;
	}
	void  Setv123(Index v123[]){
		v1 = v123[0] ;
		v2 = v123[1] ;
		v3 = v123[2] ;
	}
	void  Gete123(Index e123[]){
		e123[0] = e1;
		e123[1] = e2;
		e123[2] = e3;
	}
	void  Sete123(Index e123[]){
		e1 = e123[0] ;
		e2 = e123[1] ;
		e3 = e123[2] ;
	}
	 void Getp123(PlanePara p123[]){
		p123[0] = p1;
		p123[1] = p2;
		p123[2] = p3;
	}
	void Setp123(PlanePara p123[]){
		p1 = p123[0] ;
		p2 = p123[1] ;
		p3 = p123[2] ;
	}
	void SetArea(const double TriAera) { Area = TriAera; }
	const double GetArea(void){ return Area;}

	void ResetMark(void){ 		this->mark = false;	}
	bool IsMarked(void){		return mark;	}
	void GetMarked(void){    this->mark = true; }
	void ClearAdjFace(void){	adjFace.clear();	}

	const Index GetFirstAdjFace(void){	return adjFace.front();	}
	const Index GetLastAdjFace(void){	return adjFace.back();	}
	const Index GetAdjFaceSize(void){	return adjFace.size();	}		
	void AddadjFace(const Index faceIndex){	this->adjFace.push_back(faceIndex);	}

};

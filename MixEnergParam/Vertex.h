/*************************************************************
	三角网格的基本元素 :
	三维顶点: Vertex 
	声明: Vertex.h
	@by Itair
***********************************************************/
#pragma once
#include "BasicDS.h"
#include "sVector.h"
#include <vector>
class sVector;
class Vertex {
private:
	Index index;
	sVector pos;
	bool mark; //多用途标记....
	vector <Index> adjFace;
	vector <Index> adjEgde;
public:
	void SetPos(const double  dx, const double dy,const double dz){
		pos.SetsVector(dx, dy, dz);
	}
	const sVector GetPos(){		return pos; }	
	void GetXYZ (double xyz[]) { 	pos.GetXYZ(xyz);	}

	const Index GetIndex(){		return index;	}
	void SetIndex(const Index i){		this->index=i;	}	

	void AddadjEdge(const Index EdgeIndex){	this->adjEgde.push_back(EdgeIndex);	}

	void ResetMark(void){		this->mark = false;	}
	void GetMarked(void){ this->mark = true; }
	bool IsMarked(void){		return mark;	}

	bool IsAdjFaceEmpty(void){		return adjFace.empty();	}
	void ClearAdjFace(void){	adjFace.clear();	}
	const Index GetFirstAdjFace(void){	return adjFace.front();	}
	const Index GetLastAdjFace(void){	return adjFace.back();	}
	const Index GetAdjFaceSize(void){	return adjFace.size();	}
	void AddadjFace(const Index faceIndex){	this->adjFace.push_back(faceIndex);	}
	void DeladjFace(void){	this->adjFace.pop_back();	}
	const vector<Index> GetadjFace(void){		return adjFace;	}

};
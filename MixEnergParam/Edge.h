/*************************************************************
	三角网格的基本元素 :
	三维边 : Edge 
	声明: Edge.h
	@by Itair
***********************************************************/
#pragma once
#include "BasicDS.h"
#include <vector>

class Edge{
private:
	Index index;
	Index v1,v2;
	vector <Index> adjFace;
	bool mark; //多用途标记....
public:
	
	Edge(){
		mark=false;	
	}
	
	void SetVexIndex(const Index a, const Index b){
		this->v1 = ( a < b ) ? a : b; 
		this->v2 = ( a >= b )? a : b;		
	}
	
// 	bool operator == (const Edge& edge){
// 		return (v1==edge.v1 && v2==edge.v2);
// 	}
	bool IsEqual(const const Edge& edge ){
		return (v1==edge.v1 && v2==edge.v2 );
	}

	void ClearAdjFace(void){	adjFace.clear();	}
	const Index GetFirstAdjFace(void){	return adjFace.front();	}
	const Index GetLastAdjFace(void){	return adjFace.back();	}
	const Index GetAdjFaceSize(void){	return adjFace.size();	}
	void AddadjFace(const Index faceIndex){	this->adjFace.push_back(faceIndex);	}
	void DeladjFace(void){	this->adjFace.pop_back();	}
	bool IsAdjFaceEmpty(void){		return adjFace.empty();	}

	void  Getv12(Index v12[]){		v12[0] = v1;		v12[1] = v2;			}
	 void Setv12(Index v12[]){		v1 = v12[0] ;		v2 = v12[1] ;		}

	const Index GetIndex(){		return index;	}
	void SetIndex(const Index i){		this->index=i;	}

	void ResetMark(void){		this->mark = false;	}
	void GetMarked(void){ this->mark = true; }
	bool IsMarked(void){		return mark;	}
};

class Border{
	Index v1;   // 第1 点索引 
	Index v2;   //第2 点索引 
	Index faceIndex;    // 所属面片索引  
};
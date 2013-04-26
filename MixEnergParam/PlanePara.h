/*************************************************************
	三角网格的基本元素 :
	平面网格: PlanePara 
	声明: PlanePara.h
	@by Itair
***********************************************************/
#pragma once
#include "BasicDS.h"

class PlanePara{
private:
	Index index; //对应顶点的 索引
	double u;
	double v;
public:
	void Getuv(double uv[]){		uv[0]=this->u;	uv[1]=this->v;	}

	void Setuv(const double uu, const double vv){ u = uu; v = vv; }

	const Index  GetIndex(){ return index;	}

	void SetIndex ( const Index i) { index = i ;}
	
	void Clear(void) { index = 0, u=0.0 ;v=0.0 ;}
};
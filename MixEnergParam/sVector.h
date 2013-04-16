/*************************************************************
	三角网格的基本元素 :
	三维向量 : sVector 
	声明: sVector.h
	@by Itair
***********************************************************/
#pragma once
#include "BasicDS.h"

class sVector{
private:
	double x,y,z;
	//重载
public:
	const sVector operator+ (const sVector &sv) {
		sVector res;
		res.x = x + sv.x;
		res.y = y + sv.y;
		res.z = z + sv.z;
		return res;
	}
	const sVector operator- (const sVector &sv) {
		sVector res;
		res.x = x - sv.x;
		res.y = y - sv.y;
		res.z = z - sv.z;
		return res;
	}
	const sVector operator/ (const double real){
		sVector res;
		res.x = x / real;
		res.y = y / real;
		res.z = z / real;
		return res;
	}
	double operator* (const sVector &sv){
		double xx,yy,zz;
		xx = x * sv.x;
		yy = y * sv.y;
		zz = z * sv.z;
		return  (xx*xx + yy*yy+ zz*zz);
	}
	sVector operator* (const double real){
		sVector res;
		res.x = x * real;
		res.y = y * real;
		res.z = z * real;
		return res;
	}
	//方法
public:
	void SetsVector(const double xx, const double yy, const double zz) {
		x=xx; y=yy; z=zz;
	}

	void	GetXYZ (double xyz[]){ xyz[0]=x; xyz[1]=y; xyz[2]=z; }

	const sVector Cross(const sVector &ab, const sVector &bc){
		sVector norm;		
		norm.x = (ab.y * bc.z) - (ab.z * bc.y);
		norm.y = -((ab.x * bc.z) - (ab.z * bc.x));
		norm.z = (ab.x * bc.y) - (ab.y * bc.x);
		return norm;
	}	
	const double More( void){		return sqrt(x*x + y*y + z*z);	}

};
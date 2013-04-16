/*************************************************************
	�Զ��帴����ƫ������ �� ��������
	����ƫ������: PDerivative
	����: PDerivative.h
	@by Itair
	***********************************************************/
#pragma once
#include "BasicDS.h"

class PDerivative{
private:
	Complex fu;
	Complex fv;
public:
	const double  DotProduct(const Complex fu, const Complex fv ){
		double Real = real(fu) * real(fv) - imag(fu)*imag(fv);
		return Real;
	}
	void Setfu(Complex ffu) { fu = ffu; }
	void Setfv(Complex ffv) { fv = ffv; }
	const Complex Getfu(void) {return fu;}
	const Complex Getfv(void) {return fv; }


}; 


#pragma once

class GenerateMatrix {
private:
	int size;
	double** matrixA;
	double* vectorB;
	double** matrixL;
	double** matrixU;
	int a1, a2, a3;
	int f;
public:
	GenerateMatrix(int N, int a1, int a2, int a3, int f);
	~GenerateMatrix();
	void generate();
	int getSize();
	double** getMatrixA();
	double* getVectorB();
	void printMatrixA();
	void printVectorB();
	double getNormUsingResidual(double* x);
	double** getMatrixL();
	double** getMatrixU();
	void performLUFactorization();
};

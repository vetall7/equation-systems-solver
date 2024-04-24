#pragma once
#include "GenerateMatrix.h"
#include <vector>
using namespace std;

class SolveEquation {
	GenerateMatrix* matrix_generator;
	int size;
	double* x;
	vector<double> norms;
	double time;
public:
	SolveEquation(GenerateMatrix* matrix_generator);
	~SolveEquation();
	double* solveJacobi();
	double* solveGaussSeidel();
	double* solveFactorizationLU();
	vector<double> getNorms();
	double getTime();
private:
	double* forwardSubstituition(double** L, double* b);
	void backwardSubstituition(double** U, double* b);
};
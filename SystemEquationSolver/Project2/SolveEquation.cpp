#include "SolveEquation.h"
#include <math.h>
#include <iostream>
#include <chrono>
using namespace std;

SolveEquation::SolveEquation(GenerateMatrix* matrix_generator) {
	this->matrix_generator = matrix_generator;
	this->size = matrix_generator->getSize();
	this->x = new double[size] {};
}


double* SolveEquation::solveJacobi() {
	double** A = this->matrix_generator->getMatrixA();
	double* b = this->matrix_generator->getVectorB();
	double* x_new = new double[size];
	double eps = pow(10, -9);
	double norm = 0;
	double sum = 0;
	int iter = 0;
	int max_iter = 100;
	auto start = chrono::high_resolution_clock::now();
	do {
		for (int i = 0; i < size; i++) {
			sum = 0;
			for (int j = 0; j < size; j++) {
				if (i != j) {
					sum += A[i][j] * x[j];
				}
			}
			x_new[i] = (b[i] - sum) / A[i][i];
		}
		norm = this->matrix_generator->getNormUsingResidual(x_new);
		this->norms.push_back(norm);
		for (int i = 0; i < size; i++) {
			x[i] = x_new[i];
		}
		iter++;
	} while (norm > eps && iter < max_iter);
	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
	this->time = duration.count();
	cout << "Time Jaccobi: " << time << " ms" << endl;
	cout << "Number of iterations Jaccobi : " << iter << endl;
	cout << "Norm Jaccobi : " << norm << endl;
	delete [] x_new;
	return x;
}

double* SolveEquation::solveGaussSeidel() {
	double** A = this->matrix_generator->getMatrixA();
	double* b = this->matrix_generator->getVectorB();
	this->norms.clear();
	double* x_new = new double[size] {};
	double eps = pow(10, -9);
	double norm = 0;
	double sum = 0;
	int iter = 0;
	auto start = chrono::high_resolution_clock::now();
	int max_iter = 100;
	do {
		for (int i = 0; i < size; i++) {
			sum = 0;
			for (int j = 0; j < size; j++) {
				if (i != j) {
					sum += A[i][j] * x_new[j];
				}
			}
			x_new[i] = (b[i] - sum) / A[i][i];
		}
		norm = this->matrix_generator->getNormUsingResidual(x_new);
		this->norms.push_back(norm);
		iter++;
	} while (norm > eps && iter < max_iter);
	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
	this->time = duration.count();
	cout << "Time Gauss-Seidel : " << time << " ms" << endl;
	cout << "Number of iterations Gauss-Seidel : " << iter << endl;
	cout << "Norm Gauss-Seidel : " << norm << endl;
	for (int i = 0; i < size; i++) {
		x[i] = x_new[i];
	}
	delete[] x_new;
	return x;
}


double* SolveEquation::solveFactorizationLU() {
	this->norms.clear();
	double** A = this->matrix_generator->getMatrixA();
	double* b = this->matrix_generator->getVectorB();
	auto start = chrono::high_resolution_clock::now();
	this->matrix_generator->performLUFactorization();
	double** L = this->matrix_generator->getMatrixL();
	double** U = this->matrix_generator->getMatrixU();

	double* y = forwardSubstituition(L, b);
	backwardSubstituition(U, y);
	delete[] y;
	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

	double norm = this->matrix_generator->getNormUsingResidual(x);
	this->time = duration.count();
	cout << "Time LU : " << time << " ms" << endl;
	cout << "Norm LU : " << norm << endl;
	return this->x;
}

double* SolveEquation::forwardSubstituition(double** L, double* b){
	double* y = new double[size] {};
	for (int i = 0; i < size; i++) {
		y[i] = b[i];
		for (int j = 0; j < i; j++) {
			y[i] -= L[i][j] * y[j];
		}
	}
	return y;
}

void SolveEquation::backwardSubstituition(double** U, double* b) {
	for (int i = size - 1; i >= 0; i--) {
		this->x[i] = b[i];
		for (int j = size - 1; j > i; j--) {
			this->x[i] -= U[i][j] * this->x[j];
		}
		this->x[i] /= U[i][i];
	}
}

double SolveEquation::getTime() {
	return time;
}

vector<double> SolveEquation::getNorms() {
	return norms;
}

SolveEquation::~SolveEquation() {
	delete[] x;
}
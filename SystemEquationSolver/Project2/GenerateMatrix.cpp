#include <iostream>
#include <cmath>
#include "GenerateMatrix.h"

GenerateMatrix::GenerateMatrix(int N, int a1, int a2, int a3, int f) {
	this->size = N;
	this->a1 = a1;
	this->a2 = a2;
	this->a3 = a3;
	this->f = f;
	this->matrixA = new double* [N];
	for (int i = 0; i < N; i++) {
		matrixA[i] = new double[N] {};
	}
	this->vectorB = new double[N];
}

void GenerateMatrix::generate() {
	for (int i = 0; i < this->size; i++) {
		this->matrixA[i][i] = a1;
		if (i + 1 < this->size) {
			this->matrixA[i][i + 1] = a2;
			this->matrixA[i + 1][i] = a2;
		}
		if (i + 2 < this->size) {
			this->matrixA[i][i + 2] = a3;
			this->matrixA[i + 2][i] = a3;
		}
	}
	for (int i = 0; i < this->size; i++) {
		this->vectorB[i] = sin(i*(this->f+1));
	}
}

int GenerateMatrix::getSize() {
	return this->size;
}

double** GenerateMatrix::getMatrixA() {
	return this->matrixA;
}

double* GenerateMatrix::getVectorB() {
	return this->vectorB;
}

void GenerateMatrix::printMatrixA() {
	for (int i = 0; i < this->size; i++) {
		for (int j = 0; j < this->size; j++) {
			std::cout << this->matrixA[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void GenerateMatrix::printVectorB() {
	for (int i = 0; i < this->size; i++) {
		std::cout << this->vectorB[i] << " ";
	}
	std::cout << std::endl;
}

double GenerateMatrix::getNormUsingResidual(double* x) {
	double* residual = new double[size] {};
	double norm = 0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			residual[i] += matrixA[i][j] * x[j];
		}
		residual[i] -= vectorB[i];
		norm += residual[i] * residual[i];
	}
	norm = sqrt(norm);
	delete[] residual;
	return norm;
}

double** GenerateMatrix::getMatrixL() {
	return this->matrixL;
}

double** GenerateMatrix::getMatrixU() {
	return this->matrixU;
}

void GenerateMatrix::performLUFactorization() {
	// Allocate memory for matrices L and U
	matrixL = new double* [size];
	matrixU = new double* [size];
	for (int i = 0; i < size; ++i) {
		matrixL[i] = new double[size] {};
		matrixU[i] = new double[size] {};
	}

	// Copy matrix A to U
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			matrixU[i][j] = matrixA[i][j];
		}
	}

	// Perform LU factorization with partial pivoting
	for (int k = 0; k < size - 1; ++k) {
		for (int i = k + 1; i < size; ++i) {
			double factor = matrixU[i][k] / matrixU[k][k];
			matrixL[i][k] = factor;
			for (int j = k; j < size; ++j) {
				matrixU[i][j] -= factor * matrixU[k][j];
			}
		}
	}

	// Set diagonal of L to 1
	for (int i = 0; i < size; ++i) {
		matrixL[i][i] = 1.0;
	}
}

GenerateMatrix::~GenerateMatrix() {
	for (int i = 0; i < size; i++) {
		delete[] matrixA[i];
	}
	delete[] matrixL;
	delete[] matrixU;
	delete[] matrixA;
	delete[] vectorB;
}
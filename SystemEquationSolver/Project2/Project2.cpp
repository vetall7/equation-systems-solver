#include <iostream>
#include "GenerateMatrix.h"
#include "SolveEquation.h"
#include <fstream>
using namespace std;

void plot_norms(vector<double>& normsJacobi, vector<double>& normsGaussSeidel) {
	ofstream dataFile("data.txt", ios::out | ios::trunc);
    if (dataFile.is_open()) {
		dataFile << "Iteration\tJacobi Norm\tGauss-Seidel Norm\n";
		size_t maxIterations = max(normsJacobi.size(), normsGaussSeidel.size());
        for (size_t i = 0; i < maxIterations; ++i) {
			dataFile << i + 1 << "\t";
			if (i < normsJacobi.size())
				dataFile << normsJacobi[i];
			dataFile << "\t";
			if (i < normsGaussSeidel.size())
				dataFile << normsGaussSeidel[i];
			dataFile << "\n";
		}
		dataFile.close();
	}
    else {
		cout << "Unable to open file";
		return;
	}

	// Call Python script
	system("python plot_norm.py data.txt");
}


void plot_times(vector<double>& timesJacobi, vector<double>& timesGaussSeidel, vector<double>& timesLU) {
	ofstream dataFile("data.txt", ios::out | ios::trunc);
	if (dataFile.is_open()) {
		dataFile << "Size\tJacobi Time\tGauss-Seidel Time\tLU Time\n";
		for (size_t i = 0; i < timesJacobi.size(); ++i) {
			dataFile << i + 1 << "\t" << timesJacobi[i] << "\t" << timesGaussSeidel[i] << "\t" << timesLU[i] << "\n";
		}
		dataFile.close();
	}
	else {
		cout << "Unable to open file";
		return;
	}

	// Call Python script
	system("python plot_time.py data.txt");
}


int main()
{
    int a1 = 5 + 6;
    int a2 = -1, a3 = -1;
    int f = 6;
    int N = 933;
    vector<int> sizes = { 100, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000 };
	vector<double> timesJacobi;
	vector<double> timesGaussSeidel;
	vector<double> timesLU;
    for (int size : sizes) {
		GenerateMatrix matrixes(size, a1, a2, a3, f);
		matrixes.generate();
		SolveEquation iteration_methods(&matrixes);
		
		iteration_methods.solveJacobi();
		timesJacobi.push_back(iteration_methods.getTime());

		iteration_methods.solveGaussSeidel();
		timesGaussSeidel.push_back(iteration_methods.getTime());
		
		iteration_methods.solveFactorizationLU();
		timesLU.push_back(iteration_methods.getTime());
	}
	plot_times(timesJacobi, timesGaussSeidel, timesLU);

    return 0;
}

//============================================================================
/// Name        : 2D_AdvectionDiffusion.h
/// Author      : Mason Biamonte
/// Version     : 1.0
/// Copyright   :
/// Description : This code numerically propagates a linear advective-diffusive system with
/// the Peaceman-Rachford Alternating Direction Implicit scheme in 2D using Neumann boundary conditions.
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <math.h>
#include <time.h>
#include <algorithm>
#include "CLHEP/Matrix/defs.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/GenMatrix.h"
#include "Eigen/Dense"

#include <blitz/array.h>

#define PI 3.1415926535897932385
#define qE 1.602176487E-19 ///fundamental charge (C)
#define mE 9.11E-31 ///electron mass (kg)
#define kB 1.3806503E-23 ///Boltzmann constant (J/K)
#define eps0 8.854187817620E-12 ///permitivity of free space (F/m)

using namespace std;
using namespace CLHEP;
using namespace Eigen;


void OpenICFile (ifstream &infile){

	infile.open("/home/mason/Analysis/Diffusion/Debug/InitialChargeDist_Alpha_Am241_1.txt",ifstream::in);
}

void ClearVec(HepVector &v){

	/**
	 * This function fills vector v with zeros.
	 */


	int size = v.num_row();

	for(int i = 0; i < size; i++){

		v[i] = 0.0;

	}

}

void ClearVec(VectorXf &v){

	int size = v.rows();

	for (int i = 0; i < size; i++) v(i) = 0.0;
}

double minimum (HepVector &v){

	double min = fabs(v[0]);
	int N = v.num_row();

	for (int i = 1; i < N; i++){

		if (fabs(v[i]) <= min) min = fabs(v[i]);

	}

	return min;

}

double minimum (VectorXf &v){

	double temp = 0.0;

	temp = v(0);

	double min = temp;
	int N = v.rows();

	temp = 0.0;

	for (int i = 1; i < N; i++){

		temp = v(i);

		if (temp <= min) min = temp;

		temp = 0.0;

	}

	return min;

}

double maximum (HepVector &v) {

	double max = fabs(v[0]);
	int N = v.num_row();

	for (int i = 1; i < N; i++){

		if (fabs(v[i]) >= max) max = fabs(v[i]);
	}

	return max;

}

double maximum (VectorXf &v) {

	double max = v(0);
	int N = v.rows();

	for (int i = 1; i < N; i++){

		if (v(i) >= max) max = v(i);
	}

	return max;

}

int maximum (VectorXi &v) {

	int max = v(0);
	int N = v.rows();

	for (int i = 1; i < N; i++){

		if (v(i) >= max) max = v(i);
	}

	return max;

}

int maxElement (VectorXi &v){

	int maxElem = 0;
	int max = v(0);
	int N = v.rows();

	for (int i = 1; i < N; i++){

		if (v(i) >= max) {

			max = v(i);
			maxElem = i;
		}

	}

		return maxElem;


}


void PrintHepVec (HepVector &v){

	/**
	 * This function prints the vector v.
	 */

	for (int i = 0; i < v.num_row(); i++)
	{
		printf("%.03d %.08f\n", i, v[i]);
	}

	cout << endl;

}

void PrintEigenVec (VectorXf &v){

	/**
	 * This function prints the vector v.
	 */

	int N = v.rows();

	for (int i = 0; i < N; i++)
	{
		printf("%.03d %.18f\n", i, v[i] );
	}

	cout << endl;

}

void PrintEigenVec (VectorXi &v){

	/**
	 * This function prints the vector v.
	 */

	int N = v.rows();

	for (int i = 0; i < N; i++)
	{
		printf("%.03d %.03d\n", i, v(i) );
	}

	cout << endl;

}

void PrintMatToFile (HepMatrix &u, ofstream * fd){

	fd->precision(10);

	for (int i = 0; i < u.num_row(); i++){

		for (int j = 0; j < u.num_col(); j++){

			if (u[i][j] != 0.0)	*fd << (double)(i)*1E-06 << '\t' << (double)(j)*1E-06 << '\t' << u[i][j] << endl;

		}
	}
}

double SumVec (HepVector &v){

	double sum = 0.0;

	for (int i = 0; i < v.num_row(); i++){

		sum += v[i];
	}

	return sum;
}

int NumColData (ifstream &infile){

	string line;
	int col_counter = 0;

	if (infile.is_open()) {

		while (getline(infile,line)){

			col_counter++;
		}

		infile.close();
	}

	return col_counter;


}


void DataColumnToVector (ifstream &infile, int reqCol, int dPnts, HepVector &v)
{

	///This function grabs the n_th column from the data file in "infile" and stores it in vector v


	ifstream ifs ( "./InitialChargeDist_Alpha_Am241_1.txt", ifstream::in );

	int nCol = 5;

	double val;
	int col_count = 0;
	int row_count = 0;


	vector<double> tempVec;

	while (ifs.good()){

		ifs >> val;

		if (ifs.eof()) break;

		if (col_count % nCol == 0) col_count = 0;
		if (col_count == reqCol) tempVec.push_back(val);

		col_count++;

		if (col_count == nCol) row_count++;

		int numRows = (int) tempVec.size();

		for (int i = 0; i < numRows; i++) v[i] = tempVec[i];

	}



}

int NumRow (ifstream &infile, string fileName){


	string tab("\t");
	string newline("\n");

	int row_count = 0;

	infile.open(fileName.c_str());

	while (!(getline(infile,newline))){

		while (getline(infile,tab)){

				row_count++;
		}
	}

	infile.close();

	return row_count;

}



void DataToMatrix (string fileName, MatrixXf &outMat){

	ifstream ifs (fileName.c_str(), ifstream::in);

	string tab("\t");
	string newline("\n");

	int M = 2000;
	int row_count = 0;

	while (!(getline(ifs,newline))){

		if (getline(ifs,tab)) row_count++;
	}

	outMat.resize(M,row_count);





}



void DataColumnToVector (ifstream &infile, int reqCol, int dPnts, VectorXf &v)
{

	///This function grabs the n_th column from the data file in "infile" and stores it in vector v


	ifstream ifs ( "./InitialChargeDist_Alpha_Am241_1.txt", ifstream::in );

	int nCol = 5;

	double val;
	int col_count = 0;
	int row_count = 0;


	vector<double> tempVec;

	while (ifs.good()){

		ifs >> val;

		if (ifs.eof()) break;

		if (col_count % nCol == 0) col_count = 0;
		if (col_count == reqCol) tempVec.push_back(val);

		col_count++;

		if (col_count == nCol) row_count++;

		int numRows = (int) tempVec.size();

		for (int i = 0; i < numRows; i++) v(i) = tempVec[i];

	}

}


void diffVec (HepVector &v, HepVector &dv){

	int N = v.num_row();

	for (int i = 0; i < N-1; i++){

		dv[i] = fabs(v[i+1] - v[i]);
	}
}

void diffVec (VectorXf &v, VectorXf &dv){

	int N = v.rows();

	double temp = 0.0;

	for (int i = 0; i < N-1; i++){

		temp = v(i+1) - v(i);
		dv(i) = fabs(temp);

		temp = 0.0;
	}
}

void diffVec (vector<double> &v, vector<double> &dv){

	int N = v.size();

	double temp = 0.0;

	for (int i = 0; i < N-1; i++){

		temp = v[i+1] - v[i];
		dv[i] = fabs(temp);

		temp = 0.0;
	}
}

void dot (vector<double> &v, vector<double> &w, vector<double> &z){

	int N = v.size();

	for (int i = 0; i < N-1; i++){

		z[i] = v[i]*w[i];

	}
}

double Integrate (vector<double> &f, vector<double> &dx){

	///Integrates a discretized function using a non-uniform Simpson's rule

	int N = f.size();
	double sum = 0.0;

	sum += f[0]*dx[0] + f[N-1]*dx[N-1];

	for (int j = 0; j < N/2-1; j++) sum += 2.0*f[2*j]*dx[2*j];
	for (int j = 1; j < N/2; j++) sum += 4.0*f[2*j-1]*dx[2*j-1];

	sum /= 3.0;

	return sum;
}

double Integrate (VectorXf &f, VectorXf &dx){

	///Integrates a discretized function using a non-uniform Simpson's rule

	int N = f.rows();
	double sum = 0.0;

	sum += f(0)*dx(0) + f(N-1)*dx(N-1);

	for (int j = 0; j < N/2-1; j++) sum += 2.0*f(2*j)*dx(2*j);
	for (int j = 1; j < N/2; j++) sum += 4.0*f(2*j-1)*dx(2*j-1);

	sum /= 3.0;

	return sum;
}


double max (VectorXf &v){

	int N = v.rows();

	double temp = v(0);

	for (int i = 0; i < N; i++){


		if (v(i) >= temp ) temp = v(i);

	}

	return temp;
}


double min (VectorXf &v){

	int N = v.rows();

	double temp = v(0);

	for (int i = 0; i < N-1; i++){

		if (v(i) <= temp) temp = v(i);
	}

	return temp;
}



void TransformXCoord (HepVector &Xdata, HepVector &Xsim, double deltaGD){

	int J = Xsim.num_row();

	int m  = 0;

	for (int j = 0; j < J; j++){

		m = J-j-1;

		Xsim[m] = 1E-06*Xdata[j] - deltaGD;

		Xdata[j] = 0.0;

	}
}


void TransformYCoord (HepVector &Ydata, HepVector &Ysim, int numPix, double pitch){

	int J = Ysim.num_row();
	double halfway = (double) (numPix*pitch)/2.0;

	for (int j = 0; j < J; j++){

		Ysim[j] = halfway + Ydata[j]*1.0E-06;

	}
}


void TransformXCoord (VectorXf &Xdata, VectorXf &Xsim, double deltaGD){

	int J = Xsim.rows();

	int m  = 0;

	for (int j = 0; j < J; j++){

		m = J-j-1;

		Xsim(m) = (1E-06*Xdata(j) - deltaGD);
		Xdata(j) = 0.0;

	}
}


void TransformYCoord (VectorXf &Ydata, VectorXf &Ysim, int numPix, double pitch){

	int J = Ysim.rows();

	double halfway = (double) (numPix*pitch)/2.0;

	for (int j = 0; j < J; j++){

		Ysim(j) = halfway + Ydata(j)*1.0E-06;

	}
}

void EnergyLossToParticles (HepVector &dEdata, HepVector &EHPsim, double EHP_Si){

	int J = EHPsim.num_row();

	for (int j = 0; j < J; j++){

		EHPsim[j] = floor((dEdata[j]/EHP_Si)*1E03);

	}
}

void EnergyLossToParticles (VectorXf &dEdata, VectorXf &EHPsim, double EHP_Si){

	int J = EHPsim.rows();

	double temp = 0.0;

	for (int j = 0; j < J; j++){

		temp = dEdata(j)/EHP_Si;
		temp *= 1.0E03;

		EHPsim(j) = round(temp);

		temp = 0.0;

	}
}

void ParticlesToDensity (HepVector &EHPsim, HepVector &u0sim, double dA){

	int N = EHPsim.num_row();

	for (int j = 0; j < N; j++){

		u0sim[j] = EHPsim[j]/(dA);
	}
}

void ParticlesToDensity (VectorXf &EHPsim, VectorXf &u0sim, double dA){

	int N = EHPsim.rows();

	for (int j = 0; j < N; j++){

		u0sim(j) = EHPsim(j)/(dA);
	}
}

void FillInitialMatrix (HepVector &xIntIndex, HepVector &yIntIndex, HepVector &u0sim, HepMatrix &u_0){

	int N = xIntIndex.num_row();

	for (int i = 0; i < N; i++)	u_0[xIntIndex[i]][yIntIndex[i]] += u0sim[i];

}

void FillInitialMatrix (VectorXf &Xsim, VectorXf &Ysim, VectorXf &X_muGrid, VectorXf &Y_muGrid, VectorXf &EHPsim, MatrixXf &u_0, double dx, double dy){

	int dPnts = Xsim.rows();

	double xmin = 0.0;
	double xmax = 0.0;
	double ymin = 0.0;
	double ymax = 0.0;

	double num;
	double dA = dx*dy;

	int i_min = (round(min(Xsim)/dx) - 1) - 200;
	int i_max = (round(max(Xsim)/dx) + 1) - 200;

	int j_min = round(min(Ysim)/dy) - 1;
	int j_max = round(max(Ysim)/dy) + 1;

	int count = 0;



	for (int i = i_min; i < i_max; i++){
		for (int j = j_min; j < j_max; j++){

			xmin = (X_muGrid(i)) - dx/2.0;
			xmax = (X_muGrid(i)) + dx/2.0;

			ymin = Y_muGrid(j) - dy/2.0;
			ymax = Y_muGrid(j) + dy/2.0;


			for (int k = 0; k < dPnts; k++) {


				if (Xsim(k) > xmin && Xsim(k) <= xmax && Ysim(k) > ymin && Ysim(k) <= ymax) num += EHPsim(k);


			}

			if (num >= 1.0) u_0(i,j) = num/dA;

			num = 0.0;

			count++;

		}
	}









}

void FillInitialChargeIndex (VectorXi &XsimIndex, VectorXi &YsimIndex, int M, int J){

	int Nx = XsimIndex.rows();
	int Ny = YsimIndex.rows();


	int k = -(Ny-1)/2;

	for (int i = 0; i < Nx; i++) {

		XsimIndex(i) = M - Nx + i;
	}


	for (int j = 0; j < Ny; j++) {

		YsimIndex(j) = (J-1)/2 + k;
		k++;
	}



}


double Average (HepVector &v){

	int N = v.num_row();

	double result = 0.0;

	for (int i = 0; i < N; i++){

		result += v[i];

	}

	result *= (1.0)/N;

	return result;
}

double Average (VectorXf &v){

	int N = v.rows();

	double result = 0.0;

	for (int i = 0; i < N; i++){

		result += v[i];

	}

	result *= (1.0)/N;

	return result;
}

double GlobalChargeDensity(HepMatrix &u, double dx, double dy){

	double result = 0.0;
	double dA = dx*dy;

	int M = u.num_row();
	int N = u.num_col();

	double totalArea = (M-1)*(N-1)*dA;

	for (int i = 0; i < M; i++){

		for (int j = 0; j < N; j++){

			result += u[i][j]*dA;

		}
	}

	result /= totalArea;

	return result;
}

double GlobalChargeDensity(MatrixXf &u, double dx, double dy){

	double result = 0.0;
	double dA = dx*dy;

	int M = u.rows();
	int N = u.cols();

	double totalArea = (M-1)*(N-1)*dA;

	for (int i = 0; i < M; i++){

		for (int j = 0; j < N; j++){

			result += u(i,j)*dA;

		}
	}

	result /= totalArea;

	return result;
}


double dArea (int n, double D_n, double u0, double dx, double dy, double t){

	return n*n*D_n*u0*t*dx*dy*(dx*dx+dy*dy+48.0*D_n*t)
			*(1.0/(sqrt(dx*dy*(dx*dx+24.0*D_n*t)*(dy*dy+24.0*D_n*t))))*t;
}

int dGridPoints (int n, double D_n, double u0, double dx, double dy, double t){

	return abs(ceil(dArea(n,D_n,u0,dx,dy,t)/(dx*dy)));
}

void ClearMat (HepMatrix &u){

	int M = u.num_row();
	int J = u.num_col();

	for (int i = 0; i < M; i++){

		for (int j = 0; j < J; j++){

			u[i][j] = 0.0;

		}
	}


}

void ClearMat (MatrixXf &u){

	int M = u.rows();
	int J = u.cols();

	for (int i = 0; i < M; i++){

		for (int j = 0; j < J; j++){

			u(i,j) = 0.0;

		}
	}


}



HepMatrix AddRingsSolnMat (int k, HepMatrix &inMat){

	int M0 = inMat.num_row();
	int J0 = inMat.num_col();

	int M = M0 + 2*k;
	int J = J0 + 2*k;

	HepMatrix outMat(M,J,0.0);

	for (int i = 0; i < M0; i++){

		for (int j = 0; j < J0; j++){

			outMat[i+k][j+k] = inMat[i][j];

		}
	}

	return outMat;


}




void TriDiagMatVecMult (HepMatrix &B, HepVector &v, HepVector &y){

	int size = v.num_row();
	int j_start = 0;
	int j_end = 0;

	for (int i = 0; i < size; i++)
	{
		j_start = max(1,i-1);
		j_end = min(size,i+1);

		for (int j = j_start; j <= j_end; j++)
		{
			y[i] += B[i][j]*v[j];
		}
	}

}

void TriDiagSolve (HepVector &a, HepVector &b, HepVector &c, HepVector &r, HepVector &v, int &tCount){

	double num = b[0];
	int n = b.num_row();

	HepVector gamma(n);

	if (num == 0.0) cout << "Error 1 in TriDiagSolve" << endl;

	v[0] = r[0]/num;


	for (int j = 1; j < n; j++){

		gamma[j] = c[j-1]/num;
		num = b[j]-a[j]*gamma[j];

		if (num == 0.0) cout << "Error 2 in TriDiagSolve" << endl;

		v[j] = (r[j]-a[j]*v[j-1])/num;
	}

	for (int j = (n-1); j>=0; j--){

		v[j] -= gamma[j+1]*v[j+1];
	}

	ClearVec(gamma);//output to file

	//tCount++;

}

void TriDiagSolve2D (HepVector &a, HepVector &b, HepVector &c, HepVector &r, HepVector &v){

	double num = b[0];
	int n = b.num_row();

	HepVector gamma(n);

	if (num == 0.0) cout << "Error 1 in TriDiagSolve" << endl;

	v[0] = r[0]/num;


	for (int j = 1; j < n; j++){

		gamma[j] = c[j-1]/num;
		num = b[j]-a[j]*gamma[j];

		if (num == 0.0) cout << "Error 2 in TriDiagSolve" << endl;

		v[j] = (r[j]-a[j]*v[j-1])/num;
	}

	for (int j = (n-1); j>=0; j--){

		v[j] -= gamma[j+1]*v[j+1];
	}

	ClearVec(gamma);//output to file

}

void TriDiagSolve2D (VectorXf &a, VectorXf &b, VectorXf &c, VectorXf &r, VectorXf &v){

	double num = b(0);
	int n = b.rows();

	//cout << "dim(b) = " << n << endl;
	//cout << "dim(v) = " << v.rows() << endl;

	double temp = 0.0;

	VectorXf gamma(n);

	if (num == 0.0) cout << "Error 1 in TriDiagSolve" << endl;

	v(0) = r(0)/num;


	for (int j = 1; j < n; j++){

		gamma(j) = c(j-1)/num;
		num = b(j)-a(j)*gamma(j);

		if (num == 0.0) cout << "Error 2 in TriDiagSolve" << endl;

		v(j) = (r(j)-a(j)*v(j-1))/num;
	}

	for (int j = (n-2); j>=0; j--){

		temp = gamma(j+1)*v(j+1);

		v(j) -= temp;
	}

	gamma.setZero(n);

}


void PrintChargeDensity (double dx, double dy, HepVector &Xsim, HepVector &Ysim, HepMatrix &u, ofstream * fd){

	int M = u.num_row();
	int N = u.num_col();

	double x0 = minimum(Xsim);
	double y0 = minimum(Ysim);

	for (int i = 0; i < M; i++){

		for (int j = 0; j < N; j++){

			*fd << x0 + (i*dx) << '\t' << y0 + (j*dy) << '\t' << u[i][j] << endl;

		}
	}

}

void PrintEHPData (VectorXf &X, VectorXf &Y, VectorXf &EHP, ofstream * fd){

	int dPnts = X.rows();


	for (int i = 0; i < dPnts; i++) {

		*fd << Y(i) << '\t' << X(i) << '\t' << EHP(i) << endl;

	}
}

void PrintGrid (VectorXf &X, VectorXf &Y, ofstream * fd){

	int M = X.rows();
	int J = Y.rows();

	for (int i = 0; i < M; i++){

		*fd  << X(i) << endl;

	}

	*fd << endl;


	for (int j = 0; j < J; j++){

		*fd << Y(j) <<  endl;

	}

}

void PrintChargeDensity (double dx, double dy, VectorXf &X, VectorXf &Y, MatrixXf &u, bool carrierType, ofstream * fd){

	int M = u.rows();
	int J = u.cols();

	//double dA = dx*dy;

	int count = 0;

	for (int i = 0; i < M; i++){

		for (int j = 0; j < J; j++){

			if (carrierType == true) {

				if (u(i,j) > 0.0 && count%2==0) *fd << X(i) << '\t' << Y(j) << '\t' << -qE*u(i,j) << endl;

			}

			else {

				if (u(i,j) > 0.0 && count%2==0) *fd << X(i) << '\t' << Y(j) << '\t' << qE*u(i,j) << endl;

			}


			count++;

		}

	}

}

void PrintNumberDensity (double dx, double dy, VectorXf &X, VectorXf &Y, MatrixXf &u, bool carrierType, ofstream * fd){

	int M = u.rows();
	int J = u.cols();

	double dA = dx*dy;

	double numEPS = 1.0/dA;

	int count = 0;

	for (int i = 0; i < M; i++){

		for (int j = 0; j < J; j++){


			if (u(i,j) >= numEPS && count%2==0) *fd << X(i) << '\t' << Y(j) << '\t' << u(i,j) << endl;

			count++;

		}

	}

}

void PrintVector (VectorXf &v, double dy, ofstream * fd){

	int J = v.rows();
	fd->precision(14);

	for (int j = 0; j < J; j++){

		*fd << (double) (j*dy) << "\t" << v(j) << endl;
	}
}

void PrintVector (vector<int> &v, double dt, ofstream *fd){

	int J = v.size();
	fd->precision(12);

	for (int j = 0; j < J; j++){

		*fd << (double) (j*dt) << "\t" << v[j] << endl;

	}
}

void PrintVector (vector<double> &v, double dt, ofstream *fd){

	int J = v.size();
	fd->precision(12);

	for (int j = 0; j < J; j++){

		*fd << (double) (j*dt) << "\t" << v[j] << endl;

	}

}


void PrintMatrix (MatrixXf &A, ofstream * fd){

	int M = A.rows();
	int J = A.cols();

	for (int i = 0; i < M; i++){

		for (int j = 0; j < J; j++){

			*fd << A(i,j) << "\t";


		}

		*fd << endl;
	}

}


void FillAppliedField_x (HepMatrix &Ex, double E0){

	/**
	 * This function fills the matrix Ey with the y-component
	 * of the electric field resulting from the applied voltage
	 * for the 2D algorithm.
	 */

	int M = Ex.num_row();
	int J = Ex.num_col();

	for (int m = 0; m < M; m++){

		for (int j = 0; j < J; j++){

			Ex[m][j] = E0;

		}
	}



}

void FillAppliedField_x (MatrixXf &Ex, double E0, double dx){

	/**
	 * This function fills the matrix Ey with the y-component
	 * of the electric field resulting from the applied voltage
	 * for the 2D algorithm.
	 */

	int M = Ex.rows();
	int J = Ex.cols();

	double x = 0.0;

	for (int m = 0; m < M; m++){

		x = m*dx;

		for (int j = 0; j < J; j++){

			if (x <= 0.6E-03) Ex(m,j) = E0;

		}
	}



}

void FillAppliedField_y (HepMatrix &Ey, double E0){

	/**
	 * This function fills the matrix Ey with the y-component
	 * of the electric field resulting from the applied voltage
	 * for the 2D algorithm.
	 */

	int M = Ey.num_row();
	int J = Ey.num_col();

	for (int m = 0; m < M; m++){

		for (int j = 0; j < J; j++){

			Ey[m][j] = 0.0;

		}
	}



}

void FillAppliedField_y (MatrixXf &Ey, double E0){

	/**
	 * This function fills the matrix Ey with the y-component
	 * of the electric field resulting from the applied voltage
	 * for the 2D algorithm.
	 */

	int M = Ey.rows();
	int J = Ey.cols();

	for (int m = 0; m < M; m++){

		for (int j = 0; j < J; j++){

			Ey(m,j) = 0.0;

		}
	}



}

void FillDiffField_x (HepMatrix &Ex, HepMatrix &diffEx) {

	/**
	 * This function fills the matrix diffEy with the y-component
	 * of the difference in the electric field vector between
	 * neighboring grid points for the purposes of differentiation.
	 */

	///HARD-CODED FOR CONSTANT APPLIED BIAS FIELD!!!

	int M = Ex.num_row();
	int J = Ex.num_col();

	for (int m = 0; m < M; m++) {

		for (int j = 0; j < J; j++){

			diffEx[m][j] = 0.0;
		}
	}




}

void FillDiffField_x (MatrixXf &Ex, MatrixXf &diffEx) {

	/**
	 * This function fills the matrix diffEy with the y-component
	 * of the difference in the electric field vector between
	 * neighboring grid points for the purposes of differentiation.
	 */

	///HARD-CODED FOR CONSTANT APPLIED BIAS FIELD!!!

	int M = Ex.rows();
	int J = Ex.cols();

	for (int m = 0; m < M; m++) {

		for (int j = 0; j < J; j++){

			diffEx(m,j) = 0.0;
		}
	}




}

void FillDiffField_y (HepMatrix &Ey, HepMatrix &diffEy) {

	/**
	 * This function fills the matrix diffEy with the y-component
	 * of the difference in the electric field vector between
	 * neighboring grid points for the purposes of differentiation.
	 */

	///HARD-CODED FOR CONSTANT APPLIED BIAS FIELD!!!

	int M = Ey.num_row();
	int J = Ey.num_col();

	for (int m = 0; m < M; m++) {

		for (int j = 0; j < J; j++){

			diffEy[m][j] = 0.0;

		}
	}



}

void FillDiffField_y (MatrixXf &Ey, MatrixXf &diffEy) {

	/**
	 * This function fills the matrix diffEy with the y-component
	 * of the difference in the electric field vector between
	 * neighboring grid points for the purposes of differentiation.
	 */

	///HARD-CODED FOR CONSTANT APPLIED BIAS FIELD!!!

	int M = Ey.rows();
	int J = Ey.cols();

	for (int m = 0; m < M; m++) {

		for (int j = 0; j < J; j++){

			diffEy(m,j) = 0.0;

		}
	}



}



void CNVectorsPopulate2Dx (HepVector &aAx, HepVector &bAx, HepVector &cAx,
		double rx) {

	/**
	 * This function populates the matrix on the left-hand side of the
	 * tri-diagonal system of equations to be solved in the first step of
	 * the 1D Crank-Nicolson scheme.
	 */


	int M_int = bAx.num_row();

	bAx[0] = 1.0 + 0.5*rx;
	bAx[M_int-1] = bAx[0];

	for (int m = 1; m < M_int-1; m++){

		bAx[m] = 1.0 + rx;
		aAx[m] = -0.5*rx;
	}

	for (int m = 0; m < M_int-1; m++) cAx[m] = -0.5*rx;

}


void CNVectorsPopulate2Dy (HepVector &aAx, HepVector &bAx, HepVector &cAx,
		double ry, bool schemeType, bool carrierType) {


	int J_int = bAx.num_row();

	if (carrierType == true){

		bAx[0] = 1.0 + 0.5*ry;
		bAx[J_int-1] = 1.0 + 0.5*ry;

		for (int m = 1; m < J_int-1; m++){

			bAx[m] = 1.0 + ry;
			aAx[m] = -0.5*ry;
		}

		for (int m = 0; m < J_int-1; m++) cAx[m] = -0.5*ry;

	}

	else {


	}

}

void CNVectorsPopulate2D_Implicit (VectorXf &aA, VectorXf &bA, VectorXf &cA, double r, bool carrierType) {

	int N = bA.rows();

	///implicit scheme

	if (carrierType == true){

		///scheme for electrons

		bA(0) = 1.0 + 0.5*r;
		bA(N-1) = bA(0);

		for (int l = 1; l < N-1; l++){

			bA(l) = 1.0 + r;
			aA(l) = -0.5*r;

		}

		for (int l = 0; l < N-1; l++) cA(l) = -0.5*r;

		/*

		if (component == true){ ///X-direction

			bA(0) = 1.0 + r + 0.5*s*E(1,j);
			bA(N-1) = 1.0 + r - 0.5*s*E(N-2,j);

			for (int l = 1; l < N-1; l++){

				bA(l) = 1.0 + 2.0*r;
				aA(l) =  -r + 0.5*s*E(l,j);
			}

			for (int l = 0; l < N-1; l++) cA(l) = -r - 0.5*s*E(l,j);

		}

		else { ///Y-direction

			bA(0) = 1.0 + r + 0.5*s*E(j,1);
			bA(N-1) = 1.0 + r - 0.5*s*E(j,N-2);

			for (int l = 1; l < N-1; l++){

				bA(l) = 1.0 + 2.0*r;
				aA(l) =  -r + 0.5*s*E(j,l);
			}

			for (int l = 0; l < N-1; l++) cA(l) = -r - 0.5*s*E(j,l);



		}

		 */


	}

	else {

		///scheme for holes


	}




}


void CNVectorsPopulate2D_SemiImplicit (VectorXf &aA, VectorXf &bA, VectorXf &cA,
		double r, bool carrierType){


	int N = bA.rows();


	///semi-implicit scheme

	if (carrierType == true){

		///scheme for electrons

		bA(0) = 1.0 + 0.5*r;
		bA(N-1) = 1.0 + 0.5*r;

		for (int l = 1; l < N-1; l++){

			bA(l) = 1.0 + r;
			aA(l) =  -0.5*r;
		}

		for (int l = 0; l < N-1; l++) cA(l) = -0.5*r;



	}

	else {

		///scheme for holes


	}



}



void CNVectorsPopulate2Dx (VectorXf &aAx, VectorXf &bAx, VectorXf &cAx,
		double rx, bool schemeType, bool carrierType) {

	/**
	 * This function populates the matrix on the left-hand side of the
	 * tri-diagonal system of equations to be solved in the first step of
	 * the 1D Crank-Nicolson scheme.
	 */

	int M_int = bAx.rows();

	if (schemeType == true){

		///implicit diffusion scheme

		if (carrierType == true){

			///scheme for electrons

			bAx(0) = 2.0*(1.0 + 0.5*rx);
			bAx(M_int-1) = 2.0*(1.0 + 0.5*rx);

			for (int m = 1; m < M_int-1; m++){

				bAx(m) = 2.0*(1.0 + rx);
				aAx(m) =  -rx;

			}

			for (int m = 0; m < M_int-1; m++) cAx(m) = -rx;


		}

		else {

			///scheme for holes



		}




	}

	else {

		///semi-implicit diffusion scheme

		if (carrierType == true){

			///scheme for electrons

			bAx(0) = 1.0 + 0.5*rx;
			bAx(M_int-1) = 1.0 + 0.5*rx;

			for (int m = 1; m < M_int-1; m++){

				bAx(m) = 1.0 + rx;
				aAx(m) =  -0.5*rx;

			}

			for (int m = 0; m < M_int-1; m++) cAx(m) = -0.5*rx;

		}

		else {

			///scheme for holes


		}

	}

}



void CNVectorsPopulate2Dy (VectorXf &aAy, VectorXf &bAy, VectorXf &cAy,
		double ry, bool schemeType, bool carrierType) {

	int J_int = bAy.rows();

	if (schemeType == true){

		///implicit diffusion scheme

		if (carrierType == true){

			///scheme for electrons

			bAy(0) = 2.0*(1.0 + 0.5*ry);
			bAy(J_int-1) = 2.0*(1.0 + 0.5*ry);

			for (int j = 1; j < J_int-1; j++){

				bAy(j) = 2.0*(1.0 + ry);
				aAy(j) =  -ry;

			}

			for (int j = 0; j < J_int-1; j++) cAy(j) = -ry;

		}

		else {

			///scheme for holes




		}




	}

	else {

		///semi-implicit scheme

		if (carrierType == true){

			///scheme for electrons

			bAy(0) = 1.0 + 0.5*ry;
			bAy(J_int-1) = 1.0 + 0.5*ry;

			for (int j = 1; j < J_int-1; j++){

				bAy(j) = 1.0 + ry;
				aAy(j) = -0.5*ry;
			}

			for (int j = 0; j < J_int-1; j++) cAy(j) = -0.5*ry;

		}

		else {

			///scheme for holes




		}

	}

}



void RHSMatPopulate2D(HepMatrix &B, double r, double s, HepMatrix &E, HepMatrix &diffE, int p){

	/**
	 * This function populates the matrix on the right-hand side of the
	 * tri-diagonal system of equations to be solved in the first step of
	 * the 2D Crank-Nicolson scheme.
	 */

	int K_int = B.num_row();

	B[0][0] = 1.0 - r + s*(diffE[1][p]-E[1][p]);
	B[0][1] = r + s*E[1][p];
	B[K_int-1][K_int-2] = r - s*E[K_int-1][p];
	B[K_int-1][K_int-1] = 1 - r + s*(diffE[K_int-1][p]+E[K_int-1][p]);

	for (int i = 1; i < K_int-1; i++){

		for (int k = 0; k < K_int; k++){

			if (k == i) B[i][k] = 1.0 - 2.0*r + s*(diffE[i][p]);
			else if (k == i-1) B[i][k] = r - s*E[i][p];
			else if (k == i+1) B[i][k] = r + s*E[i][p];
			else B[i][k] = 0.0;
		}

	}

}

void RHSMatPopulate2D_Implicit(MatrixXf &B, double s, MatrixXf &E, int p, bool component, bool carrierType){

	int K_int = B.rows();

	if (carrierType == true){

		if (component == true) {

		B(0,0) = 1.0 - s*E(1,p);
		B(0,1) = s*E(1,p);
		B(K_int-1,K_int-2) = -s*E(K_int-1,p);
		B(K_int-1,K_int-1) = 1.0 + s*E(K_int - 1, p);

		for (int i = 1; i < K_int-1; i++){

			for (int k = 0; k < K_int; k++){

				if (k == i) B(i,k) = 1.0;
				else if (k == i-1) B(i,k) = -s*E(i,p);
				else if (k == i+1) B(i,k) =  s*E(i,p);
				else B(i,k) = 0.0;
			}

		}

		}

		else {


			B(0,0) = 1.0 - s*E(p,1);
			B(0,1) = s*E(p,1);
			B(K_int-1,K_int-2) = -s*E(p,K_int-1);
			B(K_int-1,K_int-1) = 1.0 + s*E(p,K_int - 1);

			for (int i = 1; i < K_int-1; i++){

				for (int k = 0; k < K_int; k++){

					if (k == i) B(i,k) = 1.0;
					else if (k == i-1) B(i,k) = -s*E(p,i);
					else if (k == i+1) B(i,k) =  s*E(p,i);
					else B(i,k) = 0.0;
				}

			}


		}

	}

}

void RHSMatPopulate2D_SemiImplicit(MatrixXf &B, double r, double s, MatrixXf &E, MatrixXf &diffE, int p, bool component, bool carrierType){

	/**
	 * This function populates the matrix on the right-hand side of the
	 * tri-diagonal system of equations to be solved in the 2D Crank-Nicolson scheme.
	 */

	int K_int = B.rows();



	if (carrierType == true){

		if (component == true){

		B(0,0) = 1.0 + s*(diffE(1,p)-E(1,p));
		B(0,1) = s*E(1,p);
		B(K_int-1,K_int-2) = -s*E(K_int-1,p);
		B(K_int-1,K_int-1) = 1.0 + s*(diffE(K_int-1,p)+E(K_int-1,p));


		for (int i = 1; i < K_int-1; i++){

			for (int k = 0; k < K_int; k++){

				if (k == i) B(i,k) = 1.0 + s*diffE(i,p);
				else if (k == i-1) B(i,k) = -s*E(i,p);
				else if (k == i+1) B(i,k) =  s*E(i,p);
				else B(i,k) = 0.0;
			}

		}

		}

		else {

			B(0,0) = 1.0 + s*(diffE(p,1)-E(p,1));
			B(0,1) = s*E(p,1);
			B(K_int-1,K_int-2) = -s*E(p,K_int-1);
			B(K_int-1,K_int-1) = 1.0 + s*(diffE(p,K_int-1)+E(p,K_int-1));


			for (int i = 1; i < K_int-1; i++){

				for (int k = 0; k < K_int; k++){

					if (k == i) B(i,k) = 1.0 + s*diffE(p,i);
					else if (k == i-1) B(i,k) = -s*E(p,i);
					else if (k == i+1) B(i,k) =  s*E(p,i);
					else B(i,k) = 0.0;
				}

			}



		}

	}

	else {


	}


}


void RHSMatPopulate2Dx(MatrixXf &B, double r, double s, MatrixXf &E, MatrixXf &diffE, int p, bool schemeType, bool carrierType){

	/**
	 * This function populates the matrix on the right-hand side of the
	 * tri-diagonal system of equations to be solved in the first step of
	 * the 2D Crank-Nicolson scheme.
	 */

	int K_int = B.rows();

	if (schemeType == true){

		///implicit scheme

		if (carrierType == true){

			///scheme for electrons

			B(0,0) = 1.0 + s*(diffE(1,p) - E(1,p));
			B(0,1) = s*E(1,p);

			B(K_int-1,K_int-2) = -s*E(K_int-1,p);
			B(K_int-1,K_int-1) = 1.0 + s*(diffE(K_int-1,p)+E(K_int-1,p));

			for (int i = 1; i < K_int-1; i++){

				for (int k = 0; k < K_int; k++){

					if (k == i) B(i,k) = 1.0 + s*diffE(i,p);
					else if (k == i-1) B(i,k) = -s*E(i,p);
					else if (k == i+1) B(i,k) =  s*E(i,p);
					else B(i,k) = 0.0;
				}

			}

		}

		else {

			///scheme for holes
		}


	}

	else {

		///semi-implicit scheme

		if (carrierType == true){

			///scheme for electrons

			B(0,0) = 1.0 - r + s*(diffE(1,p)-E(1,p));
			B(0,1) = r + s*E(1,p);

			B(K_int-1,K_int-2) = r - s*E(K_int-1,p);
			B(K_int-1,K_int-1) = 1.0 - r + s*(diffE(K_int-1,p)+E(K_int-1,p));


			for (int i = 1; i < K_int-1; i++){

				for (int k = 0; k < K_int; k++){

					if (k == i) B(i,k) = 1.0 - 2.0*r + s*diffE(i,p);
					else if (k == i-1) B(i,k) = r - s*E(i,p);
					else if (k == i+1) B(i,k) = r + s*E(i,p);
					else B(i,k) = 0.0;
				}

			}

		}

		else {

			///scheme for holes
		}


	}

}

void RHSMatPopulate2Dy(MatrixXf &B, double r, double s, MatrixXf &E, MatrixXf &diffE, int p, bool schemeType, bool carrierType){

	/**
	 * This function populates the matrix on the right-hand side of the
	 * tri-diagonal system of equations to be solved in the first step of
	 * the 2D Crank-Nicolson scheme.
	 */

	int K_int = B.rows();

	if (schemeType == true){

		///implicit scheme

		if (carrierType == true){

			///scheme for electrons


		}

		else {

			///scheme for holes


		}

	}

	else {

		///semi-implicit scheme

		if (carrierType == true){

			///scheme for electrons

			B(0,0) = 1.0 - r + s*(diffE(p,1)-E(p,1));
			B(0,1) = r + s*E(p,1);
			B(K_int-1,K_int-2) = r - s*E(p,K_int-1);
			B(K_int-1,K_int-1) = 1.0 - r + s*(diffE(p,K_int-1)+E(p,K_int-1));


			for (int i = 1; i < K_int-1; i++){

				for (int k = 0; k < K_int; k++){

					if (k == i) B(i,k) = 1.0 - 2.0*r + s*diffE(p,i);
					else if (k == i-1) B(i,k) = r - s*E(p,i);
					else if (k == i+1) B(i,k) = r + s*E(p,i);
					else B(i,k) = 0.0;
				}

			}

		}

		else {

			///scheme for holes


		}


	}
}






void AddRingsSolnMat (int k, MatrixXf &inMat){

	int M0 = inMat.rows();
	int J0 = inMat.cols();

	int M = M0 + 2*k;
	int J = J0 + 2*k;

	MatrixXf outMat(M,J);
	outMat.setZero(M,J);

	for (int i = 0; i < M0; i++){

		for (int j = 0; j < J0; j++){

			outMat(i+k,j+k) = inMat(i,j);

		}
	}

	inMat.setZero(M0,J0);
	inMat.resize(M,J);
	inMat.setZero(M,J);

	inMat = outMat;

	outMat.setZero(M,J);

}

HepMatrix AddRingsMat (int k, HepMatrix &inMat){

	int M0 = inMat.num_row();
	int J0 = inMat.num_col();

	int M = M0 + 2*k;
	int J = J0 + 2*k;

	HepMatrix outMat(M,J,0.0);

	for (int i = 0; i < M; i++){

		for (int j = 0; j < J; j++){

			outMat[i][j] = 0.0;

		}
	}

	return outMat;


}

void AddRingsMat (int k, MatrixXf &inMat){

	int M0 = inMat.rows();
	int J0 = inMat.cols();

	int M = M0 + 2*k;
	int J = J0 + 2*k;

	MatrixXf outMat(M,J);
	outMat.setZero(M,J);

	inMat.resize(M,J);

	inMat = outMat;


}

void AddRowsMat (int k, MatrixXf &inMat){

	///Call in conjunction with AddRingsVec for Xsim

	int M0 = inMat.rows();
	int J0 = inMat.cols();

	int M = M0 + k;

	MatrixXf outMat(M,J0);
	outMat.setZero(M,J0);

	inMat.resize(M,J0);

	inMat = outMat;

}

void AddColsMat (int k, MatrixXf &inMat){

	///Call in conjunction with AddRingsVec for Ysim

	int M0 = inMat.rows();
	int J0 = inMat.cols();

	int J = J0 + k;

	MatrixXf outMat(M0,J);
	outMat.setZero(M0,J);

	inMat.resize(M0,J);

	inMat = outMat;

}



HepMatrix AddRingsFieldMat (int k, HepMatrix &inMat, double val){

	int M0 = inMat.num_row();
	int J0 = inMat.num_col();

	int M = M0 + 2*k;
	int J = J0 + 2*k;

	HepMatrix outMat(M,J,0.0);

	for (int i = 0; i < M; i++){

		for (int j = 0; j < J; j++){

			outMat[i][j] = val;

		}
	}

	return outMat;


}

void AddRingsPosVec (int k, VectorXf &v_in, double dz){

	int M0 = v_in.rows();
	int M = M0 + 2*k;

	VectorXf v_out(M);

	v_out.setZero(M);


	for (int j = 0; j < k; j++){

		v_out(j) = v_in(0) - (k-j)*dz;
		v_out(M0+k+j) = v_in(M0-1) + (j+1)*dz;
	}

	for (int j = 0; j < M0; j++){

		v_out(j+k) = v_in(j);
	}

	v_in.resize(M);

	v_in = v_out;

}

HepMatrix AddRingsVec (int k, HepVector &v_in){

	int M0 = v_in.num_row();
	int M = M0 + 2*k;

	HepVector v_out(M,0.0);

	for (int i = 0 ; i < M0; i++){

		v_out[i+k] = 0.0;

	}

	return v_out;

}

void AddRingsVec (int k, VectorXf &v_in){

	int M0 = v_in.rows();
	int M = M0 + 2*k;

	VectorXf v_out(M);
	v_out.setZero(M);

	v_in.resize(M);

	v_in = v_out;


}

double MaxMatVal (HepMatrix &A){

	double result = A[0][0];

	int M = A.num_row();
	int J = A.num_col();

	for (int i = 0; i < M; i++){

		for (int j = 0; j < J; j++){

			if (A[i][j] >= result) result = A[i][j];

		}
	}

	return result;

}

double MaxMatVal (MatrixXf &A){

	double result = A(0,0);

	int M = A.rows();
	int J = A.cols();

	for (int i = 0; i < M; i++){

		for (int j = 0; j < J; j++){

			if (A(i,j) >= result) result = A(i,j);

		}
	}

	return result;

}

int NumGrdPntAdd (double u_max, int stdDev, double D, double dx, double dy, double dt){

	double dA = dx*dy;
	int result = 0;

	result = round((stdDev*stdDev*D*u_max*dt*(dx*dx + dy*dy + 48.0*D*dt))/(sqrt(dA*(dx*dx + 24.0*D*dt)*(dy*dy + 24.0*D*dt))));

	return result;

}

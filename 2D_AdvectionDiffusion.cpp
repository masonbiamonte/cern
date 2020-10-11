//========================================================================================================
/// Name        : 2D_AdvectionDiffusion.cpp
/// Author      : Mason Biamonte
/// Version     : 1.0
/// Copyright   :
/// Description : This code numerically propagates a linear advective-diffusive system with
/// the Peaceman-Rachford Alternating Direction Implicit scheme in 2D using Neumann boundary conditions.
/// The algorithm terminates when the charge density distribution has passed through 300 um of the silicon
/// medium, using a model for the charge collection mechanism of the Timepix device.
//========================================================================================================

#include "2D_AdvectionDiffusion.h"

using namespace CLHEP;
using namespace std;
using namespace Eigen;

//BZ_USING_NAMESPACE(blitz)



int main(int argc, char ** argv)
{


	///RUN TEMPLATE: ./InitialChargeDistRead (bias_type) (voltage)
	// (time_step_neg_power) (interaction_type) (carrier_type)

	double V;  ///applied voltage [V]
	double dt; ///time step [s]

	double T; ///ambient Terature [K]


	//int timeSteps = 10;

	bool interactionType;    ///true for Coulomb interactions
	bool carrierType; ///true for electrons

	if (string(argv[1]) == "f") V = atof(argv[2]);
	else if (string(argv[1]) == "b") V = -atof(argv[2]);

	//dt = 1.0E-11;
	//interactionType = true;
	//carrierType = true;

	dt = 1.0*pow(10.0,-atof(argv[3]));
	if (string(argv[4]) == "free") interactionType = false;
	else if (string(argv[4]) == "coulomb") interactionType = true;

	if (string(argv[5]) == "electron") carrierType = true;
	else if (string(argv[5]) == "hole") carrierType = false;

	T = atof(argv[6]);


	///---------INPUT FILE-----------///

	ifstream IC_inFile;
	OpenICFile(IC_inFile);

	///------------------------------///

	///-------OUTPUT FILES----------///


	ofstream * ICDensitySimMat = new ofstream;
	ofstream * ICDensityDataMat = new ofstream;
	ofstream * ICEHPSim = new ofstream;
	ofstream * InitialGrid = new ofstream;
	ofstream * CoulombField = new ofstream;

	///-------FREE---------///

	ofstream * outFileDiffusionTerm_Free = new ofstream;
	ofstream * outFileAdvectionTerm_Free = new ofstream;
	ofstream * outFileChargeFluct_Free = new ofstream;
	ofstream * outFileChargeLeak_Free = new ofstream;
	ofstream * outFileTotalCharge_Free = new ofstream;
	ofstream * outFileCloudSize_Free = new ofstream;
	ofstream * outFileTOT_Free = new ofstream;
	ofstream * outFileChargeCollect_Free = new ofstream;
	ofstream * outFileFINAL_Free = new ofstream;
	ofstream * outFileINTERMEDIATE_Free = new ofstream;
	ofstream * outFileCCX_Trajectory_Free = new ofstream;
	ofstream * outFileCollectionTime_Free = new ofstream;
	ofstream * outFileLinChrgDens_Free = new ofstream;
	ofstream * outFileDeltaXIntg_Free = new ofstream;
	ofstream * outFileData_Free = new ofstream;

	///------COULOMB-------///

	ofstream * outFileDiffusionTerm_Coulomb = new ofstream;
	ofstream * outFileAdvectionTerm_Coulomb = new ofstream;
	ofstream * outFileChargeFluct_Coulomb = new ofstream;
	ofstream * outFileChargeLeak_Coulomb = new ofstream;
	ofstream * outFileTotalCharge_Coulomb = new ofstream;
	ofstream * outFileCloudSize_Coulomb = new ofstream;
	ofstream * outFileTOT_Coulomb = new ofstream;
	ofstream * outFileChargeCollect_Coulomb = new ofstream;
	ofstream * outFileFINAL_Coulomb = new ofstream;
	ofstream * outFileINTERMEDIATE_Coulomb = new ofstream;
	ofstream * outFileCCX_Trajectory_Coulomb = new ofstream;
	ofstream * outFileCollectionTime_Coulomb = new ofstream;
	ofstream * outFileLinChrgDens_Coulomb  = new ofstream;
	ofstream * outFileDeltaXIntg_Coulomb = new ofstream;
	ofstream * outFileData_Coulomb = new ofstream;

	vector<ofstream *> arrivalTimeSliceFiles;

	//int N = 64;
	// blitz::Array<double,3> A(N,N,N), B(N,N,N);





	///------------------------------///

	///---------------INITIAL PARAMETERS-------------///


	int tCount = 0;
	//int tCountFinal = 0;

	//double t = 0.0;
	///total propagation time

	double dt_half = 0.5*dt;
	///half time-step [s]

	double dx = 1.0E-06;
	/// desired spatial step in the X direction

	double dy = 1.0E-06;
	/// desired spatial step in the Y direction

	double dA = dx * dy;

	double totalCharge0 = 0.0;
	double totalChargeF = 0.0;

	double totalCharge = 0.0;
	double chargeCenterX = 0.0;

	int numPix;

	///Two extra pixels of silicon
	///for extra diffusion space

	if (V == 87.9) numPix = 8;
	if (V == 50.0) numPix = 9;
	if (V == 30.4) numPix = 10;
	if (V == 14.8) numPix = 10;
	if (V == 6.0) numPix = 9;


	int dPnts = NumColData(IC_inFile);
	///# of entries in initial
	///charge distribution file

	OpenICFile(IC_inFile);


	double depth = 0.6*1E-03;
	///depth of computational domain [m]

	double bottom = 0.2*1E-03;
	///value of the depth at the
	///bottom of the computational domain [m]

	double topGEANT4 = 0.95E-03;
	///valude of depth at the top
	///of the GEANT4 computational domain [m]

	//double thickness = 0.3E-03;
	///thickness of the sensor [m]

	double E0 = V/depth;
	///magnitude of applied electric field

	double pitch = 55.0E-06;
	///width of pixel [m]

	double EHP_Si = 3.66;
	///effective energy required
	/// to create an EHP in Si [eV]

	double epsR_Si = 11.68;
	///relative permittivity of silicon

	//double T = 300.0;

	double A_Si_n = 1.43E+09;
	///electron mobility coefficient
	///parameter for silicon near
	///room Terature [cm^2*K^gamma*V^-1*s^-1]

	double gamma_Si_n = 2.42;
	///electron mobility Terature
	///power parameter for silicon near room Terature

	double mu_n = A_Si_n * pow(T, -gamma_Si_n) * 1.0E-04;
	///electron mobility in silicon in SI units (m^2/V*s)
	///Terature dependence comes from acoustic phonon scattering

	double D_n = (mu_n*kB*T / qE);
	///diffusion constant for electrons in (m^2/s)

	double xLength = depth - bottom;
	double yLength = (double) (numPix)*pitch;

	double xSteps = xLength/dx;
	double ySteps = yLength/dy;

	int M = floor(xSteps) + 21; ///number of gridpoints in the X direction with a buffer for diffusion at the top of the computational domain
	int J = floor(ySteps) + 1;  ///number of gridpoints in the Y direction, adjusted by 1 to be an odd number

	int M_int = M - 2;
	int J_int = J - 2;


	int depthIndex = 101;


	double deltaGD = topGEANT4 - depth;
	///length difference between height
	///of GEANT4 and Diffusion computational domains [m]

	///------------------------------///

	///-------------INTIAL VECTORS/MATRICES-----------------///

	//-------Eigen---------//

	VectorXf Xdata(dPnts);
	VectorXf Ydata(dPnts);
	VectorXf dEdata(dPnts);

	VectorXf Xsim(dPnts);
	VectorXf Ysim(dPnts);

	VectorXf dXsim(dPnts);
	VectorXf dYsim(dPnts);

	VectorXf EHPsim(dPnts);
	VectorXf u0sim(dPnts);

	VectorXf u0sim_mu;

	VectorXf TOTcounter_n (J);
	VectorXf TOTcounter_p (J);

	VectorXf ChargeCounter_n (J);
	VectorXf ChargeCounter_p (J);

	MatrixXf u_0(M,J);

	MatrixXf u_n(M,J);
	MatrixXf u_p(M,J);




	///------------------------------///

	///-----------------------------///

	DataColumnToVector(IC_inFile, 1, dPnts, Ydata);
	OpenICFile(IC_inFile);

	DataColumnToVector(IC_inFile, 2, dPnts, Xdata);
	OpenICFile(IC_inFile);

	DataColumnToVector(IC_inFile, 3, dPnts, dEdata);

	TransformXCoord(Xdata, Xsim, deltaGD);
	TransformYCoord(Ydata, Ysim, numPix, pitch);

	EnergyLossToParticles(dEdata, EHPsim, EHP_Si);

	diffVec(Xsim, dXsim);
	diffVec(Ysim, dYsim);


	double rx_n = D_n*dt_half / (dx*dx);
	double ry_n = D_n*dt_half / (dy*dy);

	double sx_n = (mu_n * dt_half) / dx;
	double sy_n = (mu_n * dt_half) / dy;

	double threshE = 3.66;
	///threshold energy value [eV]

	double EPS = qE/dA; ///one carrier per area element [C/m^2]
	double threshDens = (threshE/EHP_Si)*EPS; ///charge density

	double renormCharge = 0.0;


	///-----------------------------///


	VectorXf X_muGrid(M);
	VectorXf Y_muGrid(J);

	for (int i = 0; i < M; i++) X_muGrid(i) = (i*dx) + bottom;
	for (int j = 0; j < J; j++) Y_muGrid(j) = j*dy;

	FillInitialMatrix(Xsim,Ysim,X_muGrid,Y_muGrid,EHPsim,u_0,dx,dy);

	u_n = u_0; ///sets electron solution matrix to initial distribution obtained from GEANT4
	u_p = u_0;  ///sets hole solution matrix to initial distribution obtained from GEANT4



	for (int i = 0; i < M; i++){
		for (int j = 0; j < J; j++){

			if (carrierType == true) {

				totalCharge0 += -u_n(i,j)*dA*qE;
				totalChargeF += -u_n(i,j)*dA*qE;

				chargeCenterX += -u_n(i,j)*dA*(i*dx)*qE;
			}

			else {

				totalCharge0 += u_p(i,j)*dA*qE;
				totalChargeF += u_p(i,j)*dA*qE;

				chargeCenterX += u_p(i,j)*dA*(i*dx)*qE;
			}



		}
	}

	chargeCenterX /= totalCharge0;
	chargeCenterX += bottom;


	///---CN Left-Hand Side Matrix Vectors---///

	///--------Electrons---------///


	VectorXf aAx_n(M_int);
	VectorXf bAx_n(M_int);
	VectorXf cAx_n(M_int);

	VectorXf bx_n(M_int);
	///resultant right-hand side vector
	///in first sub-step of ADI scheme for electrons

	VectorXf aAy_n(J_int);
	VectorXf bAy_n(J_int);
	VectorXf cAy_n(J_int);

	VectorXf by_n(J_int);
	///resultant right-hand side vector
	/// in second sub-step of ADI scheme for electrons


	//-------------------*/

	///----------Holes-----------///

	//-------------------*/

	//------Eigen--------//


	VectorXf aAx_p(M_int);
	VectorXf bAx_p(M_int);
	VectorXf cAx_p(M_int);

	VectorXf bx_p(M_int);
	///resultant right-hand side vector
	///in first sub-step of ADI scheme for holes

	VectorXf aAy_p(J_int);
	VectorXf bAy_p(J_int);
	VectorXf cAy_p(J_int);

	VectorXf by_p(J_int);
	///resultant right-hand side vector
	///in second sub-step of ADI scheme for holes



	//-------------------*/


	MatrixXf uINT_n(M_int, J_int);
	///interior region solution for electrons

	MatrixXf uINT_p(M_int, J_int);
	///interior region solution for holes

	VectorXf uINTx_n(M_int);
	///intermediate column vector for first sub-step

	VectorXf uINTx_n_prev(M_int);
	VectorXf uINTx_n_curr(M_int);
	VectorXf uINTx_n_next(M_int);

	VectorXf uINTy_n(J_int);
	///intermediate row vector for second sub-step

	VectorXf uINTy_n_prev(J_int);
	VectorXf uINTy_n_curr(J_int);
	VectorXf uINTy_n_next(J_int);

	MatrixXf Bx_n(M_int, M_int);
	///right-hand side matrix

	MatrixXf By_n(J_int, J_int);
	///right-hand side matrix_


	MatrixXf Ex(M, J);
	///longitudinal component of electric field [V/m]

	MatrixXf Ey(M, J);
	///tranverse component of electric field [V/m]

	MatrixXf diffEx(M, J);
	MatrixXf diffEy(M, J);

	MatrixXf AppliedField_x(M, J);
	MatrixXf AppliedField_y(M, J);

	MatrixXf CoulombField_x(M, J);
	MatrixXf CoulombField_y(M, J);

	FillAppliedField_x(AppliedField_x, E0, dx);
	FillAppliedField_y(AppliedField_y, E0);

	time_t start;
	time_t end;

	clock();
	time(&start);

	double totalTime = 0.0;

	bool schemeType = true;
	///true for Peaceman-Rachford implicit scheme


	int distCount = 0;

	double x = 0.0;
	double y = 0.0;
	double r = 0.0;
	//double r_max = 0.0;
	double theta = 0.0;

	vector<double> distCoord_x;
	vector<double> distCoord_y;



	///------------------TRACKERS/COUNTERS-----------------------------///


	int threshCount = 0;
	//double collectionTime = 0.0;
	///time at which threshold is
	///first exceeded by charge cloud
	///at measurement depth

	VectorXf collectionTime(J);
	///stores values of the collection time
	///as a function of the transverse length


	vector<int> CldSz_Coulomb, CldSz_Free;
	///store the number of grid points the charge
	///density cloud occupies at each time step
	///(above threshold)

	vector<double> ChrgLk_Coulomb, ChrgLk_Free;
	///store the global charge over time

	vector<double> ChrgFluct_Coulomb, ChrgFluct_Free;
	///stores fluctuations in global charge

	vector<double> totalChrg_Coulomb, totalChrg_Free;
	///stores value of global charge

	vector<double> CCX_Free, CCX_Coulomb;
	///store the depth of the center of charge

	vector<double> CCXdelta_Free, CCXdelta_Coulomb;

	MatrixXf linChrgDensX_Free(2000,J), linChrgDensX_Coulomb(2000,J);
	///store the linear charge density in the direction of drift

	cout << "J = " << J << endl;




	vector<double> DiffusionTerm, AdvectionTerm;
	///stores the magnitude of the contribution to time rate of change
	///of the number density from diffusion and advection, respectively

	double diffusion;
	double advection;
	///placeholders for values inserted in vectors above

	//double globalChargeLoss;
	///placeholder for charge lost through domain boundaries

	VectorXf timeDurationCount(J);

	VectorXf tA(J); ///vector of time indices denoting arrival of charge cluster
	VectorXf tB(J); ///vector of time indices denoting departure of charge cluster

	VectorXi xIntDim(J); ///array of integers specifying length over which the charge integration takes place

	VectorXf bottomVec(J);
	VectorXf topVec(J);
	///monitor charge loss at top and bottom of domain

	VectorXf leftVec(M);
	VectorXf rightVec(M);
	///monitor charge loss along right and left of domain


	VectorXf measureVec(J);

	int passDepthCheck = 0;
	///determines existence of
	///charge cloud above or below
	///measurement depth

	bool depthReach = false;
	///boolean that is false before
	///the charge distribution reaches
	///the measurement depth

	///-----------------------------------------------------------------///

	///----------------OUTPUT FILENAME STRING I/O--------------///


	ostringstream timeConvert;
	ostringstream voltageConvert;
	ostringstream threshConvert;

	voltageConvert << V;
	threshConvert << threshE;

	string timeStamp ("_tStep=");
	string timeString;

	string voltageStamp ("_V=");
	string voltageString;

	string threshStamp ("_thresh=");
	string threshString;

	string outFileFormat(".txt");

	threshString = threshConvert.str();
	threshStamp.append(threshString);
	threshStamp.append("eV");

	voltageString = voltageConvert.str();
	voltageStamp.append(voltageString);
	voltageStamp.append(threshStamp);
	voltageStamp.append(outFileFormat);

	string simMat ("ICDensitySimMat");
	string dataMat ("ICDensityDataMat");
	string EHPSim ("ICEHPCountSim");
	string IGrid ("InitialGrid");

	string Diffusion_Free("DiffusionTerm_Free");
	string Advection_Free("AdvectionTerm_Free");
	string Final_Free ("2D_DiffusionFINAL_Free");
	string Intermediate_Free("2D_DiffusionINTERMEDIATE_Free");
	string CFluct_Free("ChargeFluct_Free");
	string CLeak_Free ("ChargeLeak_Free");
	string totalCharge_Free("totalCharge_Free");
	string Cloud_Free("CloudSize_Free");
	string TOT_Free("TOT_Free");
	string ChargeCollect_Free("ChargeCollect_Free");
	string TrajectoryX_Free("TrajectoryX_Free");
	string CollectionTime_Free ("CollectionTime_Free");
	string LinChrgDensMat_Free ("LinChrgDensMat_Free");
	string DeltaXIntg_Free ("DeltaXIntg_Free");
	string Data_Free("Data_Free");

	string Diffusion_Coulomb("DiffusionTerm_Coulomb");
	string Advection_Coulomb("AdvectionTerm_Coulomb");
	string Final_Coulomb ("2D_DiffusionFINAL_Coulomb");
	string Intermediate_Coulomb("2D_DiffusionINTERMEDIATE_Coulomb");
	string CFluct_Coulomb("ChargeFluct_Coulomb");
	string CLeak_Coulomb ("ChargeLeak_Coulomb");
	string totalCharge_Coulomb ("totalCharge_Coulomb");
	string Cloud_Coulomb("CloudSize_Coulomb");
	string TOT_Coulomb("TOT_Coulomb");
	string ChargeCollect_Coulomb ("ChargeCollect_Coulomb");
	string TrajectoryX_Coulomb("TrajectoryX_Coulomb");
	string CollectionTime_Coulomb ("CollectionTime_Coulomb");
	string LinChrgDensMat_Coulomb ("LinChrgDensMat_Coulomb");
	string DeltaXIntg_Coulomb ("DeltaXIntg_Coulomb");
	string Data_Coulomb("Data_Coulomb");

	string CField("CoulombField");

	vector<string> outputFileNames;
	vector<string> outputFileNames_Free;
	vector<string> outputFileNames_Coulomb;

	outputFileNames.push_back(simMat);
	outputFileNames.push_back(dataMat);
	outputFileNames.push_back(EHPSim);
	outputFileNames.push_back(IGrid);
	outputFileNames.push_back(CField);

	outputFileNames_Free.push_back(Diffusion_Free);
	outputFileNames_Free.push_back(Advection_Free);
	outputFileNames_Free.push_back(Final_Free);
	outputFileNames_Free.push_back(CFluct_Free);
	outputFileNames_Free.push_back(CLeak_Free);
	outputFileNames_Free.push_back(totalCharge_Free);
	outputFileNames_Free.push_back(Cloud_Free);
	outputFileNames_Free.push_back(TOT_Free);
	outputFileNames_Free.push_back(ChargeCollect_Free);
	outputFileNames_Free.push_back(TrajectoryX_Free);
	outputFileNames_Free.push_back(CollectionTime_Free);
	outputFileNames_Free.push_back(LinChrgDensMat_Free);
	outputFileNames_Free.push_back(DeltaXIntg_Free);
	outputFileNames_Free.push_back(Data_Free);
	outputFileNames_Free.push_back(Intermediate_Free);

	outputFileNames_Coulomb.push_back(Diffusion_Coulomb);
	outputFileNames_Coulomb.push_back(Advection_Coulomb);
	outputFileNames_Coulomb.push_back(Final_Coulomb);
	outputFileNames_Coulomb.push_back(CFluct_Coulomb);
	outputFileNames_Coulomb.push_back(CLeak_Coulomb);
	outputFileNames_Coulomb.push_back(totalCharge_Coulomb);
	outputFileNames_Coulomb.push_back(Cloud_Coulomb);
	outputFileNames_Coulomb.push_back(TOT_Coulomb);
	outputFileNames_Coulomb.push_back(ChargeCollect_Coulomb);
	outputFileNames_Coulomb.push_back(TrajectoryX_Coulomb);
	outputFileNames_Coulomb.push_back(CollectionTime_Coulomb);
	outputFileNames_Coulomb.push_back(LinChrgDensMat_Coulomb);
	outputFileNames_Coulomb.push_back(DeltaXIntg_Coulomb);
	outputFileNames_Coulomb.push_back(Data_Coulomb);
	outputFileNames_Coulomb.push_back(Intermediate_Coulomb);

	for (int i = 0; i < (signed int) (outputFileNames.size()); i++) outputFileNames[i].append(voltageStamp);
	for (int i = 0; i < (signed int) (outputFileNames_Free.size()); i++) outputFileNames_Free[i].append(voltageStamp);
	for (int i = 0; i < (signed int) (outputFileNames_Coulomb.size()); i++) outputFileNames_Coulomb[i].append(voltageStamp);


	///---------------INTERACTION INDEPENDENT OUTPUT FILES-------------------///

	ICDensitySimMat->                open((char *) outputFileNames[0].c_str());
	ICDensityDataMat->               open((char *) outputFileNames[1].c_str());
	ICEHPSim->                       open((char *) outputFileNames[2].c_str());
	InitialGrid ->                   open((char *) outputFileNames[3].c_str());
	CoulombField->                   open((char *) outputFileNames[4].c_str());



	if (interactionType == true){

		outFileDiffusionTerm_Coulomb->   open((char *) outputFileNames_Coulomb[0].c_str());
		outFileAdvectionTerm_Coulomb->   open((char *) outputFileNames_Coulomb[1].c_str());
		outFileFINAL_Coulomb->           open((char *) outputFileNames_Coulomb[2].c_str());
		outFileChargeFluct_Coulomb->     open((char *) outputFileNames_Coulomb[3].c_str());
		outFileChargeLeak_Coulomb->      open((char *) outputFileNames_Coulomb[4].c_str());
		outFileTotalCharge_Coulomb->     open((char *) outputFileNames_Coulomb[5].c_str());
		outFileCloudSize_Coulomb->       open((char *) outputFileNames_Coulomb[6].c_str());
		outFileTOT_Coulomb->             open((char *) outputFileNames_Coulomb[7].c_str());
		outFileChargeCollect_Coulomb->   open((char *) outputFileNames_Coulomb[8].c_str());
		outFileCCX_Trajectory_Coulomb->  open((char *) outputFileNames_Coulomb[9].c_str());
		outFileCollectionTime_Coulomb->  open((char *) outputFileNames_Coulomb[10].c_str());
		outFileLinChrgDens_Coulomb->     open((char *) outputFileNames_Coulomb[11].c_str());
		outFileData_Coulomb->            open((char *) outputFileNames_Coulomb[12].c_str());

		outFileData_Coulomb->precision(20);

		*outFileData_Coulomb << "Threshold Energy: " << threshE << " eV" << endl;
		*outFileData_Coulomb << "Threshold Density: " << threshDens << " C/m^2" << endl;
		*outFileData_Coulomb << "Initial Total Charge: " << totalCharge0 << " C" << endl;
		*outFileData_Coulomb << "Initial Depth of Center of Charge: " << chargeCenterX << " m" << endl;
		*outFileData_Coulomb << "Applied Bias Voltage : V = " << V <<  "V" << endl;
		*outFileData_Coulomb << "Temperature: T = " << T << " K" << endl;

	}

	else {

		outFileDiffusionTerm_Free->      open((char *) outputFileNames_Free[0].c_str());
		outFileAdvectionTerm_Free->      open((char *) outputFileNames_Free[1].c_str());
		outFileFINAL_Free->              open((char *) outputFileNames_Free[2].c_str());
		outFileChargeFluct_Free->        open((char *) outputFileNames_Free[3].c_str());
		outFileChargeLeak_Free->         open((char *) outputFileNames_Free[4].c_str());
		outFileTotalCharge_Free->        open((char *) outputFileNames_Free[5].c_str());
		outFileCloudSize_Free->          open((char *) outputFileNames_Free[6].c_str());
		outFileTOT_Free->                open((char *) outputFileNames_Free[7].c_str());
		outFileChargeCollect_Free->      open((char *) outputFileNames_Free[8].c_str());
		outFileCCX_Trajectory_Free->     open((char *) outputFileNames_Free[9].c_str());
		outFileCollectionTime_Free->     open((char *) outputFileNames_Free[10].c_str());
		outFileLinChrgDens_Free->        open((char *) outputFileNames_Free[11].c_str());
		outFileData_Free->               open((char *) outputFileNames_Free[12].c_str());

		outFileData_Free->precision(20);

		*outFileData_Free << "Threshold Energy: " << threshE << " eV" << endl;
		*outFileData_Free << "Threshold Density: " << threshDens << " C/m^2" << endl;
		*outFileData_Free << "Initial Total Charge: " << totalCharge0 << " C" << endl;
		*outFileData_Free << "Initial Depth of Center of Charge: " << chargeCenterX << " m" << endl;
		*outFileData_Free << "Applied Bias Voltage : V = " << V <<  "V" << endl;
		*outFileData_Free << "Temperature: T = " << T << " K" << endl;

	}


	cout << "Threshold Density: " << threshDens << " C/m^2" << endl;
	cout << "Initial Total Charge: " << totalCharge0 << " C" << endl;
	cout << "Initial Depth of Center of Charge: " << chargeCenterX << " m" << endl;
	cout << endl;


	PrintEHPData(Xsim,Ysim,EHPsim,ICEHPSim);
	PrintGrid(X_muGrid,Y_muGrid,InitialGrid);
	PrintChargeDensity(dx,dy,X_muGrid,Y_muGrid,u_0,carrierType,ICDensitySimMat);

	chargeCenterX = 0.0;




	///-----------------TIME STEP ALGORITHM----------------------------///


	if (schemeType == true)  {

		CNVectorsPopulate2D_Implicit(aAx_n, bAx_n, cAx_n, rx_n, carrierType);
		CNVectorsPopulate2D_Implicit(aAy_n, bAy_n, cAy_n, ry_n, carrierType);

	}

	else {

		CNVectorsPopulate2D_SemiImplicit(aAx_n, bAx_n, cAx_n, rx_n, carrierType);
		CNVectorsPopulate2D_SemiImplicit(aAy_n, bAy_n, cAy_n, ry_n, carrierType);

	}


	do {

		/*

		timeStamp.clear();
		timeStamp= "_tStep=";

		timeConvert << tCount;
		timeString = timeConvert.str();
		timeStamp.append(timeString);
		timeStamp.append(outFileFormat);

		 */

		///for now, we are not producing all of the output plots of final distribution as it crosses the measurement plane
		///we need to figure out how to produce all of these objects on the fly -> vector of ofstream objects?


		passDepthCheck = 0;

		bool arrivalPrev = false;
		bool arrivalCurr = false;


		///booleans that monitor the
		///arrival of the charge cloud
		///at the measurement depth



		///-------STEP #1: Solves J_int (M_int) x (M_int) tridiagonal systems.------///

		///---loop over the columns of the solution grid for the first sub-step---///

		if (interactionType == true){

			distCount = 0;
			passDepthCheck = 0;

			x = 0.0;
			y = 0.0;
			r = 0.0;
			//r_max = 2.0*sqrt(dx*dx+dy*dy); ///nearest-neighbor interactions
			//r_max = 50.0E-06;
			theta = 0.0;

			for (int i = 0; i < M; i++){
				for (int j = 0; j < J; j++){

					if (fabs((double)(u_n(i,j)*qE)) > EPS){

						distCoord_x.push_back(i);
						distCoord_y.push_back(j);
						distCount++;
					}

				}
			}


			for (int l = 0; l < distCount; l++){
				for (int k = 0; k < distCount; k++){

					if (k != l){

						x = dx*(distCoord_x[k]-distCoord_x[l]);
						y = dy*(distCoord_y[k]-distCoord_y[l]);

						r = sqrt(x*x+y*y);
						theta = atan2(x,y);

						//cout << "Theta = " << theta << endl;

						if (carrierType == true){

							CoulombField_x(distCoord_x[l],distCoord_y[l]) += -qE*dA*u_n(distCoord_x[k],distCoord_y[k])*sin(theta)/(4.0*PI*eps0*epsR_Si*r*r);
							CoulombField_y(distCoord_x[l],distCoord_y[l]) += -qE*dA*u_n(distCoord_x[k],distCoord_y[k])*cos(theta)/(4.0*PI*eps0*epsR_Si*r*r);

						}

						else {

							CoulombField_x(distCoord_x[l],distCoord_y[l]) += qE*dA*u_n(distCoord_x[k],distCoord_y[k])*sin(theta)/(4.0*PI*eps0*epsR_Si*r*r);
							CoulombField_y(distCoord_x[l],distCoord_y[l]) += qE*dA*u_n(distCoord_x[k],distCoord_y[k])*cos(theta)/(4.0*PI*eps0*epsR_Si*r*r);

						}



					}

				}
			}

		}

		Ex = AppliedField_x + CoulombField_x;
		Ey = AppliedField_y + CoulombField_y;

		FillDiffField_x(Ex, diffEx);
		FillDiffField_y(Ey, diffEy);


		for (int j = 0; j < J_int; j++)
		{

			if (schemeType == true) RHSMatPopulate2D_Implicit(Bx_n, sx_n, Ex, j, true, carrierType);
			else RHSMatPopulate2D_SemiImplicit(Bx_n, rx_n, sx_n, Ex, diffEx, j, true, carrierType);
			///right-hand side matrix multiplying previous time step solution

			if (j == 0) {

				for (int m = 0; m < M_int; m++)
				{

					uINTx_n_curr(m) = u_n(m + 1, j + 1);
					uINTx_n_next(m) = u_n(m + 1, j + 2);

					if (schemeType == true){

						if (interactionType == true){

							bx_n(m) = (1.0 - ry_n - 2.0*sy_n * Ey(m + 1, j + 1) - 2.0*mu_n*qE*uINTx_n_curr(m)/(eps0*epsR_Si*dy)) * uINTx_n_curr(m)
											+ (ry_n +  2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_next(m);

						}

						else {

							bx_n(m) = (1.0 - ry_n - 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_curr(m)
											+ (ry_n +  2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_next(m);

						}

					}

					else {

						bx_n(m) = (1.0 - 2.0 * ry_n + 2.0*sy_n * (diffEy(m + 1, j + 1)
								- Ey(m + 1, j + 1))) * uINTx_n_curr(m)
								+ (2.0*ry_n + 2.0*sy_n*Ey(m + 1, j + 1)) * uINTx_n_next(m);

					}

				}


				bx_n += Bx_n*uINTx_n_curr; ///application of the RHS matrix
				bx_n *= 0.5;


				TriDiagSolve2D(aAx_n, bAx_n, cAx_n, bx_n, uINTx_n);
				///updates the (j+1)_th = 1_th column of the solution grid; storing in uINTx_n

				for (int m = 0; m < M_int; m++)
				{

					u_n(m + 1, j + 1) = uINTx_n(m);
					///fills (j+1)_th column of solution grid with uINTx_n

					u_n(m + 1, j) = u_n(m + 1, j + 1);
					///Neumann boundary conditions along left edge of domain
				}

				///---Neumann Boundary Conditions---///

				///----endpoints of column j = 1---///

				u_n(0, j + 1) = u_n(1, j + 1);
				u_n(M - 1, j + 1) = u_n(M_int, j + 1);

				///-----left-side corners--------///

				u_n(0, j) = u_n(0, j + 1);
				u_n(M - 1, j) = u_n(M_int, j + 1);

			}

			else if (j == J_int - 1) {

				for (int m = 0; m < M_int; m++)
				{

					uINTx_n_prev(m) = u_n(m + 1, j);
					uINTx_n_curr(m) = u_n(m + 1, j + 1);

					if (schemeType == true){

						if (interactionType == true) {

							bx_n(m) = (ry_n - 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_prev(m)
											+ (1.0 - ry_n +  2.0*sy_n * Ey(m + 1, j + 1) - (2.0*qE/(eps0*epsR_Si*dy))*mu_n*uINTx_n(m)) * uINTx_n_curr(m);

						}

						else {

							bx_n(m) = (ry_n - 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_prev(m)
											+ (1.0 - ry_n + 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_curr(m);

						}

					}

					else {

						bx_n(m) = (2.0 * ry_n - 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_prev(m)
										+ (1.0 - 2.0 * ry_n + 2.0*sy_n * (diffEy(m + 1, j + 1)
												+ Ey(m + 1, j + 1))) * uINTx_n_curr(m);

					}

				}



				bx_n += Bx_n * uINTx_n_curr;
				bx_n *= 0.5;




				TriDiagSolve2D(aAx_n, bAx_n, cAx_n, bx_n, uINTx_n);
				///updates the (j+1)_th = 1_th column of the solution grid; storing in uINTx_n

				for (int m = 0; m < M_int; m++)
				{

					u_n(m + 1, j + 1) = uINTx_n(m);
					///fills (j+1)_th = (J-2)_th column of solution grid with uINTx_n

					u_n(m + 1, j + 2) = u_n(m + 1, j + 1);
					///Neumann boundary conditions along right edge of domain

				}

				///---Neumann Boundary Conditions---///

				///----endpoints of column j = J-2 ---///

				u_n(0, j + 1) = u_n(1, j + 1);
				u_n(M - 1, j + 1) = u_n(M_int, j + 1);

				///-----right-side corners--------///

				u_n(0, j + 2) = u_n(0, j + 1);
				u_n(M - 1, j + 2) = u_n(M - 1, j + 1);

			}

			else
			{

				for (int m = 0; m < M_int; m++)
				{

					uINTx_n_prev(m) = u_n(m + 1, j);
					uINTx_n_curr(m) = u_n(m + 1, j + 1);
					uINTx_n_next(m) = u_n(m + 1, j + 2);

					if (schemeType == true){

						if (interactionType == true){

							bx_n(m) = (ry_n - 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_prev(m)
											+ (1.0 - 2.0 * ry_n - (2.0*qE/(eps0*epsR_Si*dy))*mu_n*uINTx_n_curr(m)) * uINTx_n_curr(m)
											+ (ry_n +  2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_next(m);

						}

						else {

							bx_n(m) = (ry_n - 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_prev(m)
							    			+ (1.0 - 2.0 * ry_n) * uINTx_n_curr(m)
							    			+ (ry_n + 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_next(m);
						}

					}

					else {

						bx_n(m) = (2.0 * ry_n - 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_prev(m)
										+ (1.0 - 4.0 * ry_n + 2.0*sy_n* diffEy(m + 1, j + 1)) * uINTx_n_curr(m)
										+ (2.0 * ry_n + 2.0*sy_n * Ey(m + 1, j + 1)) * uINTx_n_next(m);

					}


				}


				bx_n += Bx_n * uINTx_n_curr;
				bx_n *= 0.5;


				TriDiagSolve2D(aAx_n, bAx_n, cAx_n, bx_n, uINTx_n);
				///updates the (j+1)_th = 2_th column of the solution grid; storing in uINTx_n

				for (int m = 0; m < M_int; m++) u_n(m + 1, j + 1) = uINTx_n(m);


				///---Neumann Boundary Conditions---///

				///----endpoints of column j+1 --///

				u_n(0, j + 1) = u_n(1, j + 1);
				u_n(M - 1, j + 1) = u_n(M_int, j + 1);

			}

		}



		for (int i = 1; i < M-1; i++){
			for (int j = 1; j < J-1; j++){

				diffusion += D_n*((u_n(i-1,j) - 2.0*u_n(i,j) + u_n(i+1,j))/(dx*dx)
						+ (u_n(i,j-1) - 2.0*u_n(i,j) + u_n(i,j+1))/(dy*dy));

				advection += mu_n*((u_n(i+1,j)-u_n(i-1,j))*Ex(i,j)/(2.0*dx)
						+ (u_n(i,j+1)-u_n(i,j-1))*Ey(i,j)/(2.0*dy)
						+ u_n(i,j)*((Ex(i+1,j)-Ex(i-1,j))/(2.0*dx)
								+ (Ex(i,j+1)-Ex(i,j-1))/(2.0*dy)));

			}
		}

		diffusion /= (double) (M*J);
		advection /= (double) (M*J);

		diffusion = abs(diffusion);
		advection = abs(advection);

		DiffusionTerm.push_back(diffusion);
		AdvectionTerm.push_back(advection);



		ClearMat(CoulombField_x);
		ClearMat(CoulombField_y);
		ClearMat(Ex);
		ClearMat(Ey);
		ClearMat(diffEx);
		ClearMat(diffEy);

		distCoord_x.clear();
		distCoord_y.clear();

		distCount = 0;

		diffusion = 0.0;
		advection = 0.0;


		///-------STEP #2: Set up M_int (J_int) x (J_int) tridiagonal systems; Solve them------///

		///---loop over the rows of the solution grid for the second sub-step---///

		///--------Initialization----------///

		if (interactionType == true){

			distCount = 0;
			passDepthCheck = 0;

			x = 0.0;
			y = 0.0;
			r = 0.0;
			//r_max = 2.0*sqrt(dx*dx+dy*dy); ///two-neighbor interactions
			//r_max = 50.0E-06;
			theta = 0.0;

			for (int i = 0; i < M; i++){
				for (int j = 0; j < J; j++){

					if (fabs((double)(u_n(i,j))*qE) > EPS){

						distCoord_x.push_back(i);
						distCoord_y.push_back(j);
						distCount++;
					}

				}
			}


			for (int l = 0; l < distCount; l++){
				for (int k = 0; k < distCount; k++){

					if (k != l){

						x = dx*(distCoord_x[k]-distCoord_x[l]);
						y = dy*(distCoord_y[k]-distCoord_y[l]);

						r = sqrt(x*x+y*y);
						theta = atan2(x,y);

						//cout << "Theta = " << theta << endl;


						if (carrierType == true){

							CoulombField_x(distCoord_x[l],distCoord_y[l]) += -qE*dA*u_n(distCoord_x[k],distCoord_y[k])*sin(theta)/(4.0*PI*eps0*epsR_Si*r*r);
							CoulombField_y(distCoord_x[l],distCoord_y[l]) += -qE*dA*u_n(distCoord_x[k],distCoord_y[k])*cos(theta)/(4.0*PI*eps0*epsR_Si*r*r);

						}

						else {

							CoulombField_x(distCoord_x[l],distCoord_y[l]) += qE*dA*u_n(distCoord_x[k],distCoord_y[k])*sin(theta)/(4.0*PI*eps0*epsR_Si*r*r);
							CoulombField_y(distCoord_x[l],distCoord_y[l]) += qE*dA*u_n(distCoord_x[k],distCoord_y[k])*cos(theta)/(4.0*PI*eps0*epsR_Si*r*r);

						}


					}

				}
			}

		}

		Ex = AppliedField_x + CoulombField_x;
		Ey = AppliedField_y + CoulombField_y;

		FillDiffField_x(Ex, diffEx);
		FillDiffField_y(Ey, diffEy);



		for (int m = 0; m < M_int; m++)
		{

			if (schemeType == true) RHSMatPopulate2D_Implicit(By_n, sy_n, Ey, m, false, carrierType);
			else  RHSMatPopulate2D_SemiImplicit(By_n, ry_n, sy_n, Ey, diffEy, m, false, carrierType);
			///right-hand side matrix multiplying previous time sub-step solution


			if (m == 0)
			{

				for (int j = 0; j < J_int; j++)
				{

					uINTy_n_curr(j) = u_n(m + 1, j + 1);
					uINTy_n_next(j) = u_n(m + 2, j + 1);

					if (schemeType == true){

						if (interactionType == true){

							by_n(j) = (1.0 - rx_n - 2.0*sx_n * Ex(m + 1, j + 1) - 2.0*qE*mu_n*uINTy_n_curr(j)/(eps0*epsR_Si*dx)) * uINTy_n_curr(j)
											+ (rx_n + 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_next(m);

						}

						else {

							by_n(j) = (1.0 - rx_n - 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_curr(j)
											+ (rx_n + 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_next(j);

						}

					}

					else {

						by_n(j) = (1.0 - 2.0 * rx_n + 2.0*sx_n * (diffEx(m + 1, j + 1)
								- Ex(m + 1, j + 1))) * uINTy_n_curr(j)
								+ (2.0 * rx_n + 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_next(j);

					}

				}


				by_n += By_n * uINTy_n_curr;
				by_n *= 0.5;



				TriDiagSolve2D(aAy_n, bAy_n, cAy_n, by_n, uINTy_n);
				///updates the (m+1)_th = 1_th row of the solution grid; storing in uINTy_n

				for (int j = 0; j < J_int; j++)
				{

					u_n(m + 1, j + 1) = uINTy_n(j);
					u_n(m, j + 1) = u_n(m + 1, j + 1);
				}

				///---Neumann Boundary Conditions---///

				///----endpoints of column m = 1---///

				u_n(m + 1, 0) = u_n(m + 1, 1);
				u_n(m + 1, J - 1) = u_n(m + 1, J_int);

				///-----bottom corners--------///

				u_n(m, 0) = u_n(m + 1, 0);
				u_n(m, J - 1) = u_n(m + 1, J - 1);

			}

			else if (m == M_int - 1)
			{

				for (int j = 0; j < J_int; j++)
				{

					uINTy_n_prev(j) = u_n(m, j + 1);
					uINTy_n_curr(j) = u_n(m + 1, j + 1);

					if (schemeType == true){

						if (interactionType == true){

							by_n(j) = (rx_n - 2.0*sx_n * Ex(m + 1, j + 1))* uINTy_n_prev(j)
											+ (1.0 - rx_n + 2.0*sx_n * Ex(m + 1, j + 1) - 2.0*qE*mu_n*uINTy_n_curr(j)/(eps0*epsR_Si*dx)) * uINTy_n_curr(j);

						}

						else {

							by_n(j) = (rx_n - 2.0*sx_n * Ex(m + 1, j + 1))* uINTy_n_prev(j)
											+ (1.0 - rx_n + 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_curr(j);

						}

					}

					else {

						by_n(j) = (2.0 * rx_n - 2.0*sx_n * Ex(m + 1, j + 1))* uINTy_n_prev(j)
				    					+ (1.0 - 2.0 * rx_n + 2.0*sx_n* (diffEx(m + 1, j + 1)
				    							+ Ex(m + 1, j + 1))) * uINTy_n_curr(j);

					}

				}


				by_n += By_n * uINTy_n_curr;
				by_n *= 0.5;


				TriDiagSolve2D(aAy_n, bAy_n, cAy_n, by_n, uINTy_n); ///updates the (m+1)_th row of the solution grid; storing in uINTx_n

				for (int j = 0; j < J_int; j++)
				{

					u_n(m + 1, j + 1) = uINTy_n(j); ///fills (m+1)_th row of solution grid with uINTx_n
					u_n(m + 2, j + 1) = u_n(m + 1, j + 1); ///Neumann boundary conditions along top edge of domain

				}

				///---Neumann Boundary Conditions---///

				///----endpoints of row m = M-2 ---///

				u_n(m + 1, 0) = u_n(m + 1, 1);
				u_n(m + 1, J - 1) = u_n(m + 1, J_int);

				///-----top corners--------///

				u_n(m + 2, 0) = u_n(m + 1, 0);
				u_n(m + 2, J - 1) = u_n(m + 1, J - 1);

			}

			else
			{

				for (int j = 0; j < J_int; j++)
				{

					uINTy_n_prev(j) = u_n(m, j + 1);
					uINTy_n_curr(j) = u_n(m + 1, j + 1);
					uINTy_n_next(j) = u_n(m + 2, j + 1);

					if (schemeType == true){

						if (interactionType == true){

							by_n(j) = (rx_n - 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_prev(j)
											+ (1.0 - 2.0 * rx_n - 2.0*qE*mu_n*uINTy_n_curr(j)/(eps0*epsR_Si*dx)) * uINTy_n_curr(j)
											+ (rx_n + 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_next(j);


						}

						else {

							by_n(j) = (rx_n - 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_prev(j)
											+ (1.0 - 2.0 * rx_n) * uINTy_n_curr(j)
											+ (rx_n + 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_next(j);

						}

					}

					else {

						by_n(j) = (2.0 * rx_n - 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_prev(j)
										+ (1.0 - 4.0 * rx_n + 2.0*sx_n * diffEx(m + 1, j + 1)) * uINTy_n_curr(j)
										+ (2.0 * rx_n + 2.0*sx_n * Ex(m + 1, j + 1)) * uINTy_n_next(j);

					}

				}


				by_n += By_n * uINTy_n_curr;
				by_n *= 0.5;


				TriDiagSolve2D(aAy_n, bAy_n, cAy_n, by_n, uINTy_n);
				///updates the (m+1)_th row of the solution grid; storing in uINTy_n

				for (int j = 0; j < J_int; j++)	u_n(m + 1, j + 1) = uINTy_n(j); ///fills (m+1)_th row interior of solution grid;


				///---Neumann Boundary Conditions---///

				///----endpoints of row m+1 --///

				u_n(m + 1, 0) = u_n(m + 1, 1);
				u_n(m + 1, J - 1) = u_n(m + 1, J_int);

			}

		}


		for (int i = 1; i < M-1; i++){
			for (int j = 1; j < J-1; j++){

				diffusion += D_n*((u_n(i-1,j) - 2.0*u_n(i,j) + u_n(i+1,j))/(dx*dx)
						+ (u_n(i,j-1) - 2.0*u_n(i,j) + u_n(i,j+1))/(dy*dy));

				advection += mu_n*((u_n(i+1,j)-u_n(i-1,j))*Ex(i,j)/(2.0*dx)
						+ (u_n(i,j+1)-u_n(i,j-1))*Ey(i,j)/(2.0*dy)
						+ u_n(i,j)*((Ex(i+1,j)-Ex(i-1,j))/(2.0*dx)
								+ (Ex(i,j+1)-Ex(i,j-1))/(2.0*dy)));

			}
		}

		diffusion /= (double) (M*J);
		advection /= (double) (M*J);

		diffusion = abs(diffusion);
		advection = abs(advection);

		DiffusionTerm.push_back(diffusion);
		AdvectionTerm.push_back(advection);


		ClearMat(CoulombField_x);
		ClearMat(CoulombField_y);
		ClearMat(Ex);
		ClearMat(Ey);
		ClearMat(diffEx);
		ClearMat(diffEy);

		distCoord_x.clear();
		distCoord_y.clear();

		distCount = 0;

		diffusion = 0.0;
		advection = 0.0;


		for (int j = 0; j < J; j++){

			if (u_n(depthIndex,j)*qE > threshDens){

				timeDurationCount(j)++;

				TOTcounter_n(j) += dt;
				ChargeCounter_n(j) += u_n(depthIndex,j)*dA*qE; ///understood to be negative charge

				if (interactionType == false) linChrgDensX_Free(timeDurationCount(j)-1,j) += u_n(depthIndex,j)*dy*qE;
				else linChrgDensX_Coulomb(timeDurationCount(j)-1,j) += u_n(depthIndex,j)*dy*qE;

				arrivalCurr = true;

				if (interactionType == true){

					PrintChargeDensity(dx,dy,X_muGrid,Y_muGrid,u_n,carrierType,outFileINTERMEDIATE_Coulomb);
					outFileINTERMEDIATE_Coulomb->close();

				}

				else {

					PrintChargeDensity(dx,dy,X_muGrid,Y_muGrid,u_n,carrierType,outFileINTERMEDIATE_Free);
					outFileINTERMEDIATE_Free->close();

				}

				if (arrivalCurr == true && arrivalPrev == false) {

					collectionTime(j) = (double) (dt*tCount);
					depthReach = true;

					tA(j) = tCount;
					tB(j) = tCount-1;

				}

				arrivalPrev = arrivalCurr;
				arrivalCurr = false;

				tB(j)++;

			}



		}

		for (int i = depthIndex; i < M; i++){
			for (int j = 0; j < J; j++){

				if (u_n(i,j)*qE > threshDens ) passDepthCheck++;
			}
		}

		CldSz_Free.push_back(passDepthCheck);


		for (int i = 0; i < M; i++){
			for (int j = 0; j < J; j++){

				if (u_n(i,j)*qE >= threshDens) threshCount++;

				if (carrierType == true) totalCharge += -qE*u_n(i,j)*dA;
				else totalCharge += qE*u_p(i,j)*dA;


			}
		}

		for (int i = 0; i < M; i++){

			leftVec(i) = u_n(i,1);
			rightVec(i) = u_n(i,M-2);

			//*if (abs((double) u_n(i, 1)*qE) > EPS)*/ totalChargeF += leftVec(i)*dA*qE;
			//*if (abs((double) u_n(i, J-2)*qE) > EPS)*/ totalChargeF += rightVec(i)*dA*qE;

		}

		for (int j = 0; j < J; j++){

			bottomVec(j) = u_n(1,j);
			topVec(j) = u_n(M-2,j);

			//*if (abs((double) u_n(1,j)*qE) > EPS)*/ totalChargeF += bottomVec(j)*dA*qE;
			//*if (abs((double) u_n(M-2,j)*qE) > EPS)*/ totalChargeF += topVec(j)*dA*qE;

			totalChargeF += (bottomVec(j) + topVec(j))*qE*dy*dt*(mu_n*E0 + sqrt(D_n/dt));


			///reduces the value of the total (global) charge by the amount leaving the computational domain

		}


		///-----RENORMALIZATION COMPENSATION FOR CHARGE LEAKAGE-------///

		if (tCount != 0){

			renormCharge = abs(totalChargeF/totalCharge);
			if (carrierType == true) u_n *= renormCharge;
			else u_p *= renormCharge;

		}

		totalCharge = 0.0;

		for (int i = 0; i < M; i++){
			for (int j = 0; j < J; j++){

				if (carrierType == true){

					totalCharge += -u_n(i,j)*dA*qE;
					chargeCenterX += -u_n(i,j)*dA*(i*dx)*qE;

				}

				else {

					totalCharge += u_p(i,j)*dA*qE;
					chargeCenterX += u_p(i,j)*dA*(i*dx)*qE;

				}


			}
		}


		chargeCenterX /= totalCharge;
		chargeCenterX += bottom;

		if (depthReach == true) {

			///create output files for displaying cluster as it passes measurement plane


		}

		printf("Renormalization factor: %.16f\n",renormCharge);
		cout << "Grid points over threshold above measurement depth: " << passDepthCheck << endl;
		cout << "Depth of Center of Charge: " << chargeCenterX << endl;
		cout << "Fractional Change in Total Charge: " << (totalCharge - totalChargeF)/totalChargeF  << endl;
		cout << "Total Charge: " << totalCharge << " C" << endl;
		cout << "Time Step = " << tCount << endl;
		cout << "t = " << dt*tCount << " s" << endl;
		cout << endl;

		if (interactionType == true){

			ChrgLk_Coulomb.push_back(totalCharge0-totalChargeF);
			CCX_Coulomb.push_back(chargeCenterX);
			totalChrg_Coulomb.push_back(totalChargeF);

		}

		else {

			ChrgLk_Free.push_back(totalCharge0-totalChargeF);
			CCX_Free.push_back(chargeCenterX);
			totalChrg_Free.push_back(totalChargeF);

		}

		threshCount = 0;
		totalCharge = 0.0;
		chargeCenterX = 0.0;

		ClearVec(bx_n);
		ClearVec(by_n);

		tCount++;

	} while (passDepthCheck != 0);


	///--------------------------END LOOP----------------------------///


	//tCountFinal = tCount;
	/*

	for (int j = 0; j < J; j++) xIntDim(j) = tB(j)-tA(j);

	int MX = maximum(xIntDim);
	int mX = maxElement(xIntDim);

	VectorXf dX(MX);
	VectorXf temp(MX); ///stores columns of linChrgDens matrices for integration
	VectorXf chargeCollected(J);

	*/

	time(&end);
	totalTime = difftime(end, start);

	if (interactionType == true){

		/*

		diffVec(CCX_Coulomb,CCXdelta_Coulomb);

		for (int i = 0; i < MX; i++) dX(i) = CCXdelta_Coulomb[tA(mX)+i];

		for (int j = 0; j < J; j++){
			for (int i = 0; i < MX; i++){

					if (linChrgDensX_Coulomb(i,j) != 0.0) temp(i) = linChrgDensX_Coulomb(i,j);
					else temp(i) = 0.0;
				}

				chargeCollected(j) = Integrate(temp,dX);
				ClearVec(temp);
		}

		*/

		//PrintChargeDensity(dx,dy,X_muGrid,Y_muGrid,u_n,carrierType,outFileFINAL_Coulomb);
		PrintNumberDensity(dx,dy,X_muGrid,Y_muGrid,u_n,carrierType,outFileFINAL_Coulomb);

		PrintVector(TOTcounter_n,dy,outFileTOT_Coulomb);
		//PrintVector(chargeCollected, dy, outFileChargeCollect_Coulomb);
		PrintVector(collectionTime, dy, outFileCollectionTime_Coulomb);
		PrintVector(CldSz_Coulomb,dt,outFileCloudSize_Coulomb);
		PrintVector(ChrgLk_Coulomb,dt,outFileChargeLeak_Coulomb);
		PrintVector(totalChrg_Coulomb,dt,outFileTotalCharge_Coulomb);
		PrintVector(CCX_Coulomb,dt,outFileCCX_Trajectory_Coulomb);
		PrintVector(DiffusionTerm,dt,outFileDiffusionTerm_Coulomb);
		PrintVector(AdvectionTerm,dt,outFileAdvectionTerm_Coulomb);
		PrintMatrix(linChrgDensX_Coulomb, outFileLinChrgDens_Coulomb);

		*outFileData_Coulomb << "Final Depth of Center of Charge: " << CCX_Coulomb[CCX_Coulomb.size()-1] << " m" << endl;
		//*outFileData_Coulomb << "Collection Time (Free): " << collectionTime << " s" << endl;
		*outFileData_Coulomb << "Total Compuational Time: " << totalTime << " s" << endl;



		//printf("Collection Time (Free): %.12f\n", collectionTime);
		printf("Total Computational Time: %.12f\n", totalTime);

		cout << "Final Depth of Center of Charge: " << CCX_Coulomb[CCX_Coulomb.size()-1] << " m" << endl;
		cout << "Threshold Density: " << threshDens << " C/m^2" << endl;
		cout << "Initial Total Charge: " << totalCharge0 << " C" << endl;
		cout << "Final Total Charge: " << totalChargeF << " C" << endl;


		CCX_Coulomb.clear();
		ChrgLk_Coulomb.clear();
		totalChrg_Coulomb.clear();
		ClearVec(TOTcounter_n);


		outFileDiffusionTerm_Coulomb->close();
		outFileAdvectionTerm_Coulomb->close();
		outFileFINAL_Coulomb->close();
		outFileCloudSize_Coulomb->close();
		outFileTOT_Coulomb->close();
		outFileCCX_Trajectory_Coulomb->close();
		outFileChargeLeak_Coulomb->close();
		outFileTotalCharge_Coulomb->close();
		outFileLinChrgDens_Coulomb->close();
		outFileDeltaXIntg_Coulomb->close();
		outFileData_Coulomb->close();


	}


	else {

		/*

		diffVec(CCX_Free,CCXdelta_Free);

		for (int i = 0; i < MX; i++) dX(i) = CCXdelta_Free[tA(mX)+i];

		for (int j = 0; j < J; j++){
			for (int i = 0; i < MX; i++){

				if (linChrgDensX_Free(i,j) != 0.0) temp(i) = linChrgDensX_Free(i,j);
							else temp(i) = 0.0;
			}

				chargeCollected(j) = Integrate(temp,dX);
				ClearVec(temp);
		}

		*/

		//PrintChargeDensity(dx,dy,X_muGrid,Y_muGrid,u_n,carrierType,outFileFINAL_Free);
		PrintNumberDensity(dx,dy,X_muGrid,Y_muGrid,u_n,carrierType,outFileFINAL_Free);

		PrintVector(TOTcounter_n,dy,outFileTOT_Free);
		//PrintVector(chargeCollected, dy, outFileChargeCollect_Free);
		PrintVector(collectionTime, dy, outFileCollectionTime_Free);
		PrintVector(CldSz_Free,dt,outFileCloudSize_Free);
		PrintVector(ChrgLk_Free,dt,outFileChargeLeak_Free);
		PrintVector(totalChrg_Free,dt,outFileTotalCharge_Free);
		PrintVector(CCX_Free,dt,outFileCCX_Trajectory_Free);
		PrintVector(DiffusionTerm,dt,outFileDiffusionTerm_Free);
		PrintVector(AdvectionTerm,dt,outFileAdvectionTerm_Free);
		PrintMatrix(linChrgDensX_Free,outFileLinChrgDens_Free);

		*outFileData_Free << "Final Total Charge: " << totalChargeF << endl;
		*outFileData_Free << "Final Depth of Center of Charge: " << CCX_Free[CCX_Free.size()-1] << " m" << endl;
		//*outFileData_Free << "Collection Time (Free): " << collectionTime << " s" << endl;
		*outFileData_Free << "Total Compuational Time: " << totalTime << " s" << endl;


		//printf("Collection Time (Free): %.12f\n", collectionTime);
		printf("Total Computational Time: %.12f\n", totalTime);

		cout << "Final Depth of Center of Charge: " << CCX_Free[CCX_Free.size()-1] << " m" << endl;
		cout << "Threshold Density: " << threshDens << " C/m^2" << endl;
		cout << "Initial Total Charge: " << totalCharge0 << " C" << endl;


		CCX_Free.clear();
		ChrgLk_Free.clear();
		totalChrg_Free.clear();
		ClearVec(TOTcounter_n);


		outFileDiffusionTerm_Free->close();
		outFileAdvectionTerm_Free->close();
		outFileFINAL_Free->close();
		outFileCloudSize_Free->close();
		outFileTOT_Free->close();
		outFileCCX_Trajectory_Free->close();
		outFileChargeLeak_Free->close();
		outFileTotalCharge_Free->close();
		outFileLinChrgDens_Free->close();
		outFileDeltaXIntg_Free->close();
		outFileData_Free->close();


	}

	ICDensitySimMat->close();
	ICDensityDataMat->close();
	InitialGrid->close();

	return 0;

}

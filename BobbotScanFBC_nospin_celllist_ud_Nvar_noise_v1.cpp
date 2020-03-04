// Author: Bahnisikha Dutta
// Modified from: Code by Shengkai Li
//Date created

//Libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <list>
#include <algorithm>
#include <chrono>

//Initializations
using namespace std;
double norm(double x, double y);
double cross(double x1, double y1, double x2, double y2);
vector<double> pairForce0(double x10, double y10, double phi10,
	double x20, double y20, double phi20,
	double R, double FM0, double RB, double dphi);
vector<double> loosePos(double x1, double y1, double x2, double y2, double phi10, double phi20, double dphi, double R);

int main(int argc, char *argv[]) {
	double FM = atof(argv[1]);//Magnetic Strength
	double v0 = atof(argv[2]);// Velocity of Bots
	int n = atoi(argv[3]); //iteration number also used as rng seed
	double T = atof(argv[4]); //Total time of simulation
	double dt0 = atof(argv[5]); //Sampling duration
	const int nBots = atoi(argv[6]);// Number of robots
	double A =atof(argv[7]);// 3*0.1;// Translational noise
	double AP =atof(argv[8]);// 3*0.001;//Rotational noise	
	double FA = 0;// 0.1*nBots*0.06; //2*0.5437;//0.06/exp(-1); // air flow force (equivalent to drive force, in Newtons) BD: needs to be exponential decay and not inform across a cross section and over time
	double air_k0=1;//log(nBots);
		//printf ("test1, air_k0=%f\n", air_k0);
		//printf ("test5, N2=%d\n", N2);

	string fn = "result_FM_" + string(argv[1]) + "_VEL_" + string(argv[2]) + "_N_" + string(argv[6]) + "_A_" + string(argv[7]) + "_AP_" + string(argv[8]) + "_Exp_" + string(argv[3]) + ".dat"; //Name of result variable 
	ofstream myfile(fn);//stream the file
	srand(n);//initiate rng
	const double pi = 3.14159265358979323846;

	// SI unit throughout the code: m, s, kg, N=kg*m*s^-s
	double dt = 1e-3; // dt = 1e-4  integration timestep
	int L = floor(T / dt);// Number of simulation steps
	int L0 = floor(0.1*L); //??
	int L1 = floor(0.95*L);//??
	int idS = floor(dt0 / dt);//Number of samples to save
	// double A = (1e-3) * pow(dt, -0.5);


	// Parameters

	double D = 1; //??
	double Dr = 1; //??
	double miu = 0.143; // bot-bot coefficient of friction
	double miuW = 0.143; // bot-wall coefficient of friction
	double R0 = 0.03; // radius of the robots
	double AT = 2 * R0; // air flow active size
	double RA = 1.1*R0;//R0 + AT; // air flow active size considering the radius of the bot
	double RC = 0.025; //radius of curvature of trajectory
	double phi0 = 0.25*pi; // center orientation of each magnet slot
	double n_temp=nBots;
	double B1 = sqrt(n_temp/400);//0.50; //Boundary dimension 1
	double B2 = sqrt(n_temp/400); //Boundary dimension 2
		//printf ("test2, B1=%f\n", B1);
		//printf ("test3, B2=%f\n", B2);
	double kb = 1000; //??
	double kbb = 1000; //??
	double kbbd = 10; //??
	double m0 = 0.060; //??
	double MI = 0.5*m0*R0*R0; //??
	double eta = 1.0; //??
	double etaP = 0.0003; //??
	double FM0 = (FM / 1000)*9.8; // 0.020*9.8
	double RB = 0.004;//??
	double FD = eta*v0;//Drive force?
	double TD = etaP*v0 / RC;//Torque?
	double dphi = 15 * pi / 180; // slot range of loose magnet: phi0+-dphi
	const double eps = 1e-15; // a small value

	// cell-list info
	int NC_temp=floor(2*B1*100/8); 
	const int NC1 = NC_temp; //12 Number of Cells in dimension 1
	const int NC2 = NC_temp; //12 round(10 * B2 / B1) Number of Cells in dimension 2
		//printf ("test4, NC1=%d\n", NC1);
		//printf ("test5, NC2=%d\n", NC2);
	double DC1 = 2 * B1 / NC1; // Number of cells in x direction
	double DC2 = 2 * B2 / NC2; // Number of cells in y direction
	struct map { vector<int> ind; } g[NC1][NC2]; //what you trying to do ???
	vector<int> di = { 1,1,0,-1 }; // neighbor cell additive coordinates x 1st and 2nd quadrant
	vector<int> dj = { 0,1,1,1 }; // neighbor cell additive coordinates y 1st and 2nd quadrant
	vector<int> id1, id2; //??
	int id10, id20; //??
	int idx, idy, i1, i2, j1, j2;// ??

	// Initialization
	// Initialize an initial configuration of the bots so that non of them are
	// not overlapping with each other.
	double dd; //??
	double xT[nBots], yT[nBots]; // Final assignment of positions for one time step
	double xQ, yQ; // temp variable
	double x[2][nBots], y[2][nBots], phi[2][nBots];// position information of robots for multiple timesteps
	double xt[2][nBots], yt[2][nBots], phit[2][nBots];// velocity information of robots for multiple timesteps
	double xtt[nBots], ytt[nBots], phitt[nBots];// ??
	double Fbx[nBots], Fby[nBots]; // ??
	bool f;// ??
	
    int c = 0; // data sampling counter
	double rx, ry, d;//??
	double nx, ny;//??
	double vr[2], ff[2];//??
	double FM1, FM2, T1, T2, Tf; //??
	double fbx, fby;//??
	double Fbb, Fx, Fy;// ??
	double Fbbx[nBots], Fbby[nBots], Tor[nBots]; // bob-bot interaction forces and torque ??
	double vj[2][nBots], vk[2][nBots]; // ??
	double vn; // ??
	double d1, d2, d3, d4; // ??
	vector<double> F; // ??
	int i,j,jj,kk,k,q;//iteration dummy variables
	
	//Initialization of robot ??
	for (j = 0; j < nBots; j++) { //j: robot indices
		xt[0][j] = eps;
		yt[0][j] = eps;
		phit[0][j] = eps;
	}
    // Initialize positions while checking for overlaps
	for (j = 0; j < nBots; j++) {
		f = false;
		while (f == false) {
			xQ = 2 * (((double)rand() / RAND_MAX) - 0.5)*(B1 - R0);
			yQ = 2 * (((double)rand() / RAND_MAX) - 0.5)*(B2 - R0);
			f = true;
			for (int k = 0; k < j - 1; k++) {
				dd = norm(xQ - xT[k], yQ - yT[k]);
				if (dd < 2 * R0) {
					f = false;
					break;
				}
			}
		}
		xT[j] = xQ; //final assignment
		yT[j] = yQ;
	}
	
	// set robot credentials for 0th timestep


	for (j = 0; j < nBots; j++) {
		x[0][j] = xT[j];
		y[0][j] = yT[j];
		phi[0][j] = 2 * pi*((double)rand() / RAND_MAX);
	}
	
	// Iteration of simulation steps
	for (i = 0; i < L; i++) { //i:Iteration step
		
		//const double mean = 0.0;//B
    		//const double stddev = 0.5;//B
    		//default_random_engine generator_air;//B
		


    		//normal_distribution<double> dist_air(mean, stddev);//B
    		


		//FA += 0.1*A*(((double)rand() / RAND_MAX) - 0.5);//B add noise to airflow amplitude over time
		//FA += A*dist_air(generator_air);//B add noise to airflow amplitude over time

		// Cell registration
		for (j = 0; j < NC1; j++) {
			for (k = 0; k < NC2; k++) {
				g[j][k].ind.clear();
			}
		}
		//Repopulating g
		for (j = 0; j < nBots; j++) {
			idx = ceil((B1 + x[0][j]) / DC1) - 1;
			idy = ceil((B2 + y[0][j]) / DC2) - 1;
			idx = max(0, min(NC1 - 1, idx));
			idy = max(0, min(NC2 - 1, idy));
			g[idx][idy].ind.push_back(j);
		}

		// Check the cell-list
		/*
		for (int j = 0; j < NC; j++) {
		for (int k = 0; k < NC; k++) {
		cout << "(" << j << "," << k << "): ";
		for (int l = 0; l < g[j][k].ind.size(); l++) {
		cout << g[j][k].ind[l] << " ";
		}
		cout << "\n";
		}
		}
		*/
		
		// Bot-bot interaction initialization
		for (j = 0; j < nBots; j++) {
			Fbbx[j] = 0;
			Fbby[j] = 0;
			Tor[j] = 0;
			vj[0][j] = 0;
			vj[1][j] = 0;
			vk[0][j] = 0;
			vk[1][j] = 0;
		}

		// evaluation of bot-bot forces and torque
		for (i1 = 0; i1 < NC1; i1++) { //i1: cell x index 
			for (j1 = 0; j1 < NC2; j1++) {// j1: cell y index
				id1 = g[i1][j1].ind; // current cell robot index list
				id10 = id1.size(); // total number of robots in current cell
				// Cross-cell interaction
				for (q = 0; q < 4; q++) {
					i2 = i1 + di[q]; //i2: prospective adjacent cells x
					j2 = j1 + dj[q]; //j2: prospective adjacent cells y
					if ((i2 <= NC1 - 1) && (j2 <= NC2 - 1) && (i2 >= 0) && (j2 >= 0)) {// checking boundary conditions
						id2 = g[i2][j2].ind; //adjacent cell robot index list
						id20 = id2.size(); // total number of robots in adjacent cell
						
						for (jj = 0; jj < id10; jj++) {
							for (kk = 0; kk < id20; kk++) {
								j = id1[jj]; // pick a robot in current cell
								k = id2[kk]; // pick a robot in adjacent cell
								// Bobbot-bobbot force and torque evaluation
								{
									rx = x[0][k] - x[0][j]; // x distance
									ry = y[0][k] - y[0][j]; // y distance
									d = norm(rx, ry); // distance between centers of robots

									// magnetic force evaluation
									
									
									F = pairForce0(x[0][k], y[0][k], phi[0][k], 
										x[0][j], y[0][j], phi[0][j],
										R0, FM0, RB, dphi); //paiForce0: magnet force calculations, x,y,phi[. k]:neighbor cell robot;
															//x,y,phi [.,j]: current robot; R0: radius of bot, RB:radius of bead; dphi: range of magnet slot
									FM1 = F[0];//temp variables for forces and torques
									FM2 = F[1];
									T1 = F[2];
									T2 = F[3];
									Fbbx[j] += FM1; // Newton's third law
									Fbby[j] += FM2; 
									Fbbx[k] -= FM1;
									Fbby[k] -= FM2;
									Tor[j] += T2;
									Tor[k] += T1; // F: Magnetic interaction Fbbx,Fbby,Tor: final force on robot j and likewise for k; this keeps cumulatively adding

									// normal force and friction evaluation dem??
									
									
									if (d < (2 * R0)) {
										// normal force
										nx = rx / d;
										ny = ry / d;
										vn = (xt[0][j] - xt[0][k])*nx + (yt[0][j] - yt[0][k])*ny; // its the normal velocity component of interaction between the bots
										Fbb = kbb*(2 * R0 - d) + kbbd*vn; // Final exclusion force between robots
										Fx = Fbb*nx; // x direction exclusion force
										Fy = Fbb*ny; // y direction exclusion force
										Fbbx[j] -= Fx; // Newton's third law  
										Fbby[j] -= Fy;
										Fbbx[k] += Fx;
										Fbby[k] += Fy; //Fx,Fy: is the exclusion force

										// friction
										vj[0][j] = xt[0][j] - phit[0][j] * R0*ny;
										vj[1][j] = yt[0][j] + phit[0][j] * R0*nx;
										vk[0][k] = xt[0][k] + phit[0][k] * R0*ny;
										vk[1][k] = yt[0][k] - phit[0][k] * R0*nx;
										vr[0] = vj[0][j] - vk[0][k];
										vr[1] = vj[1][j] - vk[1][k];
										ff[0] = -miu*Fbb*vr[0] / norm(vr[0], vr[1]);
										ff[1] = -miu*Fbb*vr[1] / norm(vr[0], vr[1]);
										Tf = cross(rx, ry, ff[0], ff[1]);
										Fbbx[j] += ff[0];
										Fbby[j] += ff[1];
										Fbbx[k] -= ff[0];
										Fbby[k] -= ff[1];
										Tor[j] += Tf;
										Tor[k] += Tf;
									}
								}
							}
						}
					}
				}

				// Same-cell interaction
				for (jj = 0; jj < id10; jj++) {
					for (kk = 0; kk < id10; kk++) {
						j = id1[jj];
						k = id1[kk];
						// Bobbot-bobbot force and torque evaluation
						if (j<k) {
							rx = x[0][k] - x[0][j];
							ry = y[0][k] - y[0][j];
							d = norm(rx, ry);

							// magnetic force evaluation
							F = pairForce0(x[0][k], y[0][k], phi[0][k],
								x[0][j], y[0][j], phi[0][j],
								R0, FM0, RB, dphi);
							FM1 = F[0];
							FM2 = F[1];
							T1 = F[2];
							T2 = F[3];
							Fbbx[j] += FM1;
							Fbby[j] += FM2;
							Fbbx[k] -= FM1;
							Fbby[k] -= FM2;
							Tor[j] += T2;
							Tor[k] += T1;

							// normal force and friction evaluation
							if (d < (2 * R0)) {
								// normal force
								nx = rx / d;
								ny = ry / d;
								vn = (xt[0][j] - xt[0][k])*nx + (yt[0][j] - yt[0][k])*ny;
								Fbb = kbb*(2 * R0 - d) + kbbd*vn;
								Fx = Fbb*nx;
								Fy = Fbb*ny;
								Fbbx[j] -= Fx;
								Fbby[j] -= Fy;
								Fbbx[k] += Fx;
								Fbby[k] += Fy;

								// friction
								vj[0][j] = xt[0][j] - phit[0][j] * R0*ny;
								vj[1][j] = yt[0][j] + phit[0][j] * R0*nx;
								vk[0][k] = xt[0][k] + phit[0][k] * R0*ny;
								vk[1][k] = yt[0][k] - phit[0][k] * R0*nx;
								vr[0] = vj[0][j] - vk[0][k];
								vr[1] = vj[1][j] - vk[1][k];
								ff[0] = -miu*Fbb*vr[0] / norm(vr[0], vr[1]);
								ff[1] = -miu*Fbb*vr[1] / norm(vr[0], vr[1]);
								Tf = cross(rx, ry, ff[0], ff[1]);
								Fbbx[j] += ff[0];
								Fbby[j] += ff[1];
								Fbbx[k] -= ff[0];
								Fbby[k] -= ff[1];
								Tor[j] += Tf;
								Tor[k] += Tf;
							}
						}
					}
				}
			}
		}
		
			// boundary exclusion force 
		for (j = 0; j < nBots; j++) {
			d1 = 0; //flags
			d2 = 0;
			d3 = 0;
			d4 = 0;
			if (x[0][j] < (-B1 + R0)) d1 = 1; //hitting left boundary
			if (x[0][j] > (B1 - R0)) d2 = 1; //hitting right boundary
			if (y[0][j] < (-B2 + R0)) d3 = 1; //hitting bottom boundary
			if (y[0][j] > (B2 - R0)) d4 = 1; //hitting up boundary
			Fbx[j] = -kb*((x[0][j] - (-B1 + R0))*d1 + (x[0][j] - (B1 - R0))*d2); //Normal force in x direction
			Fby[j] = -kb*((y[0][j] - (-B2 + R0))*d3 + (y[0][j] - (B2 - R0))*d4); //Normal force in y direction
			fby = -miuW*abs(Fbx[j])*copysign(1, yt[0][j] + phit[0][j] * R0*(d2 - d1)); //Frictional force in y direction
			fbx = -miuW*abs(Fby[j])*copysign(1, xt[0][j] + phit[0][j] * R0*(d3 - d4)); //Frictional force in x direction
			Fbx[j] += fbx;
			Fby[j] += fby;
			if (x[0][j] < (-B1 + R0)) Tor[j] += -fby*R0;
			if (x[0][j] > (B1 - R0)) Tor[j] += fby*R0;
			if (y[0][j] < (-B2 + R0)) Tor[j] += fbx*R0;
			if (y[0][j] > (B2 - R0)) Tor[j] += -fbx*R0;
		}

		// air flow force : constant force not right
		for (j = 0; j < nBots; j++) {
                        //RA+=A*(((double)rand() / RAND_MAX) - 0.5);//B
			d1 = 0; //flags
			d2 = 0;
			d3 = 0;
			d4 = 0;
			if (x[0][j] < (-B1 + RA)) d1 = 1; 
			if (x[0][j] > (B1 - RA)) d2 = 1;
			if (y[0][j] < (-B2 + RA)) d3 = 1;
			if (y[0][j] > (B2 - RA)) d4 = 1;

			//Fbx[j] += FA*(d1 - d2);
			//Fby[j] += FA*(d3 - d4); //Fbx,Fby:  Final boundary interaction forces
			//Fbx[j] += (FA*exp(((x[0][j] - (-B1 + RA))*d1 - (x[0][j] - (B1 - RA))*d2)))*(d1-d2);//B
			//Fby[j] += (FA*exp(((y[0][j] - (-B2 + RA))*d3 - (y[0][j] - (B2 - RA))*d4)))*(d3-d4); //Fbx,Fby:  Final boundary interaction forces//B

		    	//const double mean = 0.0;//B
    			//const double stddev = 0.5;//B
    			//default_random_engine generator4;//B
			//default_random_engine generator5;//B


    			//normal_distribution<double> dist4(mean, stddev);//B
    			//normal_distribution<double> dist5(mean, stddev);//B

			//Fbx[j] += (FA*exp(-(x[0][j] - (-B1)-R0)/RA)*d1 + FA*exp(-(B1-x[0][j]-R0)/RA)*d2)*(d1-d2);//B
			//Fby[j] += (FA*exp(-(y[0][j] - (-B2)-R0)/RA)*d3 + FA*exp(-(B2-y[0][j]-R0)/RA)*d4)*(d3-d4);//B //Fbx,Fby:  Final boundary interaction forces//B






			Fbx[j] += (FA*exp(-air_k0*((x[0][j] - (-B1)-R0)/RA)) - FA*exp(-air_k0*((B1-x[0][j]-R0)/RA)));//B
			Fby[j] += (FA*exp(-air_k0*((y[0][j] - (-B2)-R0)/RA)) - FA*exp(-air_k0*((B2-y[0][j]-R0)/RA)));//B //Fbx,Fby:  Final boundary interaction forces//B
			//Fbx[j] += A*(((double)rand() / RAND_MAX) - 0.5);//B
			//Fby[j] += A*(((double)rand() / RAND_MAX) - 0.5);//B
			//Fbx[j] += A*dist4(generator4);//B
			//Fby[j] += A*dist5(generator5);//B

			


		}




		// Sum up forces and torques
		for (j = 0; j < nBots; j++) {
			xtt[j] = FD * cos(phi[0][j] + phi0) - eta*xt[0][j] + Fbx[j] + Fbbx[j]; // drive+ floor friction+ boundary+bot-bot
			ytt[j] = FD * sin(phi[0][j] + phi0) - eta*yt[0][j] + Fby[j] + Fbby[j];

			if (j < 10) { // play with this for higher N
				phitt[j] = TD - etaP*phit[0][j] + Tor[j];
			}
			else {
				phitt[j] = -TD - etaP*phit[0][j] + Tor[j];
			}

    			//const double mean = 0.0;//B
    			//const double stddev = 0.01;//B
    			//default_random_engine generator1;//B
			//default_random_engine generator2;//B
			//default_random_engine generator3;//B


    			//normal_distribution<double> dist1(mean, stddev);//B
			//normal_distribution<double> dist2(mean, stddev);//B
			//normal_distribution<double> dist3(mean, stddev);//B

			//xtt[j] += A*dist1(generator1);//B
			//ytt[j] += A*dist2(generator2);//B
			//phitt[j] += AP*dist3(generator3);//B


			xtt[j] += A*(((double)rand() / RAND_MAX) - 0.5); //stochastic term; change to white gaussian
			ytt[j] += A*(((double)rand() / RAND_MAX) - 0.5); // stochastic term
			phitt[j] += AP*(((double)rand() / RAND_MAX) - 0.5); //stochastic term


			xtt[j] = xtt[j] / m0;
			ytt[j] = ytt[j] / m0;
			phitt[j] = phitt[j] / MI;
		}
		
		
		// Dynamics integration (forward Euler change to RK2)
		for (j = 0; j < nBots; j++) {
			xt[1][j] = xt[0][j] + dt*xtt[j];
			yt[1][j] = yt[0][j] + dt*ytt[j];
			phit[1][j] = phit[0][j] + dt*phitt[j];

			x[1][j] = x[0][j] + dt*xt[0][j];
			y[1][j] = y[0][j] + dt*yt[0][j];
			phi[1][j] = phi[0][j] + dt*phit[0][j];

			xt[0][j] = xt[1][j];
			yt[0][j] = yt[1][j];
			phit[0][j] = phit[1][j];

			x[0][j] = x[1][j];
			y[0][j] = y[1][j];
			phi[0][j] = phi[1][j];
		}


		//saving values in file
	
		if (i % idS == 0) { // (i > L1)
		       // overall position
			for (j = 0; j < nBots; j++) {
				myfile << phi[0][j] << " ";
			}
			myfile << "\n";
			for (j = 0; j < nBots; j++) {
				myfile << x[0][j] << " ";
			}
			myfile << "\n";
			for (j = 0; j < nBots; j++) {
				myfile << y[0][j] << " ";
			}
			myfile << "\n" << "\n";
		}
}


	myfile.close();
	return 0;
}

double norm(double x, double y)
{
	double n = sqrt(x*x + y*y);
	return n;
}

double cross(double x1, double y1, double x2, double y2) // ??
{
	double c = x1*y2 - x2*y1;
	return c;
}

vector<double> magF0(double x, double y, double FM0, double RB) {
	vector<double> F;
	F.resize(2, 0);
	double r = norm(x, y);
	double dH = 2 * RB; // The hardcore model
	double d0 = 0.0015;
	double F0;
	F0 = FM0*exp(-(r - dH) / d0);
	F[0] = F0*x / r;
	F[1] = F0*y / r;
	return F;
}

vector<double> pairForce0(double x10, double y10, double phi10,
	double x20, double y20, double phi20,
	double R, double FM0, double RB, double dphi) {
	const double pi = 3.14159265358979323846;
	vector<double> F;
	F.resize(4, 0);
	double R0 = R - RB;
	double x, y;
	double r1[2], r2[2];
	vector<double> F12;
	double phi1[4], phi2[4], x1[4], y1[4], x2[4], y2[4]; // magnet coordinates
	for (int j = 0; j < 4; j++) {
		phi1[j] = phi10 + j*0.5*pi;
		phi2[j] = phi20 + j*0.5*pi;
		x1[j] = x10 + R0*cos(phi1[j]);
		y1[j] = y10 + R0*sin(phi1[j]);
		x2[j] = x20 + R0*cos(phi2[j]);
		y2[j] = y20 + R0*sin(phi2[j]);
	}
	double T1 = 0;
	double T2 = 0;

	double d = 1e6;// ??
	double dN;
	int i0 = 0;
	int j0 = 0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dN = norm(x1[i] - x2[j], y1[i] - y2[j]);
			if (dN < d) {
				i0 = i;
				j0 = j;
				d = dN;
			}
		}
	}

	vector<double> phi = loosePos(x10, y10, x20, y20, phi1[i0], phi2[j0], dphi, R0);
	x1[i0] = x10 + R0*cos(phi[0]);
	y1[i0] = y10 + R0*sin(phi[0]);
	x2[j0] = x20 + R0*cos(phi[1]);
	y2[j0] = y20 + R0*sin(phi[1]);
	x = x1[i0] - x2[j0];
	y = y1[i0] - y2[j0];
	if ((x*x + y*y) < (400 * RB*RB)) { // ??
		F12 = magF0(x, y, FM0, RB);
		F[0] += F12[0];
		F[1] += F12[1];

		r2[0] = R0*cos(phi2[j0]);
		r2[1] = R0*sin(phi2[j0]);
		T2 += cross(r2[0], r2[1], F12[0], F12[1]);

		r1[0] = R0*cos(phi1[i0]);
		r1[1] = R0*sin(phi1[i0]);
		T1 += cross(r1[0], r1[1], -F12[0], -F12[1]);
	}
	F[2] = T1;
	F[3] = T2;
	return F;
}

vector<double> loosePos(double x1, double y1, double x2, double y2, double phi10, double phi20, double dphi, double R) {
	vector<double> phi;
	phi.resize(2, 0);
	double dphi1, dphi2;
	double rr[2];
	double u1[2], u2[2];
	rr[0] = x2 - x1;
	rr[1] = y2 - y1;
	double r = norm(rr[0], rr[1]);
	rr[0] = rr[0] / r;
	rr[1] = rr[1] / r;
	u1[0] = cos(phi10);
	u1[1] = sin(phi10);
	u2[0] = cos(phi20);
	u2[1] = sin(phi20);
	double s1 = cross(rr[0], rr[1], u1[0], u1[1]);
	double DP1 = acos(rr[0] * u1[0] + rr[1] * u1[1]);
	DP1 = DP1*copysign(1, s1);
	double s2 = cross(rr[0], rr[1], -u2[0], -u2[1]);
	double DP2 = -acos(-rr[0] * u2[0] - rr[1] * u2[1]);
	DP2 = DP2*copysign(1, s2);

	if (abs(DP1) < dphi) {
		dphi2 = DP2;
	}
	else {
		if (DP1 > 0) {
			dphi2 = DP2 - asin(R *sin(DP1 - dphi) / sqrt(R*R + r*r - 2 * R*r*cos(DP1 - dphi)));
		}
		else {
			dphi2 = DP2 - asin(R*sin(DP1 + dphi) / sqrt(R*R + r*r - 2 * R*r*cos(DP1 + dphi)));
		}
	}

	if (abs(DP2) < dphi) {
		dphi1 = DP1;
	}
	else {
		if (DP2 > 0) {
			dphi1 = DP1 - asin(R*sin(DP2 - dphi) / sqrt(R*R + r*r - 2 * R*r*cos(DP2 - dphi)));
		}
		else {
			dphi1 = DP1 - asin(R*sin(DP2 + dphi) / sqrt(R*R + r*r - 2 * R*r*cos(DP2 + dphi)));
		}
	}

	dphi1 = max(-dphi, dphi1);
	dphi1 = min(dphi, dphi1);
	dphi2 = max(-dphi, dphi2);
	dphi2 = min(dphi, dphi2);
	phi[0] = phi10 - dphi1;
	phi[1] = phi20 + dphi2;

	double D = pow(x1 + R*cos(phi[0]) - x2 - R*cos(phi[1]), 2)
		+ pow(y1 + R*sin(phi[0]) - y2 - R*sin(phi[1]), 2);
	double d0 = 0.015;
	if (D > (d0*d0)) {
		phi[0] = phi10;
		phi[1] = phi20;
	}
	return phi;
}

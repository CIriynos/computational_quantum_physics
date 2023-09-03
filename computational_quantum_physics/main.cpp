#include <iostream>
#include "Wavefunction.h"
#include "NumericalGrid.h"
#include "ExportDataToMatFile.h"
#include "TDSESolverFD.h"
#include "TDSESolverFFT.h"
#include "ImagTimePropagation.h"
#include "TDSESolverFDInSH.h"
#include <cmath>
using namespace std;


void test_1d()
{
	double delta_t = 0.002;
	double t_span = 2.0;
	double left_bound = -10;
	double right_bound = 10;
	int N = 1000;

	NumericalGrid1D grid(N, right_bound - left_bound, (right_bound + left_bound) / 2);
	auto wave1 = create1DWaveByExpression(grid, makeGaussPkg1D(1, 0, 10));
	auto potiential = create1DWaveByExpression(grid, [](double) { return 0; });
	TDSESolverFFT<1> solver1(potiential, wave1, delta_t, t_span, PERIODIC_BOUNDARY_COND);
	TDSESolverFD<1> solver2(potiential, wave1, delta_t, t_span, PERIODIC_BOUNDARY_COND);
	vector< WaveMsg<1> > buffer;
	buffer.reserve(10000);

	int state = 0; int k = 0;
	do {
		state = get<0>(solver1.execute(5));
		solver2.execute(5);
		buffer.push_back(WaveMsg<1>(solver1.getCurrentState(), "fft", 1));
		buffer.push_back(WaveMsg<1>(solver2.getCurrentState(), "fd", 1));
		k++;
	} while (state != TDSE_REACHED_END);

	WaveAnimation1DToMat(buffer, REAL_MODE);
}


void test_2d() 
{
	double delta_t = 0.01;
	double t_span = 3.0;
	double x_len = 20;
	double y_len = 20;
	int Nx = 200;
	int Ny = 200;

	NumericalGrid2D grid(Nx, x_len, 0.0, Ny, y_len, 0);
	auto wave1 = create2DWaveByExpression(grid, makeGaussPkg2D(1, 1, 0, 0, 10, 0));
	auto potiential = create2DWaveByExpression(grid, [](double, double) { return 0; });
	
	TDSESolverFFT<2> solver(potiential, wave1, delta_t, t_span, PERIODIC_BOUNDARY_COND);

	vector< WaveMsg<2> > buffer;
	buffer.reserve(10000);

	int state = 0;
	do {
		buffer.push_back(WaveMsg<2>(solver.getCurrentState(), "fft", 1));
		state = get<0>(solver.execute(5));
	} while (state != TDSE_REACHED_END);

	WaveAnimation2DToMat(buffer, NORM_MODE);
}

void test_3d()
{
	int N = 75;
	double L = 20;
	double delta_t = 0.01;
	double t_span = 1.0;

	NumericalGrid3D grid(N, L, 0, N, L, 0, N, L, 0);
	auto wave1 = create3DWaveByExpression(grid, makeGaussPkg3D(1, 1, 1, 0, 0, 0, 10, 0, 0));
	auto potiential = create3DWaveByExpression(grid, [](double, double, double) { return 0; });

	TDSESolverFFT<3> solver(potiential, wave1, delta_t, t_span, PERIODIC_BOUNDARY_COND);
	vector< WaveMsg<2> > buffer;
	
	int state = 0;
	do {
		auto slice_indicator = make_pair(Z_DEGREE, (int)(grid.getCount(Z_DEGREE) / 2));
		buffer.push_back(WaveMsg<2>(solver.getCurrentState().getSlice<2>({slice_indicator}), "fft", 1));
		state = get<0>(solver.execute(5));
	} while (state != TDSE_REACHED_END);

	WaveAnimation2DToMat(buffer, NORM_MODE);
}

void test_imag_time_propagation()
{
	double omega = 1;
	NumericalGrid1D grid(1000, 20, 0);
	auto potiential = create1DWaveByExpression(grid, [omega](double x) { return 0.5 * x * x * omega * omega - 5; });

	constexpr int order = 5;
	vector< Wavefunction1D > waves;
	vector< Wavefunction1D*> waves_p;
	vector< TDSESolverFD<1> > solvers;
	vector< TDSESolverFD<1>*> solvers_p;
	for (int i = 0; i < order; i++) {
		waves.push_back(create1DWaveByExpression(grid, makeGaussPkg1D(1.0 + (double)i, 0, 1.0 + (double)i)));
		solvers.push_back(createSolverForImagTime(potiential, waves.back(), 0.01));
	}
	for (int i = 0; i < order; i++) {
		waves_p.push_back(&waves[i]);
		solvers_p.push_back(&solvers[i]);
	}

	ImagTimePropagationSolver<1, order> solver(solvers_p);

	solver.execute();

	vector< WaveMsg<1> > buffer;
	for (int i = 0; i < order; i++) {
		auto eigen_state = solver.getResult(i);
		auto energy = solver.getEnergy(i);
		cout << "Residual value = " << solver.checkResidual(i) << endl;
		cout << "Energy = " << energy << endl;
		buffer.push_back(WaveMsg<1>(eigen_state, "it", 1));
	}
	WaveAnimation1DToMat(buffer, NORM_MODE);
}

void test_imag_time_propagation_3d()
{
	int Nr = 200; int Ntheta = 200; int Nphi = 200;
	double Lr = 10; int Nl = 3; double dt = 0.01;

	auto sh_grid = createSphericalGrid(Lr, Nr, Ntheta, Nphi);
	auto po_f = [](double r, double theta, double phi) { return -1.0 / r; };
	auto potiential = create3DWaveByExpression(sh_grid, po_f);

	constexpr int M = 4;
	ImagTimePropagationSolverSH<M> solver(potiential, dt, Nl);
	solver.execute();

	vector< WaveMsg<3> > buffer;
	auto slice_indicator = make_pair(THETA_DEGREE, (int)(sh_grid.getCount(THETA_DEGREE) / 2));
	for (int i = 0; i < M; i++) {
		//auto eigen_state = solver.getResult(i);
		double energy = solver.getEnergey(i);
		cout << "Residual value sum = " << solver.checkResidualSum(i) << endl;
		cout << "Energy = " << energy << endl;
		//cout << "Norm of eigen_state = " << eigen_state.norm() << std::endl;
		//eigen_state.normalize();
		//buffer.push_back(WaveMsg<3>(eigen_state, "it", 2));
	}

	//WaveAnimation3DToMat(buffer, NORM_MODE, SPHERE_COORDINATE);
}


void test_3d_fd()
{
	int Nr = 50; int Ntheta = 50; int Nphi = 50;
	double Lr = 5; int Nl = 4;
	double dt = 0.01; double t_span = 5;

	auto sh_grid = createSphericalGrid(Lr, Nr, Ntheta, Nphi);
	auto po_f = [](double r, double theta, double phi) { return -1.0 / (r + 1e-6); };
	auto wave_f = [](double r, double theta, double phi) { return r * std::exp(-r); };

	auto potiential = create3DWaveByExpression(sh_grid, po_f);
	auto wave = create3DWaveByExpression(sh_grid, wave_f);
	wave.normalize();
	
	TDSESolverSH solver(potiential, wave, dt, t_span, Nl, PERIODIC_BOUNDARY_COND);
	vector< WaveMsg<2> > buffer;
	auto slice_indicator = make_pair(THETA_DEGREE, sh_grid.getCount(THETA_DEGREE) / 2);
	int state = 0;

	do {
		state = std::get<0>(solver.execute(10));
		buffer.push_back(WaveMsg<2>(solver.getCurrentState().getSlice<2>({ slice_indicator }), "sh", 1));
	} while (state != TDSE_REACHED_END);

	WaveAnimation2DToMat(buffer, NORM_MODE, POLAR_COORDINATE);
}

int main(int argc, char *argv[])
{
	test_imag_time_propagation_3d();
	return 0;
}
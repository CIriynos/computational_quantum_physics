#include <iostream>
#include "Wavefunction.h"
#include "NumericalGrid.h"
#include "ExportDataToMatFile.h"
#include "TDSESolverFD.h"
#include "TDSESolverFFT.h"
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
	auto wave1 = create2DWaveByExpression(grid, makeGaussPkg2D(1, 1, 0, 0, 0, 10));
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
	int N = 100;
	double L = 20;
	double delta_t = 0.01;
	double t_span = 2.0;

	NumericalGrid3D grid(N, L, 0, N, L, 0, N, L, 0);
	auto wave1 = create3DWaveByExpression(grid, makeGaussPkg3D(1, 1, 1, 0, 0, 0, 0, 10, 0));
	auto potiential = create3DWaveByExpression(grid, [](double, double, double) { return 0; });

	TDSESolverFFT<3> solver(potiential, wave1, delta_t, t_span, PERIODIC_BOUNDARY_COND);
	vector< WaveMsg<3> > buffer;
	
	int state = 0;
	do {
		buffer.push_back(WaveMsg<3>(solver.getCurrentState(), "fft", 1));
		state = get<0>(solver.execute(5));
	} while (state != TDSE_REACHED_END);

	WaveAnimation3DToMat(buffer, NORM_MODE);
}

int main(int argc, char *argv[])
{
	test_3d();
	return 0;
}
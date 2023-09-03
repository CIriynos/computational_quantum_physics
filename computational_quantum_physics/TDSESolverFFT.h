#ifndef __TDSE_SOLVER_FFT_H__
#define __TDSE_SOLVER_FFT_H__

#include "TDSESolver.h"
#include "MathOperator.h"
#include "Tools.h"

template<unsigned N>
class TDSESolverFFT : public TDSESolver<N>
{
public:
	using TDSESolver<N>::total_steps;
	using TDSESolver<N>::crt_step;
	using TDSESolver<N>::crt_state;

	TDSESolverFFT(const Wavefunction<N>& po_func, const Wavefunction<N>& is, double dt, double ts, int cond)
	{}

	virtual TDSESolverResult<N> execute(int ex_times);

private:
	Wavefunction<N> phase_factor_1;
	Wavefunction<N> phase_factor_2;
};

template<>
inline TDSESolverFFT<1>::TDSESolverFFT(const Wavefunction<1>& po_func, const Wavefunction<1>& is, double dt, double ts, int cond)
	: TDSESolver<1>(po_func, is, dt, ts, cond), phase_factor_1(is.getGrid()), phase_factor_2(is.getGrid())
{
	using namespace std::literals;
	assert(po_func.getGrid() == is.getGrid());

	auto grid = is.getGrid();
	Eigen::VectorXcd V = po_func.getSamplesView();
	auto initial_state_fft = fft(initial_state, std::vector<double>(1, 0));
	auto fft_grid = initial_state_fft.getGrid();

	phase_factor_1 = create1DWaveByExpression(fft_grid, [dt](double k) { return std::exp(-0.25i * dt * std::pow(k, 2)); });
	phase_factor_2 = createWaveByExpressionWithIndex(grid, [dt, &V](int i) { return std::exp(-1i * dt * V(i)); });
}


template<>
inline TDSESolverFFT<2>::TDSESolverFFT(const Wavefunction<2>& po_func, const Wavefunction<2>& is, double dt, double ts, int cond)
	: TDSESolver<2>(po_func, is, dt, ts, cond), phase_factor_1(is.getGrid()), phase_factor_2(is.getGrid())
{
	using namespace std::literals;
	assert(po_func.getGrid() == is.getGrid());

	auto grid = is.getGrid();
	Eigen::VectorXcd V = po_func.getSamplesView();
	auto initial_state_fft = fft(initial_state, std::vector<double>(2, 0));
	auto fft_grid = initial_state_fft.getGrid();

	phase_factor_1 = create2DWaveByExpression(fft_grid, [dt](double kx, double ky) {
		return std::exp(-0.25i * dt * kx * kx) * std::exp(-0.25i * dt * ky * ky); 
	});
	phase_factor_2 = createWaveByExpressionWithIndex(grid, [dt, &V](int i) {
		return std::exp(-1i * dt * V(i)); 
	});
}


template<>
inline TDSESolverFFT<3>::TDSESolverFFT(const Wavefunction<3>& po_func, const Wavefunction<3>& is, double dt, double ts, int cond)
	: TDSESolver<3>(po_func, is, dt, ts, cond), phase_factor_1(is.getGrid()), phase_factor_2(is.getGrid())
{
	using namespace std::literals;
	assert(po_func.getGrid() == is.getGrid());

	auto grid = is.getGrid();
	Eigen::VectorXcd V = po_func.getSamplesView();
	auto initial_state_fft = fft(initial_state, std::vector<double>(3, 0));
	auto fft_grid = initial_state_fft.getGrid();

	phase_factor_1 = create3DWaveByExpression(fft_grid, [dt](double kx, double ky, double kz) {
		return std::exp(-0.25i * dt * kx * kx) * std::exp(-0.25i * dt * ky * ky) * std::exp(-0.25i * dt * kz * kz);
	});
	phase_factor_2 = createWaveByExpressionWithIndex(grid, [dt, &V](int i) {
		return std::exp(-1i * dt * V(i));
	});
}


template<unsigned N>
inline TDSESolverResult<N> TDSESolverFFT<N>::execute(int ex_times)
{
	assert(ex_times == -1 || (ex_times >= 0 && ex_times <= total_steps));

	int state = TDSE_ONGOING;
	if (total_steps == crt_step) {
		return TDSESolverResult<N>(TDSE_REACHED_END, crt_state, 0.0);
	}
	if (ex_times == -1 || ex_times > (total_steps - crt_step)) {
		ex_times = total_steps - crt_step;
		state = TDSE_END;
	}

	std::clock_t c_start = std::clock();
	Wavefunction<N> x_space_state(crt_state);
	Wavefunction<N> k_space_state = fft(x_space_state, std::vector<double>(N, 0));

	for (int i = 1; i <= ex_times; i++) {
		k_space_state = k_space_state * phase_factor_1;
		x_space_state = ifft(k_space_state, std::vector<double>(N, 0));
		x_space_state = x_space_state * phase_factor_2;
		k_space_state = fft(x_space_state, std::vector<double>(N, 0));
		k_space_state = k_space_state * phase_factor_1;
	}
	x_space_state = ifft(k_space_state, std::vector<double>(N, 0));
	crt_state = x_space_state;

	std::clock_t c_end = std::clock();
	double duration = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
	crt_step += ex_times;

	return TDSESolverResult<N>(state, crt_state, duration);
}

#endif //__TDSE_SOLVER_FFT_H__
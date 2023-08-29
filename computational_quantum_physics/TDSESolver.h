#ifndef __TDSE_SOLVER_H__
#define __TDSE_SOLVER_H__

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "Wavefunction.h"

typedef Eigen::SparseMatrix< std::complex<double> > SpMat;

template<unsigned N>
using TDSESolverResult = std::tuple<int, Wavefunction<N>&, double>;

constexpr int TDSE_ONGOING = 1;
constexpr int TDSE_REACHED_END = 2;
constexpr int TDSE_END = 3;

constexpr int REFLECTING_BOUNDARY_COND = 1;
constexpr int PERIODIC_BOUNDARY_COND = 2;

template<unsigned N>
class TDSESolver
{
protected:
	Wavefunction<N> potiential_func;
	Wavefunction<N> initial_state;
	Wavefunction<N> crt_state;
	double delta_t;
	double time_span;

	//for step-by-step calculation
	int total_steps;
	int crt_step;

	int boundary_cond;

public:
	TDSESolver(const Wavefunction<N>& po_func, const Wavefunction<N>& is, double dt, double ts, int cond)
		: potiential_func(po_func), initial_state(is),
		crt_state(is), delta_t(dt), time_span(ts), boundary_cond(cond),
		crt_step(0), total_steps(static_cast<int>(std::floor(ts / dt)))
	{
		assert(po_func.getGrid() == is.getGrid());
		assert(dt > 0.0 && ts > 0.0);
		assert(cond == REFLECTING_BOUNDARY_COND || cond == PERIODIC_BOUNDARY_COND);
	}

	virtual TDSESolverResult<N> execute(int) = 0;

	int getCurrentStep() { return crt_step; }

	Wavefunction<N>& getCurrentState() { return crt_state; }
};


#endif //__TDSE_SOLVER_H__!

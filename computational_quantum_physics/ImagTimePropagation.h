#ifndef __IMAG_TIME_PROPAGATION_H__
#define __IMAG_TIME_PROPAGATION_H__

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "Wavefunction.h"
#include "MathOperator.h"
#include "TDSESolverFD.h"
#include "TDSESolverFDInSH.h"
#include "Tools.h"

constexpr double DEFAULT_MAX_T_SPAN = 10000.0;
constexpr double MINIMAL_RESIDUAL_ERROR = 1e-5;

template<unsigned N, unsigned M>
class ImagTimePropagationSolver
{
public:
	ImagTimePropagationSolver(const std::vector<TDSESolverFD<N>*>& solvers_)
		: solvers(solvers_),
		H(solvers_.front()->getHamiltonian()),
		minimal_error_of_residual(MINIMAL_RESIDUAL_ERROR)
	{
		assert(solvers.size() == M);
		for (int i = 0; i < M; i++) {
			eigen_waves_p.push_back(&(solvers[i]->getCurrentState()));
		}
	}

	std::complex<double> checkResidual(int order) 
	{
		TDSESolverFD<N>& solver = *solvers[order];
		auto tmp = H * solver.getCurrentState();
		//std::cout << "energy distribution : " << tmp.getSamplesHandler() << std::endl;
		auto en = tmp.innerProduct(solver.getCurrentState());
		auto residual_state = H * solver.getCurrentState() - solver.getCurrentState() * en;
		return residual_state.norm();
	}

	void execute() {
		double maxerr = 0.0;
		do {
			for (int i = 0; i < M; i++) {
				(solvers[i])->execute(1);
				(solvers[i])->getCurrentState().normalize();
			}
			schmidtOrthon(eigen_waves_p);
			maxerr = 1e-10; //set max
			for (int i = 0; i < M; i++) {
				maxerr = std::max(maxerr, checkResidual(i).real());
			}
			//for (int i = 0; i < M; i++) {
			//	std::cout << checkResidual(i).real() << " ";
			//}
			//std::cout << maxerr << std::endl;
			//getchar();
		} while (maxerr > minimal_error_of_residual);
	}

	Wavefunction<N>& getResult(int order) { return *eigen_waves_p[order]; }

	double getEnergy(int order) {
		return (H * getResult(order)).innerProduct(getResult(order)).real();
	}

	MathOperatorMatrix<N>& getHamiltonian() { return H; }

private:
	std::vector< TDSESolverFD<N>* > solvers;
	MathOperatorMatrix<N> H;
	std::vector< Wavefunction<N>* > eigen_waves_p;
	double minimal_error_of_residual;
};




template<unsigned N>
class ImagTimePropagationSolver<N, 1>
{
public:

	ImagTimePropagationSolver(TDSESolverFD<N> &solver_)
		: solver(solver_),
		H(solver_.getHamiltonian()),
		result_base_eigenstate(),
		minimal_error_of_residual(MINIMAL_RESIDUAL_ERROR)
	{}

	std::complex<double> checkResidual();

	void execute();

	Wavefunction<N>& getResult() { return result_base_eigenstate; }

	double getEnergy() { return (solver.getHamiltonian() * getResult()).innerProduct(getResult()).real(); }

	MathOperatorMatrix<N>& getHamiltonian() { return H; }

private:
	TDSESolverFD<N>& solver;
	MathOperatorMatrix<N> H;
	Wavefunction<N> result_base_eigenstate;
	double minimal_error_of_residual;
};





template<unsigned N>
inline std::complex<double> ImagTimePropagationSolver<N, 1>::checkResidual()
{
	auto tmp = H * solver.getCurrentState();
	//std::cout << "energy distribution : " << tmp.getSamplesHandler() << std::endl;
	auto en = tmp.innerProduct(solver.getCurrentState());
	auto residual_state = H * solver.getCurrentState() - solver.getCurrentState() * en;
	return residual_state.norm();
}

template<unsigned N>
inline void ImagTimePropagationSolver<N, 1>::execute()
{
	double err = 0;
	do {
		solver.execute(1);
		solver.getCurrentState().normalize();
		err = checkResidual().real();
		//std::cout << err << std::endl;
	} while (err > MINIMAL_RESIDUAL_ERROR);

	result_base_eigenstate = solver.getCurrentState();
}



template<unsigned N>
inline TDSESolverFD<N> createSolverForImagTime(const Wavefunction<N>& po_func, const Wavefunction<N>& is, double dt = 0.001)
{
	return TDSESolverFD<N>(po_func, is, dt, DEFAULT_MAX_T_SPAN, IMAG_TIME_PROPAGATION_COND);
}

inline TDSESolverFDInSH createSolverForImagTimeSH(const Wavefunction<1>& po_func, const Wavefunction<1>& is, int l, double dt = 0.001)
{
	return TDSESolverFDInSH(po_func, is, dt, DEFAULT_MAX_T_SPAN, IMAG_TIME_PROPAGATION_COND, l);
}




template<unsigned M>
class ImagTimePropagationSolverSH
{
public:

	typedef std::vector< Wavefunction1D* > CouplingWavesRef;
	typedef std::vector< TDSESolverFD<1> > CouplingSolvers;

	ImagTimePropagationSolverSH(const Wavefunction3D& potiential_func_, double dt, int lmax)
		: potiential_func(potiential_func_.getSlice<1>({ std::make_pair(THETA_DEGREE, 1), std::make_pair(PHI_DEGREE, 1) })),
		sh_grid(potiential_func_.getGrid()),
		bundle_of_result_eigen_state(M, Wavefunction3D(potiential_func_.getGrid())),
		bundle_of_tdse_solvers(M),
		bundle_of_r_states_p(M),
		r_rev_addition(create1DWaveByExpression(potiential_func.getGrid(), [](double r) { return 1.0 / r; })),
		r_addition(create1DWaveByExpression(potiential_func.getGrid(), [](double r) { return r; })),
		delta_t(dt), lmax(lmax)
	{
		int k = 0;
		for (int i = 0; i < M; i++) {
			for (int l = 0; l < lmax; l++) {
				for (int m = -l; m <= l; m++) {
					bundle_of_tdse_solvers[i].push_back(createSolverForImagTimeSH(potiential_func, init(potiential_func.getGrid(), i, l), l));
				}
			}
			k = 0;
			for (int l = 0; l < lmax; l++) {
				for (int m = -l; m <= l; m++) {
					bundle_of_r_states_p[i].push_back(&(bundle_of_tdse_solvers[i][k].getCurrentState()));
					k++;
				}
			}
		}
	}

	void normalize(int order) {
		for (auto rsp : bundle_of_r_states_p[order]) {
			(*rsp) = ((*rsp)) * (1.0 / get_norm_sum(order));
		}
	}

	void check_orthon() {
		int r_size = bundle_of_r_states_p[0].size();
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				std::complex<double> tmp = 0;
				for (int k = 0; k < r_size; k++) {
					tmp += (*bundle_of_r_states_p[i][k]).innerProduct(*bundle_of_r_states_p[j][k]);
				}
				std::cout << i << ", " << j << " inner product : " << tmp << std::endl;
			}
		}
	}

	void check_distribution(int order) {
		int r_size = bundle_of_r_states_p[order].size();
		std::complex<double> tmp = 0;
		for (int k = 0; k < r_size; k++) {
			tmp = (*bundle_of_r_states_p[order][k]).norm();
			std::cout << order << " of norm : " << tmp << std::endl;
		}
		std::cout << std::endl;
	}

	void execute() {
		double maxerr = 0;
		for (int i = 0; i < M; i++) {
			normalize(i);
		}
		int n = 0;
		int start_index = 0;
		do {
			for (int i = start_index; i < M; i++) {
				for (auto& s : bundle_of_tdse_solvers[i]) {
					s.execute(1);
				}
				normalize(i);
			}
			schmidtOrthonForCouplingWavesInSH(bundle_of_r_states_p);
			for (int i = start_index; i < M; i++) {
				normalize(i);
			}
			
			if (0) {
				std::cout << "n = " << n << std::endl;
				for (int i = 0; i < M; i++) {
					std::cout << "energy of " << i << " = " << getEnergey(i) << std::endl;
					//std::cout << "distrution of " << i << " : " << std::endl;
					//std::cout << "total_norm : " << get_norm_sum(i) << std::endl;
					//check_distribution(i);
					std::cout << "err of " << i << " = " << checkResidualSum(i).real() << std::endl;
				}
				std::cout << std::endl;
				getchar();
				//check_orthon();
			}

			maxerr = 0;
			for (int i = start_index; i < M; i++) {
				double tmp = checkResidualSum(i).real();
				maxerr = std::max(maxerr, tmp);
				if (i == start_index && tmp < 3e-3) {
					start_index += 1;
				}
				//std::cout << "err of " << i << " = " << checkResidualSum(i).real() << " ";
			}
			//std::cout << std::endl << "maxerr: " << maxerr << std::endl;
			//getchar();
			//check_orthon();
			n++;
		} while (maxerr > 3e-3 && n <= 100000);
	}

	std::complex<double> get_norm_sum(int order) {
		std::complex<double> result = 0;
		for (auto rsp : bundle_of_r_states_p[order]) {
			result += (*rsp).norm();
		}
		return std::sqrt(result);
	}

	std::complex<double> checkResidualSum(int order) {
		std::complex<double> res = 0.0;
		for (auto& s : bundle_of_tdse_solvers[order]) {
			auto crt_s = s.getCurrentState();
			auto tmp = s.getHamiltonian() * crt_s;
			auto en = tmp.innerProduct(crt_s);
			auto residual_state = s.getHamiltonian() * crt_s - crt_s * en;
			res += residual_state.norm();
		}
		return res;
	}

	Wavefunction3D& getResult(int order) {
		bundle_of_result_eigen_state[order] = convergeSHBasesToWave(sh_grid, bundle_of_r_states_p[order], lmax);
		return bundle_of_result_eigen_state[order];
	}

	double getEnergey(int order) {
		double result = 0.0;
		for (auto& s : bundle_of_tdse_solvers[order]) {
			auto crt_s = s.getCurrentState();
			auto tmp = (s.getHamiltonian() * crt_s).innerProduct(crt_s).real();
			result += tmp;
		}
		return result;
	}

private:

	Wavefunction1D init(const NumericalGrid1D& grid, int order, int l) {
		double t = (1.0 + sqrt((double)order + 1.0));
		//std::cout << order << " " << l << " " << t << std::endl;
		return create1DWaveByExpression(grid, [t, l](double r) { return r * std::exp(-r) * 0.1 / (std::pow((double)l + 1.0, 2) / t); });
	}

	Wavefunction1D potiential_func;
	NumericalGrid3D sh_grid;

	std::vector< CouplingSolvers > bundle_of_tdse_solvers;
	std::vector< CouplingWavesRef > bundle_of_r_states_p;

	Wavefunction1D r_rev_addition;
	Wavefunction1D r_addition;

	std::vector<Wavefunction3D> bundle_of_result_eigen_state;
	double delta_t;
	int lmax;
};























template<>
class ImagTimePropagationSolverSH<1>
{
public:
	ImagTimePropagationSolverSH(const Wavefunction3D& potiential_func_, double dt, int lmax)
		: potiential_func(potiential_func_.getSlice<1>({ std::make_pair(THETA_DEGREE, 1), std::make_pair(PHI_DEGREE, 1) })),
		sh_grid(potiential_func_.getGrid()),
		result_eigen_state(potiential_func_.getGrid()),
		r_rev_addition(create1DWaveByExpression(potiential_func.getGrid(), [](double r) { return 1.0 / r; })),
		r_addition(create1DWaveByExpression(potiential_func.getGrid(), [](double r) { return r; })),
		delta_t(dt), lmax(lmax)
	{
		tdse_solvers.reserve(lmax * lmax);
		for (int l = 0; l < lmax; l++) {
			for (int m = -l; m <= l; m++) {
				tdse_solvers.push_back(createSolverForImagTimeSH(potiential_func, init(potiential_func.getGrid(), l), l));
				r_states_p.push_back(&tdse_solvers.back().getCurrentState());
			}
		}
		//r_rev_addition.normalize();
	}
	
	void normalize() {
		for (auto rsp : r_states_p) {
			(*rsp) = ((*rsp)) * (1.0 / get_norm_sum());
		}
	}

	void execute(){
		double err = 0;
		normalize();
		do {
			for (auto& s : tdse_solvers) {
				s.execute(1);
			}
			normalize();
			err = checkResidualSum().real();
			std::cout << err << std::endl;
		} while (err > 1e-5);
	}

	std::complex<double> get_norm_sum() {
		std::complex<double> result = 0;
		for (auto rsp : r_states_p) {
			result += (*rsp).norm();
		}
		return std::sqrt(result);
	}

	std::complex<double> checkResidualSum() {
		std::complex<double> res = 0.0;
		for (auto& s : tdse_solvers) {
			auto crt_s = s.getCurrentState();
			auto tmp = s.getHamiltonian() * crt_s;
			auto en = tmp.innerProduct(crt_s);
			auto residual_state = s.getHamiltonian() * crt_s - crt_s * en;
			res += residual_state.norm();
		}
		return res;
	}
 
	Wavefunction3D& getResult() { 
		result_eigen_state = convergeSHBasesToWave(sh_grid, r_states_p, lmax);
		return result_eigen_state;
	}

	double getEnergey() {
		double result = 0.0;
		for (auto& s : tdse_solvers) {
			auto crt_s = s.getCurrentState() * r_rev_addition;
			result += (s.getHamiltonian() * s.getCurrentState()).innerProduct(s.getCurrentState()).real();
		}
		return result;
	}

private:

	Wavefunction1D init(const NumericalGrid1D& grid, int l) {
		return create1DWaveByExpression(grid, [](double r) { return r * std::exp(-r); });
	}

	Wavefunction1D potiential_func;
	NumericalGrid3D sh_grid;
	std::vector< TDSESolverFDInSH > tdse_solvers;
	std::vector< Wavefunction1D* > r_states_p;
	Wavefunction1D r_rev_addition;
	Wavefunction1D r_addition;

	Wavefunction3D result_eigen_state;
	double delta_t;
	int lmax;
};


#endif
#ifndef __TDSE_SOLVER_FD_H__
#define __TDSE_SOLVER_FD_H__

#include "TDSESolver.h"
#include "MathOperator.h"
#include "Tools.h"

template<unsigned N>
class TDSESolverFD : public TDSESolver<N>
{
public:
	using TDSESolver<N>::total_steps;
	using TDSESolver<N>::crt_step;
	using TDSESolver<N>::crt_state;
	using TDSESolver<N>::delta_t;
	using TDSESolver<N>::boundary_cond;
	using TDSESolver<N>::initial_state;
	using TDSESolver<N>::potiential_func;
	
	TDSESolverFD(const Wavefunction<N>& po_func, const Wavefunction<N>& is, double dt, double ts, int cond);

	virtual TDSESolverResult<N> execute(int ex_times);

	virtual void update_hamiltonian();

	const MathOperatorMatrix<N>& getHamiltonian() const { return infact_hamiltonian; }

protected:
	MathOperatorMatrix<N> infact_hamiltonian;
	MathOperatorMatrix<N> A_positive;
	MathOperatorMatrix<N> A_negative;
};




inline MathOperatorMatrix<1> createHamiltonianFD1D(const NumericalGrid<1>& grid, int cond)
{
	typedef Eigen::Triplet< std::complex<double> > T;
	assert(cond == REFLECTING_BOUNDARY_COND 
		|| cond == PERIODIC_BOUNDARY_COND
		|| cond == IMAG_TIME_PROPAGATION_COND);

	int cnt = grid.getTotalCount();
	double delta_x = grid.getDelta(0);
	std::vector<T> tripletList;
	tripletList.reserve(cnt * 4);
	double scaler = -1.0 / (2.0 * delta_x * delta_x);
	SpMat matrix(cnt, cnt);

	for (int i = 0; i < cnt; i++) {
		if (i - 1 >= 0) {
			tripletList.push_back(T(i, i - 1, scaler));
		}
		tripletList.push_back(T(i, i, -2.0 * scaler));
		if (i + 1 < cnt) {
			tripletList.push_back(T(i, i + 1, scaler));
		}
	}
	if (cond == PERIODIC_BOUNDARY_COND) {
		tripletList.push_back(T(0, cnt - 1, scaler * 1));
		tripletList.push_back(T(cnt - 1, 0, scaler * 1));
	}
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());
	return MathOperatorMatrix<1>(grid, matrix);
}


inline MathOperatorMatrix<1> createMatrixV1D(const NumericalGrid<1>& grid, const Wavefunction<1>& potiential_func)
{
	typedef Eigen::Triplet< std::complex<double> > T;
	int cnt = grid.getCount(0);
	SpMat mat(cnt, cnt);
	std::vector<T> tripletList;

	for (int i = 0; i < cnt; i++) {
		tripletList.push_back(T(i, i, potiential_func.getValueByIndex(i)));
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	MathOperatorMatrix1D op(grid, mat);
	return op;
}

inline void eliminationProcess(SpMat& A, Eigen::VectorXcd& X, Eigen::VectorXcd& B, int cnt)
{
	assert(A.cols() == X.rows() && A.cols() == B.rows());

	for (int m = 1; m < cnt; m++) {
		for (int l = m; l >= m - 1; l--) {
			A.coeffRef(m, l) -= (A.coeffRef(m, m - 1) / A.coeffRef(m - 1, m - 1)) * A.coeffRef(m - 1, l);
			B(m) -= A.coeffRef(m, m - 1) / A.coeffRef(m - 1, m - 1) * B(m - 1);
		}
	}
	X(cnt - 1) = B(cnt - 1) / A.coeffRef(cnt - 1, cnt - 1);
	for (int m = cnt - 2; m >= 0; m--) {
		X(m) = (B(m) - A.coeffRef(m, m + 1) * X(m + 1)) / A.coeffRef(m, m);
	}
}

inline void solveLinearEquationsProblemInFD1D(const MathOperatorMatrix<1>& matrix, Wavefunction<1>& X, const Wavefunction<1>& B, int cond)
{
	assert(matrix.getGrid() == B.getGrid() && X.getGrid() == B.getGrid());

	int cnt = matrix.getGrid().getCount();

	if (cond == REFLECTING_BOUNDARY_COND 
		|| cond == IMAG_TIME_PROPAGATION_COND)
	{
		SpMat mat(matrix.getMatrix());
		Eigen::VectorXcd bdata(B.getSamplesView());
		Eigen::VectorXcd xdata(cnt);
		eliminationProcess(mat, X.getSamplesHandler(), bdata, cnt);
	}
	else if (cond == PERIODIC_BOUNDARY_COND)
	{
		SpMat mat1(matrix.getMatrix());
		auto corner_num = mat1.coeffRef(0, cnt - 1);

		Eigen::VectorXcd udata = Eigen::VectorXcd::Zero(cnt);
		Eigen::VectorXcd vdata = Eigen::VectorXcd::Zero(cnt);
		udata(0) = 1;
		udata(cnt - 1) = 1;
		vdata(0) = corner_num;
		vdata(cnt - 1) = corner_num;

		mat1.coeffRef(0, 0) -= corner_num;
		mat1.coeffRef(0, cnt - 1) -= corner_num;
		mat1.coeffRef(cnt - 1, 0) -= corner_num;
		mat1.coeffRef(cnt - 1, cnt - 1) -= corner_num;

		SpMat mat2(mat1);

		Eigen::VectorXcd bdata1(B.getSamplesView());
		Eigen::VectorXcd xdata1 = Eigen::VectorXcd::Zero(cnt);
		Eigen::VectorXcd bdata2(udata);
		Eigen::VectorXcd xdata2 = Eigen::VectorXcd::Zero(cnt);

		eliminationProcess(mat1, xdata1, bdata1, cnt);
		eliminationProcess(mat2, xdata2, bdata2, cnt);

		std::complex<double> tmp1 = vdata.transpose() * xdata2;
		std::complex<double> tmp2 = vdata.transpose() * xdata1;
		auto next_func = xdata1 - xdata2 * (tmp2 / (1.0 + tmp1));

		X.getSamplesHandler() = next_func;
	}
}

template<>
inline TDSESolverFD<1>::TDSESolverFD(const Wavefunction<1>& po_func, const Wavefunction<1>& is, double dt, double ts, int cond)
	: TDSESolver<1>(po_func, is, dt, ts, cond),
	A_positive(is.getGrid()), A_negative(is.getGrid())
{
	assert(is.getGrid() == po_func.getGrid());
	this->update_hamiltonian();
}

template<>
inline void TDSESolverFD<1>::update_hamiltonian()
{
	using namespace std::literals;

	auto grid = initial_state.getGrid();
	auto delta_x = grid.getDelta(0);
	auto hamiltonian = createHamiltonianFD1D(grid, boundary_cond);
	auto I = createIdentityOperator(grid);
	auto V = createMatrixV1D(grid, potiential_func);
	infact_hamiltonian = hamiltonian + V;

	if (boundary_cond == IMAG_TIME_PROPAGATION_COND) {
		A_positive = I - infact_hamiltonian * (0.5 * delta_t);
		A_negative = I + infact_hamiltonian * (0.5 * delta_t);
	}
	else {
		auto D = hamiltonian * (-2.0);
		auto M = I + D * (delta_x * delta_x / 12);
		A_positive = M - (D * (-0.5) + M * V) * (0.5i * delta_t);
		A_negative = M + (D * (-0.5) + M * V) * (0.5i * delta_t);
	}
}


template<>
inline TDSESolverResult<1> TDSESolverFD<1>::execute(int ex_times)
{
	assert(ex_times == -1 || (ex_times >= 0 && ex_times <= total_steps));

	int state = TDSE_ONGOING;
	if (total_steps == crt_step) {
		return TDSESolverResult<1>(TDSE_REACHED_END, crt_state, 0.0);
	}
	if (ex_times == -1 || ex_times > (total_steps - crt_step)) {
		ex_times = total_steps - crt_step;
		state = TDSE_END;
	}

	std::clock_t c_start = std::clock();

	Wavefunction1D half_state(crt_state.getGrid());
	for (int i = 1; i <= ex_times; i++) {
		half_state = A_positive * crt_state;
		solveLinearEquationsProblemInFD1D(A_negative, crt_state, half_state, boundary_cond);
	}

	std::clock_t c_end = std::clock();
	double duration = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
	crt_step += ex_times;

	return TDSESolverResult<1>(state, crt_state, duration);
}

#endif // !__TDSE_SOLVER_FD_H__

#ifndef __TOOLS_OF_QUANTUM_H__
#define __TOOLS_OF_QUANTUM_H__

#include "Wavefunction.h"
#include "Util.h"
#include <cmath>

inline Wavefunction1D create1DWaveByExpression(const NumericalGrid1D& grid, std::function<std::complex<double>(double)> func)
{
	std::vector< std::complex<double> > buffer(grid.getTotalCount(), 0);

	double x = 0;
	for (int i = 0; i < grid.getTotalCount(); i++) {
		x = grid.index(i).x();
		buffer[i] = func(x);
	}
	Wavefunction1D wave(grid, buffer);
	return wave;
}

inline Wavefunction2D create2DWaveByExpression(const NumericalGrid2D& grid, std::function<std::complex<double>(double, double)> func)
{
	std::vector< std::complex<double> > buffer(grid.getTotalCount(), 0);

	double x = 0, y = 0;
	for (int i = 0; i < grid.getTotalCount(); i++) {
		x = grid.index(i).x();
		y = grid.index(i).y();
		buffer[i] = func(x, y);
	}
	Wavefunction2D wave(grid, buffer);
	return wave;
}

inline Wavefunction3D create3DWaveByExpression(const NumericalGrid3D& grid, std::function<std::complex<double>(double, double, double)> func)
{
	std::vector< std::complex<double> > buffer(grid.getTotalCount(), 0);

	double x = 0, y = 0, z = 0;
	for (int i = 0; i < grid.getTotalCount(); i++) {
		x = grid.index(i).x();
		y = grid.index(i).y();
		z = grid.index(i).z();
		buffer[i] = func(x, y, z);
	}
	Wavefunction3D wave(grid, buffer);
	return wave;
}

template<unsigned N>
inline Wavefunction<N> createWaveByExpressionWithIndex(const NumericalGrid<N>& grid, std::function<std::complex<double>(int)> func)
{
	std::vector< std::complex<double> > buffer(grid.getTotalCount(), 0);
	for (int i = 0; i < grid.getTotalCount(); i++) {
		buffer[i] = func(i);
	}
	Wavefunction<N> wave(grid, buffer);
	return wave;
}

inline std::function<std::complex<double>(double)> makeGaussPkg1D(double omega_x, double x0, double p0)
{
	using namespace std::literals;
	std::complex<double> C = 1.0 / std::pow(2 * PI * omega_x * omega_x, 0.25);
	return [C, omega_x, x0, p0](double x) {
		return C * std::exp(-std::pow((x - x0) / (2 * omega_x), 2)) * std::exp(1i * p0 * x);
	};
}

inline std::function<std::complex<double>(double, double)> makeGaussPkg2D(double omega_x, double omega_y, double x0, double y0, double px, double py)
{
	using namespace std::literals;
	auto Cx = 1.0 / std::pow(2 * PI * omega_x * omega_x, 0.25);
	auto Cy = 1.0 / std::pow(2 * PI * omega_y * omega_y, 0.25);

	return [Cx, Cy, omega_x, omega_y, x0, y0, px, py](double x, double y) {
		auto xpart = Cx * std::exp(-std::pow((x - x0) / (2 * omega_x), 2)) * std::exp(1i * px * x);
		auto ypart = Cy * std::exp(-std::pow((y - y0) / (2 * omega_y), 2)) * std::exp(1i * py * y);
		return xpart * ypart;
	};
}

inline std::function<std::complex<double>(double, double, double)> makeGaussPkg3D(double omega_x, double omega_y, double omega_z, double x0, double y0, double z0, double px, double py, double pz)
{
	using namespace std::literals;
	auto Cx = 1.0 / std::pow(2 * PI * omega_x * omega_x, 0.25);
	auto Cy = 1.0 / std::pow(2 * PI * omega_y * omega_y, 0.25);
	auto Cz = 1.0 / std::pow(2 * PI * omega_z * omega_z, 0.25);

	return [Cx, Cy, Cz, omega_x, omega_y, omega_z, x0, y0, z0, px, py, pz](double x, double y, double z) {
		auto xpart = Cx * std::exp(-std::pow((x - x0) / (2 * omega_x), 2)) * std::exp(1i * px * x);
		auto ypart = Cy * std::exp(-std::pow((y - y0) / (2 * omega_y), 2)) * std::exp(1i * py * y);
		auto zpart = Cz * std::exp(-std::pow((z - z0) / (2 * omega_z), 2)) * std::exp(1i * pz * z);
		return xpart * ypart * zpart;
	};
}

typedef NumericalGrid<3> SphericalGrid;

inline SphericalGrid createSphericalGrid(double rlength, int rN, int thetaN, int phiN)
{
	double offset_in_singular_point = (double)rlength / (double)rN;
	return NumericalGrid<3>(rN, rlength, rlength / 2 + offset_in_singular_point, thetaN, PI, PI / 2, phiN, 2 * PI, PI);
}

inline Wavefunction3D create3DWaveByCartesianExpr(const NumericalGrid3D& grid, std::function<std::complex<double>(double, double, double)> func)
{
	std::vector< std::complex<double> > buffer(grid.getTotalCount(), 0);

	for (int i = 0; i < grid.getTotalCount(); i++) {
		auto r = grid.index(i).x();
		auto theta = grid.index(i).y();
		auto phi = grid.index(i).z();
		auto x = r * std::sin(theta) * std::cos(phi);
		auto y = r * std::sin(theta) * std::sin(phi);
		auto z = r * std::cos(theta);
		buffer[i] = func(x, y, z);
	}
	Wavefunction3D wave(grid, buffer);
	return wave;
}

inline Wavefunction3D makeSHWave(const NumericalGrid3D& grid, int l, int m)
{
	using namespace std::literals;
	std::vector< std::complex<double> > buffer(grid.getTotalCount(), 0);

	auto C = std::pow(-1, (double)m) * sqrt(((2 * (double)l + 1) / (4 * PI)) * (std::tgamma(l - m + 1) / std::tgamma(l + m + 1)));
	//r theta phi
	for (int i = 0; i < grid.getTotalCount(); i++) {
		auto point = grid.index(i);
		auto r = point.x();
		auto theta = point.y();
		auto phi = point.z();
		buffer[i] = C * (double)std::assoc_legendrel(l, m, std::cos(theta)) * std::exp(1i * (double)m * phi);
	}
	Wavefunction3D wave(grid, buffer);
	//wave.normalize();
	return wave;
}

inline std::vector< Wavefunction1D > expandInSHBases(const Wavefunction3D& wave, int maxl)
{
	auto grid = wave.getGrid();
	auto r_grid = grid.subset<1>({ R_DEGREE });
	auto phi_theta_grid = grid.subset<2>({ THETA_DEGREE, PHI_DEGREE });

	std::vector< Wavefunction1D > waves(maxl * maxl, Wavefunction1D(r_grid));

	int k = 0;
	for (int l = 0; l < maxl; l++) {
		for (int m = -l; m <= l; m++) {
			auto sh_base = makeSHWave(grid, l, m).getSlice<2>({std::make_pair(R_DEGREE, 1)});
			sh_base.normalize();
			for (int r_index = 0; r_index < grid.getCount(0); r_index++) {
				auto r = r_index * grid.getDelta(0);
				auto phi_r = wave.getSlice<2>({std::make_pair(R_DEGREE, r_index)});
				waves[k].getSamplesHandler()(r_index) = phi_r.innerProduct(sh_base);
			}
			k++;
		}
	}

	return waves;
}

inline Wavefunction3D convergeSHBasesToWave(const NumericalGrid<3>& sh_grid, const std::vector< Wavefunction1D* >& bases, int maxNl)
{
	assert(bases.size() == (maxNl * maxNl));
	Wavefunction3D result(sh_grid);
	int Nr = sh_grid.getCount(R_DEGREE);
	double Lr = sh_grid.getLength(R_DEGREE);
	double Offr = sh_grid.getOffset(R_DEGREE);
	auto r_rev_f = [](double r, double theta, double phi) { return 1.0 / r; };
	auto r_rev_wave = create3DWaveByExpression(sh_grid, r_rev_f);

	int k = 0;
	for (int l = 0; l < maxNl; l++) {
		for (int m = -l; m <= l; m++) {
			auto sh_base = makeSHWave(sh_grid, l, m).getSlice<2>({ std::make_pair(R_DEGREE, 1) });
			sh_base.normalize();
			auto ex_sh_base = sh_base.upliftOnce(R_DEGREE, Nr, Lr, Offr);
			result = result + ex_sh_base * bases[k]->uplift(sh_grid) * r_rev_wave;
			//cout << "l = " << l << " , m = " << m << " : " << bases[k].norm() << endl;
			k++;
		}
	}
	return result;
}

inline Wavefunction3D convergeSHBasesToWave(const NumericalGrid<3>& sh_grid, const std::vector< Wavefunction1D >& bases, int maxNl)
{
	assert(bases.size() == (maxNl * maxNl));
	Wavefunction3D result(sh_grid);
	int Nr = sh_grid.getCount(R_DEGREE);
	double Lr = sh_grid.getLength(R_DEGREE);
	double Offr = sh_grid.getOffset(R_DEGREE);
	auto r_rev_f = [](double r, double theta, double phi) { return 1.0 / r; };
	auto r_rev_wave = create3DWaveByExpression(sh_grid, r_rev_f);

	int k = 0;
	for (int l = 0; l < maxNl; l++) {
		for (int m = -l; m <= l; m++) {
			auto sh_base = makeSHWave(sh_grid, l, m).getSlice<2>({ std::make_pair(R_DEGREE, 1) });
			sh_base.normalize();
			auto ex_sh_base = sh_base.upliftOnce(R_DEGREE, Nr, Lr, Offr);
			result = result + ex_sh_base * bases[k].uplift(sh_grid) * r_rev_wave;
			//cout << "l = " << l << " , m = " << m << " : " << bases[k].norm() << endl;
			k++;
		}
	}
	return result;
}

template<unsigned N>
inline std::vector< Wavefunction<N> > schmidtOrthon(const std::vector< Wavefunction<N> >& waves)
{
	std::vector< Wavefunction<N> > new_waves;
	new_waves.reserve(waves.size());

	for (int i = 0; i < waves.size(); i++) {
		auto tmp = waves[i];
		for (int j = 0; j < i; j++) {
			auto scaler = waves[i].inner_product(new_waves[j]) / new_waves[j].inner_product(new_waves[j]);
			tmp -= new_waves[j] * scaler;
		}
		new_waves[i] = tmp;
	}
	return new_waves;
}


template<unsigned N>
inline void schmidtOrthon(const std::vector< Wavefunction<N>* >& waves)
{
	for (int i = 0; i < waves.size(); i++) {
		auto tmp = *waves[i];
		for (int j = 0; j < i; j++) {
			//auto scaler = (*(waves[i])).innerProduct(*(waves[j])) / (*(waves[j])).innerProduct(*(waves[j]));
			auto scaler1 = (*waves[i]).innerProduct(*waves[j]);
			auto scaler2 = (*waves[j]).innerProduct(*waves[j]);
			tmp = tmp - (*(waves[j])) * (scaler1 / scaler2);
		}
		*waves[i] = tmp;
		//for (int j = 0; j < waves.size(); j++) {
		//	std::cout << (*waves[j]).norm() << std::endl;
		//}
	}
}

inline void schmidtOrthonForCouplingWavesInSH(const std::vector< std::vector< Wavefunction1D* > >& coupling_waves)
{
	for (int i = 0; i < coupling_waves.size(); i++) {
		auto next_wave = coupling_waves[i];
		for (int j = 0; j < i; j++) {
			//inner product for scaler -> yield a wavefunction! not a complex number.
			Wavefunction1D scaler1(next_wave.front()->getGrid());
			Wavefunction1D scaler2(next_wave.front()->getGrid());
			for (int k = 0; k < next_wave.size(); k++) {
				scaler1 = scaler1 + (*coupling_waves[i][k]) * (*coupling_waves[j][k]).conjugate();
				scaler2 = scaler2 + (*coupling_waves[j][k]) * (*coupling_waves[j][k]).conjugate();
			}
			for (int k = 0; k < next_wave.size(); k++) {
				(*next_wave[k]) = (*next_wave[k]) - (*coupling_waves[j][k]) * (scaler1 / scaler2);
			}
		}

		//copy
		for (int k = 0; k < next_wave.size(); k++) {
			(*coupling_waves[i][k]) = (*next_wave[k]);
		}
	}
}



#endif // !1
#ifndef __WAVE_FUNCTION_H__
#define __WAVE_FUNCTION_H__

#include <vector>
#include <complex>
#include <cmath>
#include "NumericalGrid.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <fftw3.h>
#include "Util.h"
#include <iostream>

template<unsigned N>
class Wavefunction
{
public:
	Wavefunction();
	Wavefunction(const NumericalGrid<N>&);
	Wavefunction(const NumericalGrid<N>&, const std::vector<std::complex<double> >&);
	Wavefunction(const NumericalGrid<N>&, const Eigen::VectorXcd&);

	Wavefunction(const Wavefunction<N>&);
	Wavefunction(Wavefunction<N>&&) noexcept;
	Wavefunction<N>& operator=(const Wavefunction<N>&);
	Wavefunction<N>& operator=(Wavefunction<N>&&) noexcept;

	const NumericalGrid<N>& getGrid() const;
	//std::complex<double> getValue(const GridPoint<N>&) const;
	std::complex<double> getValueByIndex(int) const;
	const Eigen::VectorXcd& getSamplesView() const;
	Eigen::VectorXcd& getSamplesHandler();
	std::complex<double>& operator[](int);

	template<unsigned _N>
	Wavefunction<_N> getSlice(const std::vector< GridPoint<N> >&) const;

	Wavefunction<N> replaceGrid(const NumericalGrid<N>&);

	std::complex<double> norm() const;
	std::complex<double> innerProduct(const Wavefunction<N>&) const;
	Wavefunction<N> operator+(const Wavefunction<N>&) const;
	Wavefunction<N> operator-(const Wavefunction<N>&) const;
	Wavefunction<N> operator*(std::complex<double>) const;
	Wavefunction<N> operator*(const Wavefunction<N>&) const;

private:
	NumericalGrid<N> grid;
	Eigen::VectorXcd samples;
};

typedef Wavefunction<1> Wavefunction1D;
typedef Wavefunction<2> Wavefunction2D;
typedef Wavefunction<3> Wavefunction3D;


template<unsigned N>
Wavefunction<N>::Wavefunction()
	: grid(NumericalGrid<N>())
{
}

template<unsigned N>
Wavefunction<N>::Wavefunction(const NumericalGrid<N>& grid_)
	: grid(grid_), samples(Eigen::VectorXcd::Zero(grid_.getTotalCount()))
{
}

template<unsigned N>
Wavefunction<N>::Wavefunction(const NumericalGrid<N>& grid_, const std::vector<std::complex<double>>& data)
	: grid(grid_), samples(Eigen::VectorXcd::Zero(grid_.getTotalCount()))
{
	assert(grid.getTotalCount() == data.size());
	for (int i = 0; i < samples.size(); i++) {
		samples(i) = data[i];
	}
}

template<unsigned N>
Wavefunction<N>::Wavefunction(const NumericalGrid<N>& grid_, const Eigen::VectorXcd& samples_)
	: grid(grid_), samples(samples_)
{
	assert(grid_.getTotalCount() == samples_.size());
}

template<unsigned N>
Wavefunction<N>::Wavefunction(const Wavefunction<N>& wave)
	: grid(wave.grid), samples(wave.samples)
{
}

template<unsigned N>
Wavefunction<N>::Wavefunction(Wavefunction<N>&& wave) noexcept
	: grid(std::move(wave.grid)), samples(std::move(wave.samples))
{
}

template<unsigned N>
Wavefunction<N>& Wavefunction<N>::operator=(const Wavefunction<N>& wave)
{
	grid = wave.grid;
	samples = wave.samples;
	return *this;
}

template<unsigned N>
inline Wavefunction<N>& Wavefunction<N>::operator=(Wavefunction<N>&& wave) noexcept
{
	grid = std::move(wave.grid);
	samples = std::move(wave.samples);
	return *this;
}

template<unsigned N>
inline const NumericalGrid<N>& Wavefunction<N>::getGrid() const
{
	return grid;
}

template<unsigned N>
inline std::complex<double> Wavefunction<N>::getValueByIndex(int id) const
{
	return samples(id);
}

template<unsigned N>
inline const Eigen::VectorXcd& Wavefunction<N>::getSamplesView() const
{
	return samples;
}

template<unsigned N>
inline Eigen::VectorXcd& Wavefunction<N>::getSamplesHandler()
{
	return samples;
}

template<unsigned N>
inline std::complex<double>& Wavefunction<N>::operator[](int id)
{
	std::complex<double>& ref = samples(id);
	return ref;
}

template<unsigned N>
inline std::complex<double> Wavefunction<N>::norm() const
{
	return innerProduct(*this);
}

template<unsigned N>
inline std::complex<double> Wavefunction<N>::innerProduct(const Wavefunction<N>& wave) const
{
	assert(grid == wave.getGrid());
	std::complex<double> res = wave.samples.adjoint() * samples;
	std::complex<double> delta = 1;
	for (int i = 0; i < N; i++) {
		delta *= grid.getDelta(i);
	}
	return res * delta;
}

template<unsigned N>
inline Wavefunction<N> Wavefunction<N>::operator+(const Wavefunction<N>& wave) const
{
	assert(grid == wave.grid);
	auto res = samples + wave.samples;
	Wavefunction<N> new_wave(grid, res);
	return new_wave;
}

template<unsigned N>
inline Wavefunction<N> Wavefunction<N>::operator-(const Wavefunction<N>& wave) const
{
	assert(grid == wave.grid);
	auto res = samples - wave.samples;
	Wavefunction<N> new_wave(grid, res);
	return new_wave;
}

template<unsigned N>
inline Wavefunction<N> Wavefunction<N>::operator*(std::complex<double> scaler) const
{
	auto res = samples * scaler;
	Wavefunction<N> new_wave(grid, res);
	return new_wave;
}

template<unsigned N>
inline Wavefunction<N> Wavefunction<N>::operator*(const Wavefunction<N>& wave) const
{
	assert(grid == wave.grid);
	auto res = (samples.array() * wave.samples.array()).matrix();
	Wavefunction<N> new_wave(grid, res);
	return new_wave;
}

template<unsigned N>
template<unsigned _N>
inline Wavefunction<_N> Wavefunction<N>::getSlice(const std::vector<GridPoint<N>>& sub_area) const
{
	return Wavefunction<_N>();
}

template<unsigned N>
inline Wavefunction<N> Wavefunction<N>::replaceGrid(const NumericalGrid<N>&)
{
	return Wavefunction<N>();
}

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

inline Wavefunction1D fft_agent_func_1d(std::vector<double>, const NumericalGrid1D&, const Wavefunction1D&, int);

template<unsigned N>
inline Wavefunction<N> fft(const Wavefunction<N>& wave, const std::vector<double>& k) {
	return fft_agent_func(k, wave.getGrid(), wave, FFTW_FORWARD);
}

template<unsigned N>
inline Wavefunction<N> ifft(const Wavefunction<N>& wave, const std::vector<double>& k) {
	return fft_agent_func(k, wave.getGrid(), wave, FFTW_BACKWARD);
}

inline Wavefunction<1> fft1D(const Wavefunction<1>& wave, const std::vector<double>& k) {
	return fft_agent_func_1d(k, wave.getGrid(), wave, FFTW_FORWARD);
}

inline Wavefunction1D fft_agent_func_1d(std::vector<double> center_k_vector, const NumericalGrid1D& grid, const Wavefunction1D& myself, int mode)
{
	assert(mode == FFTW_FORWARD || mode == FFTW_BACKWARD);
	using namespace std::literals;

	int cnt = grid.getCount();
	double sign = 1;
	double mid = grid.getOffset();
	NumericalGrid1D grid_(grid);
	Wavefunction1D origin(myself);

	if (mode == FFTW_FORWARD) {
		sign = -1;
	}
	else if (mode == FFTW_BACKWARD) {
		sign = 1;
	}

	// start FFT
	int N = grid_.getCount();
	double delta_x = grid_.getDelta();
	double Lx = delta_x * N;
	double Lk = 2 * PI * N / Lx;
	double center_k = center_k_vector[0];
	double center_x = mid;

	NumericalGrid1D new_grid = NumericalGrid1D(N, Lk - Lk / N, center_k);

	fftw_complex* in, * out;
	fftw_plan p;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_1d(N, in, out, mode, FFTW_ESTIMATE);

	int fft_index = 0, grid_index = 0;
	double crt_x = 0, crt_k = 0;
	std::complex<double> tmp;
	int split_index = (int)std::floor((double)N / 2.0);

	for (int i = 0; i < N; i++) {
		fft_index = i - split_index;
		if (i < split_index) {
			fft_index += N;
		}
		crt_x = origin.getGrid().index(i).x() - center_x;
		tmp = origin.getValueByIndex(i) * std::exp(sign * 1i * center_k * crt_x);
		in[fft_index][0] = tmp.real();
		in[fft_index][1] = tmp.imag();
	}

	fftw_execute(p);

	double real = 0, imag = 0;
	std::vector< std::complex<double> > grid_data(N, std::complex<double>(0, 0));
	for (int i = 0; i < N; i++) {
		grid_index = i + split_index;
		if ((N - i - 1) < split_index) {
			grid_index -= N;
		}
		crt_k = new_grid.index(i).x() - center_k;
		tmp = std::complex<double>(out[i][0], out[i][1]) * std::exp(sign * 1i * (crt_k + center_k) * center_x) / sqrt(N);
		grid_data[grid_index] = tmp;
	}

	//for (int i = 0; i < grid.getTotalCount(); i++) {
	//	std::cout << grid_data[i] << " ";
	//}
	//std::cout << std::endl;
	//getchar();

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);

	return Wavefunction1D(new_grid, grid_data);
}


template<unsigned N>
inline Wavefunction<N> fft_agent_func(std::vector<double> center_k_vector, const NumericalGrid<N>& grid, const Wavefunction<N>& wave, int mode)
{
	assert(mode == FFTW_FORWARD || mode == FFTW_BACKWARD);
	assert(center_k_vector.size() == N);
	using namespace std::literals;

	int cnt[N]; //x-fast index
	double xmid[N], kmid[N], length_k[N];
	int cnt_reverse[N];  //y-fast index
	for (int i = 0; i < N; i++) {
		cnt[i] = grid.getCount(i);
		xmid[i] = grid.getOffset(i);
		length_k[i] = 2 * PI / grid.getLength(i) * grid.getCount(i);
		kmid[i] = center_k_vector[i];
		cnt_reverse[N - 1 - i] = grid.getCount(i);
	}

	double sign = 1;
	if (mode == FFTW_FORWARD) {
		sign = -1;
	}
	else if (mode == FFTW_BACKWARD) {
		sign = 1;
	}

	NumericalGrid<N> new_grid = NumericalGrid<N>(cnt, length_k, kmid);
	std::vector< std::complex<double> > new_wave_data(grid.getTotalCount(), 0);
	fftw_complex* in, * out;
	fftw_plan p;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * grid.getTotalCount());
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * grid.getTotalCount());
	p = fftw_plan_dft(N, cnt_reverse, in, out, mode, FFTW_ESTIMATE);
	
	int fft_index[N], grid_index[N];
	for (int i = 0; i < grid.getTotalCount(); i++) {
		auto indice = grid.expand(i);
		for (int rank = 0; rank < N; rank++) {
			int split_index = (int)std::floor((double)cnt[rank] / 2.0);
			fft_index[rank] = indice[rank] - split_index;
			if (indice[rank] < split_index) {
				fft_index[rank] += cnt[rank];
			}
		}
		in[grid.shrink(fft_index)][0] = wave.getValueByIndex(i).real();
		in[grid.shrink(fft_index)][1] = wave.getValueByIndex(i).imag();
	}

	fftw_execute(p);

	for (int i = 0; i < grid.getTotalCount(); i++) {
		auto indice = grid.expand(i);
		for (int rank = 0; rank < N; rank++) {
			int split_index = (int)std::floor((double)cnt[rank] / 2.0);
			grid_index[rank] = indice[rank] + split_index;
			if ((cnt[rank] - 1 - indice[rank]) < split_index) {
				grid_index[rank] -= cnt[rank];
			}
		}
		new_wave_data[grid.shrink(grid_index)] = std::complex<double>(out[i][0], out[i][1]) / sqrt(grid.getTotalCount());
		//std::cout << "id: " << grid.shrink(grid_index) << std::endl;
	}

	//for (int i = 0; i < grid.getTotalCount(); i++) {
	//	//std::cout << std::complex<double>(out[i][0], out[i][1]) << " ";
	//	std::cout << new_wave_data[i] << " ";
	//}
	//std::cout << std::endl;
	//getchar();

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);

	return Wavefunction<N>(new_grid, new_wave_data);
}

#endif // !__WAVE_FUNCTION_H__
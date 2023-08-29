#ifndef __EXPORT_DATA_TO_MAT_FILE_H__
#define __EXPORT_DATA_TO_MAT_FILE_H__

#include <vector>
#include "Wavefunction.h"

#define NORM_MODE 1
#define REAL_MODE 2
#define IMAG_MODE 3

template<unsigned N>
using WaveMsg = std::tuple<Wavefunction<N>, std::string, int>;

int SingleWave1DToMat(const Wavefunction1D&, int);

int WaveAnimation1DToMat(const std::vector< WaveMsg<1> >&, int);

int WaveAnimation2DToMat(const std::vector< WaveMsg<2> >&, int);

int WaveAnimation3DToMat(const std::vector< WaveMsg<3> >&, int);

#endif // !__EXPORT_DATA_TO_MAT_FILE_H__

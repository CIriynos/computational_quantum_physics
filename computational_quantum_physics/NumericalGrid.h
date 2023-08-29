#ifndef __NUMERICAL_GRID_H__
#define __NUMERICAL_GRID_H__

#include <assert.h>

constexpr double MAX_COMPARISON_ERROR = 1e-4;

template<unsigned N>
class GridPoint {
public:

	GridPoint() {
		for (int i = 0; i < N; i++) {
			points[i] = 0;
		}
	}

	template<typename... Args>
	GridPoint(Args... parameters)
		: parameter_count(sizeof...(parameters)) {
		assert(sizeof...(parameters) == N);
		init(parameters...);
	}

	GridPoint(const GridPoint<N>& p) {
		for (int i = 0; i < N; i++) {
			points[i] = p.points[i];
		}
	}

	GridPoint<N>& operator=(const GridPoint<N>& p) {
		for (int i = 0; i < N; i++) {
			points[i] = p.points[i];
		}
		return *this;
	}

	friend bool operator==(const GridPoint<N>& a, const GridPoint<N>& b) {
		bool flag = true;
		for (int i = 0; i < N; i++) {
			flag = (a.points[i] == b.points[i]) && flag;
		}
		return flag;
	}

	bool isApproxEqual(const GridPoint<N>& p) {
		double sum = 0;
		for (int i = 0; i < N; i++) {
			sum += std::abs(points[i] - p.points[i]);
		}
		return sum < MAX_COMPARISON_ERROR;
	}

	double& get(int degree) {
		assert(degree < N);
		return points[degree];
	}

	double& x() {
		assert(N >= 1);
		return points[0];
	}

	double& y() {
		assert(N >= 2);
		return points[1];
	}

	double& z() {
		assert(N >= 3);
		return points[2];
	}

	friend std::ostream& operator<<(std::ostream& os, const GridPoint<N>& p) {
		os << "(";
		for (int i = 0; i < N; i++) {
			os << p.points[i];
			if (i != N - 1) {
				os << ", ";
			}
		}
		os << ")";
		return os;
	}

private:
	void init(double num) {
		points[N - 1] = num;
	}

	template<typename... Args>
	void init(double num, Args... rest) {
		int init_update_pos = parameter_count - (sizeof...(rest) + 1);
		points[init_update_pos] = num;
		init(rest...);
	}

	int parameter_count;
	double points[N];
};

typedef GridPoint<1> Point1D;
typedef GridPoint<2> Point2D;
typedef GridPoint<3> Point3D;

constexpr int LITTEL_INDEX_FASTRUNING = 1;
constexpr int LARGE_INDEX_FASTRUNING = 2;

//DO NOT USE M = LARGE_INDEX_FASTRUNING!!!!!!!!

template<unsigned N, unsigned M = LITTEL_INDEX_FASTRUNING>
class NumericalGrid
{
public:
	template<typename... Args>
	NumericalGrid(Args... parameters)
		: parameter_count(static_cast<int>(sizeof...(parameters) / 3))
	{
		assert(sizeof...(parameters) == (3 * N));
		init(parameters...);
		update_total_cnt();
	}

	NumericalGrid(int(&cnt)[N], double(&len)[N], double(&off)[N]) : parameter_count(N)
	{
		for (int i = 0; i < N; i++) {
			pos_count[i] = cnt[i];
			length[i] = len[i];
			offset[i] = off[i];
		}
		update_total_cnt();
	}

	NumericalGrid() : parameter_count(N) {
		for (int i = 0; i < N; i++) {
			pos_count[i] = 0;
			length[i] = 0.0;
			offset[i] = 0.0;
		}
		update_total_cnt();
	}

	NumericalGrid(const NumericalGrid<N>& grid) {
		for (int i = 0; i < N; i++) {
			pos_count[i] = grid.pos_count[i];
			length[i] = grid.length[i];
			offset[i] = grid.offset[i];
		}
		update_total_cnt();
	}

	NumericalGrid<N>& operator=(const NumericalGrid<N>& grid) {
		for (int i = 0; i < N; i++) {
			pos_count[i] = grid.pos_count[i];
			length[i] = grid.length[i];
			offset[i] = grid.offset[i];
		}
		update_total_cnt();
		return *this;
	}

	friend bool operator==(const NumericalGrid<N>& a, const NumericalGrid<N>& b) {
		bool flag = true;
		for (int i = 0; i < N; i++) {
			flag = (a.pos_count[i] == b.pos_count[i]
				&& a.length[i] == b.length[i]
				&& a.offset[i] == b.offset[i]) && flag;
		}
		return flag;
	}

	GridPoint<N> index(int id) const { return index_agent<M>(id); }

	std::vector<int> expand(int id) const 
	{
		assert(id >= 0 && id < total_count);
		std::vector<int> ans(N, 0);
		for (int j = 0; j < N; j++) {
			int num = id % pos_count[j];
			ans[j] = num;
			id = (id - num) / pos_count[j];
		}
		return ans;
	}

	int shrink(const std::vector<int>& indices) const {
		assert(indices.size() == N);
		int id = 0, tmp = 1;
		for (int j = 0; j < N; j++) {
			if (indices[j] < 0 || indices[j] >= pos_count[j]) {
				return -1;
			}
			id += indices[j] * tmp;
			tmp *= pos_count[j];
		}
		return id;
	}

	int shrink(int (&indices)[N]) const {
		int id = 0, tmp = 1;
		for (int j = 0; j < N; j++) {
			if (indices[j] < 0 || indices[j] >= pos_count[j]) {
				return -1;
			}
			id += indices[j] * tmp;
			tmp *= pos_count[j];
		}
		return id;
	}

	int getTotalCount() const { return total_count; }
	int getCount(int degree = 0) const { return pos_count[degree]; }
	double getLength(int degree = 0) const { return length[degree]; }
	double getOffset(int degree = 0) const { return offset[degree]; }
	double getDelta(int degree = 0) const { return length[degree] / pos_count[degree]; }

	friend std::ostream& operator<<(std::ostream& os, const NumericalGrid<N>& grid) {
		os << "pos_count: ";
		for (int i = 0; i < N; i++) {
			os << grid.pos_count[i] << " ";
		}
		os << std::endl;

		os << "length: ";
		for (int i = 0; i < N; i++) {
			os << grid.length[i] << " ";
		}
		os << std::endl;

		os << "offset: ";
		for (int i = 0; i < N; i++) {
			os << grid.offset[i] << " ";
		}
		os << std::endl;
	}

private:
	double get_value_by_index(int index, int degree) const {
		return offset[degree] - length[degree] / 2 + index * (length[degree] / pos_count[degree]);
	}

	void update_total_cnt() {
		total_count = 1;
		for (int i = 0; i < N; i++) {
			total_count *= pos_count[i];
		}
	}

	void init(int num, double len, double off) {
		pos_count[N - 1] = num;
		length[N - 1] = len;
		offset[N - 1] = off;
	}

	template<typename... Args>
	void init(int num, double len, double off, Args... rest) {
		int init_update_pos = parameter_count - (static_cast<int>(sizeof...(rest) / 3) + 1);
		pos_count[init_update_pos] = num;
		length[init_update_pos] = len;
		offset[init_update_pos] = off;
		init(rest...);
	}

	template<unsigned __M>
	GridPoint<N> index_agent(int id) const {
		assert(__M == LITTEL_INDEX_FASTRUNING || __M == LARGE_INDEX_FASTRUNING);
		return GridPoint<N>();
	}

	template<>
	GridPoint<N> index_agent<LITTEL_INDEX_FASTRUNING>(int id) const {
		assert(id >= 0 && id < total_count);
		GridPoint<N> ans;
		for (int j = 0; j < N; j++) {
			int num = id % pos_count[j];
			ans.get(j) = get_value_by_index(num, j);
			id = (id - num) / pos_count[j];
		}
		return ans;
	}

	template<>
	GridPoint<N> index_agent<LARGE_INDEX_FASTRUNING>(int id) const {
		assert(id >= 0 && id < total_count);
		GridPoint<N> ans;
		for (int j = N - 1; j >= 0; j--) {
			int num = id % pos_count[j];
			ans.get(j) = get_value_by_index(num, j);
			id = (id - num) / pos_count[j];
		}
		return ans;
	}

private:
	int pos_count[N];
	double length[N];
	double offset[N];
	int parameter_count;
	int total_count;
};

typedef NumericalGrid<1> NumericalGrid1D;
typedef NumericalGrid<2> NumericalGrid2D;
typedef NumericalGrid<3> NumericalGrid3D;


#endif // !__NUMERICAL_GRID_H__
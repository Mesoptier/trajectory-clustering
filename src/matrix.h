#pragma once

#include <vector>

template <typename T>
class Matrix
{
private:
	using Container = typename std::vector<T>;
	Container elements;

public:
	using Size = unsigned int;
	Size const n, m;

	struct Index {
	private:
		Size _n, _m;
		Size _i = 0;
		Size _j = 0;

	public:
		Index(Matrix const& matrix) : _n(matrix.n), _m(matrix.m) {}

		Size i() const { return _i; }
		Size j() const { return _j; }
		bool valid() const { return _i < _n && _j < _m; }
		Index const& operator++() {
			if (_i == _n-1) { _i = 0; ++_j; }
			else { ++_i; }
			return *this;
		}
		void reset() { _i = 0; _j = 0; }
	};

	Matrix(Size n, Size m, T default_value = T())
		: elements(n*m, default_value), n(n), m(m) {}

	T& operator()(Size i, Size j) { return elements[i + j*n]; }
	T& operator[](Index index) { return elements[index.i() + index.j()*n]; }
	T const& operator()(Size i, Size j) const { return elements[i + j*n]; }
	T const& operator[](Index index) const { return elements[index.i() + index.j()*n]; }
	typename Container::iterator begin() { return elements.begin(); }
	typename Container::iterator end() { return elements.end(); }
	typename Container::const_iterator begin() const { return elements.cbegin(); }
	typename Container::const_iterator end() const { return elements.cend(); }
};

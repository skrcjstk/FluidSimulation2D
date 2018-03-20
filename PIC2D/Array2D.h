#pragma once
/*
Copyright (c) 2016 Ryan L. Guy

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
claim that you wrote the original software. If you use this software
in a product, an acknowledgement in the product documentation would be
appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef __ARRAY2D_H__
#define __ARRAY2D_H__

#include <vector>
#include <stdio.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>

namespace PICArray2d
{
	struct GridIndex {
		int i, j;

		GridIndex() : i(0), j(0) {}
		GridIndex(int ii, int jj) : i(ii), j(jj) {}

		bool operator==(const GridIndex &other) const {
			return i == other.i &&
				j == other.j;
		}

		bool operator!=(const GridIndex &other) const {
			return i != other.i ||
				j != other.j;
		}

		int& operator[](unsigned int idx) {
			if (idx > 1) {
				std::string msg = "Error: index out of range.\n";
				throw std::out_of_range(msg);
			}

			return (&i)[idx];
		}
	};

	template <class T>
	class Array2d
	{
	public:
		Array2d() {
			_initializeGrid();
		}

		Array2d(int i, int j) : width(i), height(j), _numElements(i*j) {
			_initializeGrid();
		}

		Array2d(int i, int j, T fillValue) : width(i), height(j),  _numElements(i*j) {
			_initializeGrid();
			fill(fillValue);
		}

		Array2d(const Array2d &obj) {
			width = obj.width;
			height = obj.height;
			_numElements = obj._numElements;

			_initializeGrid();

			T val;
			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {
					val = obj._grid[_getFlatIndex(i, j)];
					set(i, j, val);
				}
			}

			if (obj._isOutOfRangeValueSet) {
				_outOfRangeValue = obj._outOfRangeValue;
				_isOutOfRangeValueSet = true;
			}
		}

		Array2d operator=(const Array2d &rhs) {
			delete[] _grid;

			width = rhs.width;
			height = rhs.height;
			_numElements = rhs._numElements;

			_initializeGrid();

			T val;
			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {
					val = rhs._grid[_getFlatIndex(i, j)];
					set(i, j, val);
				}
			}

			if (rhs._isOutOfRangeValueSet) {
				_outOfRangeValue = rhs._outOfRangeValue;
				_isOutOfRangeValueSet = true;
			}

			return *this;
		}

		~Array2d() {
			delete[] _grid;
		}

		void fill(T value) {
			for (int idx = 0; idx < width*height; idx++) {
				_grid[idx] = value;
			}
		}

		T operator()(int i, int j) {
			bool isInRange = _isIndexInRange(i, j);
			if (!isInRange && _isOutOfRangeValueSet) {
				return _outOfRangeValue;
			}

			if (!isInRange) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(i) + " j: " + _toString(j) + "\n";
				throw std::out_of_range(msg);
			}

			return _grid[_getFlatIndex(i, j)];
		}

		T operator()(GridIndex g) {
			bool isInRange = _isIndexInRange(g.i, g.j);
			if (!isInRange && _isOutOfRangeValueSet) {
				return _outOfRangeValue;
			}

			if (!isInRange) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(g.i) + " j: " + _toString(g.j) + "\n";
				throw std::out_of_range(msg);
			}

			return _grid[_getFlatIndex(g)];;
		}

		T operator()(int flatidx) {
			bool isInRange = flatidx >= 0 && flatidx < _numElements;
			if (!isInRange && _isOutOfRangeValueSet) {
				return _outOfRangeValue;
			}

			if (!isInRange) {
				std::string msg = "Error: index out of range.\n";
				msg += "index: " + _toString(flatidx) + "\n";
				throw std::out_of_range(msg);
			}

			return _grid[flatidx];
		}

		T get(int i, int j) {
			bool isInRange = _isIndexInRange(i, j);
			if (!isInRange && _isOutOfRangeValueSet) {
				return _outOfRangeValue;
			}

			if (!isInRange) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(i) + " j: " + _toString(j) + "\n";
				throw std::out_of_range(msg);
			}

			return _grid[_getFlatIndex(i, j)];
		}

		T get(GridIndex g) {
			bool isInRange = _isIndexInRange(g.i, g.j);
			if (!isInRange && _isOutOfRangeValueSet) {
				return _outOfRangeValue;
			}

			if (!isInRange) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(g.i) + " j: " + _toString(g.j) + "\n";
				throw std::out_of_range(msg);
			}

			return _grid[_getFlatIndex(g)];;
		}

		T get(int flatidx) {
			bool isInRange = flatidx >= 0 && flatidx < _numElements;
			if (!isInRange && _isOutOfRangeValueSet) {
				return _outOfRangeValue;
			}

			if (!isInRange) {
				std::string msg = "Error: index out of range.\n";
				msg += "index: " + _toString(flatidx) + "\n";
				throw std::out_of_range(msg);
			}

			return _grid[flatidx];
		}

		void set(int i, int j, T value) {
			if (!_isIndexInRange(i, j)) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(i) + " j: " + _toString(j) + "\n";
				throw std::out_of_range(msg);
			}

			_grid[_getFlatIndex(i, j)] = value;
		}

		void set(GridIndex g, T value) {
			if (!_isIndexInRange(g)) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(g.i) + " j: " + _toString(g.j) + "\n";
				throw std::out_of_range(msg);
			}

			_grid[_getFlatIndex(g)] = value;
		}

		void set(std::vector<GridIndex> &cells, T value) {
			for (unsigned int i = 0; i < cells.size(); i++) {
				set(cells[i], value);
			}
		}

		void set(int flatidx, T value) {
			if (!(flatidx >= 0 && flatidx < _numElements)) {
				std::string msg = "Error: index out of range.\n";
				msg += "index: " + _toString(flatidx) + "\n";
				throw std::out_of_range(msg);
			}

			_grid[flatidx] = value;
		}

		void add(int i, int j, T value) {
			if (!_isIndexInRange(i, j)) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(i) + " j: " + _toString(j) + "\n";
				throw std::out_of_range(msg);
			}

			_grid[_getFlatIndex(i, j)] += value;
		}

		void add(GridIndex g, T value) {
			if (!_isIndexInRange(g)) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(g.i) + " j: " + _toString(g.j) + "\n";
				throw std::out_of_range(msg);
			}

			_grid[_getFlatIndex(g)] += value;
		}

		void add(int flatidx, T value) {
			if (!(flatidx >= 0 && flatidx < _numElements)) {
				std::string msg = "Error: index out of range.\n";
				msg += "index: " + _toString(flatidx) + "\n";
				throw std::out_of_range(msg);
			}

			_grid[flatidx] += value;
		}

		T *getPointer(int i, int j) {
			bool isInRange = _isIndexInRange(i, j);
			if (!isInRange && _isOutOfRangeValueSet) {
				return &_outOfRangeValue;
			}

			if (!isInRange) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(i) + " j: " + _toString(j) + "\n";
				throw std::out_of_range(msg);
			}

			return &_grid[_getFlatIndex(i, j)];
		}

		T *getPointer(GridIndex g) {
			bool isInRange = _isIndexInRange(g.i, g.j);
			if (!isInRange && _isOutOfRangeValueSet) {
				return &_outOfRangeValue;
			}

			if (!isInRange) {
				std::string msg = "Error: index out of range.\n";
				msg += "i: " + _toString(g.i) + " j: " + _toString(g.j) + "\n";
				throw std::out_of_range(msg);
			}

			return &_grid[_getFlatIndex(g)];
		}

		T *getPointer(int flatidx) {
			bool isInRange = flatidx >= 0 && flatidx < _numElements;
			if (!isInRange && _isOutOfRangeValueSet) {
				return &_outOfRangeValue;
			}

			if (!isInRange) {
				std::string msg = "Error: index out of range.\n";
				msg += "index: " + _toString(flatidx) + "\n";
				throw std::out_of_range(msg);
			}

			return &_grid[flatidx];
		}

		T *getRawArray() {
			return _grid;
		}

		int getNumElements() {
			return _numElements;
		}

		void setOutOfRangeValue() {
			_isOutOfRangeValueSet = false;
		}
		void setOutOfRangeValue(T val) {
			_outOfRangeValue = val;
			_isOutOfRangeValueSet = true;
		}

		bool isOutOfRangeValueSet() {
			return _isOutOfRangeValueSet;
		}
		T getOutOfRangeValue() {
			return _outOfRangeValue;
		}

		inline bool isIndexInRange(int i, int j) {
			return i >= 0 && j >= 0 && i < width && j < height;
		}

		inline bool isIndexInRange(GridIndex g) {
			return g.i >= 0 && g.j >= 0 && g.i < width && g.j < height;
		}

		int width = 0;
		int height = 0;

	private:
		void _initializeGrid() {
			if (width < 0 || height < 0) {
				std::string msg = "Error: dimensions cannot be negative.\n";
				msg += "width: " + _toString(width) +
					"height: " + _toString(height) + "\n";
				throw std::domain_error(msg);
			}

			_grid = new T[width*height];
		}

		inline bool _isIndexInRange(int i, int j) {
			return i >= 0 && j >= 0 && i < width && j < height;
		}

		inline bool _isIndexInRange(GridIndex g) {
			return g.i >= 0 && g.j >= 0 && g.i < width && g.j < height;
		}

		inline unsigned int _getFlatIndex(int i, int j) {
			return (unsigned int)i + (unsigned int)width * ((unsigned int)j);
		}

		inline unsigned int _getFlatIndex(GridIndex g) {
			return (unsigned int)g.i + (unsigned int)width * ((unsigned int)g.j);
		}

		template<class S>
		std::string _toString(S item) {
			std::ostringstream sstream;
			sstream << item;

			return sstream.str();
		}

		T *_grid;

		bool _isOutOfRangeValueSet = false;
		T _outOfRangeValue;
		int _numElements = 0;
	};
}
#endif

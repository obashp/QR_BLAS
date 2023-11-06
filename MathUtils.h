#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <fstream>
#include <vector>

extern "C"{
	extern int dswap_(int* N, double *DX, int* INCX, double *DY, int *INCY);
	extern int dscal_(int* N, double *DA, double *DX, int* INCX);
	extern int dcopy_(int* N, double *DX, int* INCX, double *DY, int *INCY);
	extern int daxpy_(int* N, double *DA, double *DX, int *INCX, double *DY, int *INCY);
	extern double ddot_(int*N, double *DX, int *INCX, double *DY, int *INCY);
	extern double dnrm2_(int *N, double *DX, int *INCX);
}

namespace blas
{
	int SWAP(int N, double *DX, double *DY);
	int SWAP(std::vector<double> DX, std::vector<double> DY);

	int SCALE(double *Z, int N, double Alpha, double *X);
	int SCALE(const int N, const double Alpha, double *X);
	int SCALE(std::vector<double> Z, double Alpha, std::vector<double> X);
	int SCALE(const int N, const double Alpha);
}


/**
 * class Matrix - It is the class which defines a dense matrix.
 * 
 */ 
class Matrix{
	private:
		// Size of matrix
		size_t M, N;	

		// Array of matrix entries stored in a 1-D array
		double *MATRIX;

		// Vector of pointers to column start locations in MATRIX
		// We use the vector in conjunction with iterators, to return
		// a sub-column of the matrix
		std::vector<double *> columns;

	public:
		// Default constructor
		Matrix(): M(0), N(0), MATRIX(nullptr){}

		// Overloaded constructor
		Matrix(size_t iM, size_t iN);

		// Copy constructor
		Matrix(const Matrix &A);

		// Assignment operator overload
		Matrix operator=(const Matrix &A);

		// Overload round brackets to access index at i,j
		/// The first operator overload allows for accessing the value at i,j and 
		/// returns a value which can be modified. The second operator overload returns
		/// a const value.
		inline double& operator()(unsigned int i, unsigned int j) {return columns[j][i];}
		inline double operator()(unsigned int i, unsigned int j) const {return columns[j][i];}

		// getSize - returns the dimensions of the matrix
		inline const void getSize(size_t &outM, size_t &outN) const{outM = M; outN = N;}

		// getM - return the number of rows of the matrix
		inline const size_t& getM() const {return M;}

		// getN - return the number of columns of the matrix
		inline const size_t& getN() const {return N;}

		// setColumn - has different overloads for different input types and returns whether 
		// the column has been successfully set in the matrix.
		// Variants:
		// 1. setColumn(unsigned int columnIdx, double* column, int N)
		// 2. setColumn(unsigned int columnIdx, vector column)
		// 3. setColumn(unsigned int columnIdx, double* column, unsigned int startIdx, unsigned int endIdx)
		// 4. setColumn(unsigned int columnIdx, vector column, unsigned int startIdx, unsigned int endIdx)
		// 5. setColumn(unsigned int columnIdx, vector iterator start, vector iterator end, unsigned int startIdx, unsigned int endIdx)
		int setColumn(const unsigned int columnIdx, const double* column, const int N);
		int setColumn(const unsigned int columnIdx, const std::vector<double>& column);
		int setColumn(const unsigned int columnIdx, const double* column,const unsigned int idx0, const unsigned int idx1);
		int setColumn(const unsigned int columnIdx, const std::vector<double>& column, const unsigned int idx0, const unsigned int idx1);
		int setColumn(const unsigned int columnIdx, const std::vector<double>::iterator colStart, const std::vector<double>::iterator colEnd,
													const unsigned int idx0, const unsigned int idx1 );


};
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>
#include <fstream>
#include <ctime>

double getRand()
{
	int N = 1001;
	double dx = 2.0/(N-1);

	int i = rand()%(N-1);
	// std::cout<<i<<std::endl;

	double x = -1.0 + i*dx;

	return exp(-x*x);

}

void copyVector(double *out, const double* in, const int N)
{
	for(int i = 0; i < N; i++)
		out[i] = in[i];
}

void scaleVector(double *out, const double Alpha, const double* in, const int N)
{
	for(int i = 0; i < N; i++)
		out[i] = Alpha*in[i];
}

void saxpy(double *z, double a, double *x, double *y, int M)
{
	for(int i = 0; i < M; i++)
		z[i] = a*x[i] + y[i];
}

double dotVector(const double *in1, const double *in2, const int N)
{
	double dotP = 0.0;
	for(int i = 0; i < N; i++)
		dotP += (in1[i]*in2[i]);

	return dotP;
}

double norm2Vec(const double *in1, int N)
{
	double norm = 0.0;
	norm = dotVector(in1, in1, N);
	norm = sqrt(norm);
	return norm;
}

class Matrix{
	private:
		int M, N;			/// Size of matrix
		double *MATRIX;		/// Matrix entries
		double **columns;	/// Pointer to column start locations

	public:
		Matrix():M(0), N(0), MATRIX(nullptr), columns(nullptr){}
		Matrix(int iM, int iN);

		Matrix(Matrix const &A);	/// Copy constructor
		Matrix operator=(Matrix const &A);	/// Equality operator
		~Matrix();

		void resizeMatrix(int iM, int iN);
		int setColumn(int iCol, double *column, int colM);
		int setRow(int iRow, double *row, int rowN);
		void saveMatrix(std::string fileName);

		void getColumn(double *column, int iCol);
		inline int getM(){return M;}
		inline int getN(){return N;}
};



void Matrix::getColumn(double *column, int iCol)
{
	column = columns[iCol];
}

Matrix::Matrix(int iM, int iN) : M(iM), N(iN)
{
	MATRIX = (double *) malloc(M*N*sizeof(double));
	columns = (double **) malloc(N*sizeof(double*));

	/// Assign columns the column pointers of MATRIX
	/// Initialize the matrix with the NxN identity
	for(int i = 0; i < N; i++){
		columns[i] = &MATRIX[i*M];
		for(int j = 0; j < M; j++){
			MATRIX[i*M+j] = (i == j ? 1.0:0.0);
		}
	}
}	

Matrix::~Matrix()
{
	for(int i = 0; i < N; i++)
		columns[i] = nullptr;
	free(columns);
	free(MATRIX);
	MATRIX = nullptr;
	M = 0; N = 0;
}

int Matrix::setColumn(int iCol, double *column, int colM)
{
	int INFO = 0;
	if(M == 0 || N == 0)
		INFO = -1;
	else if(iCol > N)
		INFO = -2;
	else if(MATRIX == nullptr)
		INFO = -3;
	else if(colM != M)
		INFO = -4;

	if(INFO != 0)
		return INFO;

	copyVector(columns[iCol], column, M);

	return INFO;
}

int Matrix::setRow(int iRow, double *row, int rowN)
{
	int INFO = 0;
	if(M == 0 || N == 0)
		INFO = -1;
	else if(iRow > M)
		INFO = -2;
	else if(MATRIX == nullptr)
		INFO = -3;
	else if(rowN != N)
		INFO = -4;

	if(INFO != 0)
		return INFO;

	for(int i = 0; i < N; i++)
		*(columns[i]+iRow) = row[i];

	return INFO;
}

void Matrix::saveMatrix(std::string fileName)
{
	std::fstream FILE;
	FILE.open(fileName,std::ios::out);

	for(int iRow = 0; iRow < M; iRow++){
		for(int iCol = 0; iCol < N; iCol++){
			FILE<<MATRIX[iCol*M+iRow]<<"	";
		}
		FILE<<std::endl;
	}
	FILE.close();
}

void Matrix::resizeMatrix(int iM, int iN)
{
	if(MATRIX != nullptr){
		free(MATRIX);
		MATRIX = nullptr;
	}

	for(int i = 0; i < N; i++)
		columns[i] = nullptr;
	free(columns);

	M = iM, N = iN;
	MATRIX = (double *) malloc(M*N*sizeof(double));
	columns = (double **) malloc(N*sizeof(double*));

	for(int i = 0; i < N; i++)
		columns[i] = &MATRIX[i*M];
}

void computeHouseHolder(double *v, double &beta, double *x, int N)
{
	double sigma = 0.0, mu = 0.0;
	for(int i = 1; i < N; i++){
		sigma += (x[i]*x[i]);
		v[i] = x[i];
	}
	v[0] = 1.0;

	if(sigma == 0.0 && v[0] >= 0.0){
		beta = 0.0;
	} else if (sigma == 0.0 && v[0] < 0.0){
		beta = -2.0;
	} else{
		mu = sqrt(column[0]*column[0] + sigma);
		if(column[0] <= 0.0)
			v[0] = column[0] - mu;
		else
			v[0] = -sigma/(column[0] + mu);
		beta = 2.0*v[0]*v[0]/(v[0]*v[0] + sigma);
		for(int i = 1; i < N; i++)
			v[i] = v[i]/v[0];
		v[0] = 1.0;
	}
}

void getQR(Matrix &A)
{
	int M, N;
	M = A.getM(); N = A.getN();
	double *V, *X, *H;
	double beta;
	V = (double*) malloc(M*sizeof(double));
	X = (double*) malloc(M*sizeof(double));
	H = (double*) malloc(M*sizeof(double));


	/// Loop over the columns of A
	for(int iCol = 0; iCol < N; iCol++)
	{
		/// Get the ith column of A
		A.getColumn(V, iCol);

		/// Construct the householder vector from the part of V from iCol to M
		for(int i = iCol; i < M; i++)
			X[i - iCol] = V[i];

		computeHouseHolder(H, Beta, X, M-iCol);

		/// Apply the householder vector to columns of A, and create low-rank updates to A
		for(int jCol = iCol; jCol < N; jCol++){
			/// Get the jth column of A
			A.getColumn(V, jCol);	

			/// Restrict the part of the column to the length below the diagonal
			for(int i = jCol; i < M; i++)
				X[i - jCol] = V[i];

			/// Construct the Householder transformation for the jth column
			saxpy(X,-Beta*dotVector(H,X,M-jCol),H,X, M-jCol);

			/// Store the Householder transformed part back in A.
			for(int i = 0; i < jCol; j++)
				V[i] = 
		}



	}



}




int main()
{
	int INFO;
	Matrix A(10, 5);
	/// Set columns of matrix M
	srand((unsigned int)time(NULL));
	double column[10];
	for(int i = 0; i < 5; i++){
		for(int j = 0; j < 10; j++){
			column[j] = getRand();
			// std::cout<<column[j]<<std::endl;
		}
		INFO = A.setColumn(i, column, 10);
		if(INFO != 0){
			std::cout<<"Error code encountered in writing column. "<<INFO<<std::endl;
			break;
		}
	}


	/// Construct a new matrix QR, whose upper-triangular part contains R and the rows below 
	/// contain the house-holder reflectors
	Matrix Q = A;
	getQR(M);








	A.saveMatrix("M.txt");


}

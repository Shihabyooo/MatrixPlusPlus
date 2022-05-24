#include "MatrixPP_f32.hpp"


// Matrix_f32::Matrix_f32()
// {
//     Array2D<float>::Array2D();
// }

// Matrix_f32::Matrix_f32(_INDEX, _INDEX)
// {
    
// }

// Matrix_f32::Matrix_f32(const Matrix_f32 & sourceMAt)
// {

// }

// Matrix_f32::Matrix_f32(_INDEX _rows, _INDEX _columns, float defaultValue)
// {

// }

// Matrix_f32::Matrix_f32(std::vector<std::vector<float>> & sourceVec)
// {

// }

// Matrix_f32::~Matrix_f32()
// {

// }

Matrix_f32::Matrix_f32(Array2D<float> sourceArr)
{
	CopyFromArray2D(sourceArr);
}

Matrix_f32::Matrix_f32(Array2D<int> sourceArr)
{
    CopyFromArray2D(sourceArr);
}

Matrix_f32::Matrix_f32(Array2D<long> sourceArr)
{
    CopyFromArray2D(sourceArr);
}

Matrix_f32 Matrix_f32::operator*(const Matrix_f32 & mat2)
{
    return MultiplyMatricesVectorized(*this, mat2);
}

Matrix_f32 Matrix_f32::operator*(const double & scalar)
{
    return MultiplayMatrixWithScalarVectorized(*this, scalar);
}

Matrix_f32 Matrix_f32::operator+(const Matrix_f32 & mat2)
{
    return AddMatricesVectorized(*this, mat2);
}

Matrix_f32 Matrix_f32::operator-(const Matrix_f32 & mat2)
{
    return SubtractMatricesVectorized(*this, mat2);
}

Matrix_f32 Matrix_f32::GetSubMatrix(_INDEX beginRow, _INDEX noOfRows, _INDEX beginColumn, _INDEX noOfColumns)
{
    	
    #ifdef _USE_BOUNDS_CHECK
	if (content == NULL || beginRow + noOfRows > rows || beginColumn + noOfColumns > columns)  //TODO consider rearranging the arguments for this method to have the pivot coord
																			//first then lengths.
	{
		//std::cout << "ERROR! Attempting to extract submatrix with a range outside the original matrix bounds" << std::endl;
		throw std::out_of_range("ERROR! Column or Row value out of range or content set to NULL");
	}
    #endif

	if (noOfRows < 1 || noOfColumns < 1)
	{
		std::cout << "ERROR! Attempting to extract submatrix with null size in one both dimensions." << std::endl;
		return Matrix_f32();
	}
    
	Matrix_f32 subMatrix(noOfRows, noOfColumns);

	for (_INDEX i = beginRow; i < beginRow + noOfRows; i++)
	{
		for (_INDEX j = beginColumn; j < beginColumn + noOfColumns; j++)
		{
			subMatrix.SetValue(i - beginRow, j - beginColumn, content[i][j]);
		}
	}

	return subMatrix;
}

Matrix_f32 Matrix_f32::Invert()
{
    return InvertMatrix(*this);
}

double Matrix_f32::Determinant()
{
    return CalculateDeterminant(*this);
}

void Matrix_f32::Overlay(const Matrix_f32 mat2, _INDEX rowOffset, _INDEX columnOffset)
{
    if (content == NULL) //if this Array2D is uninitialized, this method acts as a simple assignment.
	{
		*this = mat2;
	}

	_INDEX _rows = rows > mat2.Rows() + rowOffset ? mat2.Rows() + rowOffset : rows;
	_INDEX _columns = columns > mat2.Columns() + columnOffset ? mat2.Columns() + columnOffset : columns;

	for (_INDEX i = rowOffset; i < _rows; i++)
	{
		for (_INDEX j = columnOffset; j < _columns; j++)
			content[i][j] += mat2.GetValue(i - rowOffset, j - columnOffset);
	}
}

Matrix_f32 ** Matrix_f32::DecomposeLU()
{
    return DecomposeLU(*this);
}

Matrix_f32 ** Matrix_f32::DecomposeLUP()
{
    return DecomposeLUP(*this);
}

bool Matrix_f32::NearlyEquals(const Matrix_f32 &mat2, double tolerance)
{
    return AreNearlyEquall(*this, mat2, tolerance);
}

bool Matrix_f32::IsSymmetric(float tolerance)
{
    return IsSymmetric(*this, tolerance);
}

Matrix_f32 Matrix_f32::Identity(_INDEX dimension)
{
    Matrix_f32 uMatrix(dimension, dimension);
	
	for (size_t i = 0; i < dimension; i++)
		uMatrix.SetValue(i, i, 1.0f);

	return uMatrix;
}

bool Matrix_f32::AreMultipliable(const Matrix_f32 &mat1, const Array2D &mat2)
{
    if (mat1.Columns() != mat2.Rows())
		return false;
	else
		return true;
}

bool Matrix_f32::IsInvertible(Matrix_f32 mat, bool checkSingular)
{
    if (mat.Rows() != mat.Columns())
		return false;
	else if (checkSingular && fabs(mat.Determinant()) <= MINDET) //Making use of the optimization where determinant computation wil be skipped if checkSingular is false
																//TODO check whether this is a standard and always the case, or is compiler/language version specific.
		return false;
	else
		return true;
}

bool Matrix_f32::IsSymmetric(const Matrix_f32 & mat, double tolerance)
{
    if (!IsSquared(mat))
		return false;
	
	tolerance = fabs(tolerance); //in case it was sent as a negative value, in-which case all matrices would turn out to be symmetric.

	for (_INDEX i = 0; i < mat.Rows(); i++)
	{
		for (_INDEX j = 0; j < mat.Rows(); j++)
		{
			if (fabs(mat.GetValue(i, j) - mat.GetValue(j, i)) >= tolerance)
				return false;
		}
	}

	return true;
}

bool Matrix_f32::AreNearlyEquall(const Matrix_f32 &mat1, const Array2D &mat2, double tolerance)
{
    tolerance = fabs(tolerance); //in case it was sent as a negative value, in-which case all matrices would turn out to be equall.
	if (!AreOfSameSize(mat1, mat2))
		return false;

	for (_INDEX i = 0; i < mat1.Rows(); i++)
	{
		for (_INDEX j = 0; j < mat1.Rows(); j++)
		{
			if (fabs(mat1.GetValue(i, j) - mat2.GetValue(i, j)) >= tolerance)
				return false;
		}
	}

	return true;
}

Matrix_f32 ** Matrix_f32::DecomposeLU(const Matrix_f32 & mat)
{
    //Based on the algorithm detailed on Introduction to Algorithms (3rd ed), Cormen, T., Lieserson, C., Rivest, R., and Stein, C.

	if (!IsSquared(mat))
	{
		std::cout << "ERROR! Cannot decompose a non-squared matrix." << std::endl;
		return NULL;
	}

	Matrix_f32 ** decomposition = new Matrix_f32 * [2]; //first is lower, second is upper.

	Matrix_f32 * lower = decomposition[0] = new Matrix_f32(mat.Rows(), mat.Rows());
	Matrix_f32 * upper = decomposition[1] = new Matrix_f32(mat.Rows(), mat.Rows());
	
	//we need a temporary matrix initially holding the same content of original matrix for the computations bellow
	Matrix_f32 * tempMatrix = new Matrix_f32(mat);

	for (int i = 0; i < mat.Rows(); i++)
	{
		//upper diagonal is same for original matrix
		upper->SetValue(i, i, tempMatrix->GetValue(i, i));

		//lower diagonal = 1
		lower->SetValue(i, i, 1.0f);

		for (int j = i + 1; j < mat.Rows(); j++)
		{
			//lower is divided by upper's pivot
			lower->SetValue(j, i, tempMatrix->GetValue(j, i) / upper->GetValue(i, i));

			//upper triangle (above diagonal) is the Schur compliment
			upper->SetValue(i, j, tempMatrix->GetValue(i, j));
		}

		//compute the Schur complement in-place in the temporary matrix
		for (int k = i + 1; k < mat.Rows(); k++)
		{
			for (int l = i + 1; l < mat.Rows(); l++)
			{
				tempMatrix->SetValue(k, l, tempMatrix->GetValue(k, l) - lower->GetValue(k, i) * upper->GetValue(i, l));
			}
		}
	}

	//cleanup and return
	delete tempMatrix;

	return decomposition;
}

Matrix_f32 ** Matrix_f32::DecomposeLUP(const Matrix_f32 &mat) //TODO finish implementing this
{
    //Based on the algorithm detailed on Introduction to Algorithms (3rd ed), Cormen, T., Lieserson, C., Rivest, R., and Stein, C.

	if (!IsSquared(mat))
	{
		std::cout << "ERROR! Cannot decompose a non-squared matrix." << std::endl;
		return NULL;
	}

	Matrix_f32 ** decomposition = new Matrix_f32 * [3]; //first is lower, second is upper, third is permutation

	return decomposition;
}

// Matrix_f32 Matrix_f32::VectorizedAddition(const Matrix_f32 & mat1, const Matrix_f32 & mat2)
// {
//     if (!AreOfSameSize(mat1, mat2))
// 	{
// 		std::cout << "ERROR! Attempting to add arrays of different sizes." << std::endl;
// 		return Matrix_f32();
// 	}

// 	Matrix_f32 result(mat1.Rows(), mat1.Columns());

// 	// _INDEX size = result.Rows() * result.Columns();
// 	// int remainder = size % _VECTOR_SIZE;
// 	// int vectorizedSpan = size - remainder;
// 	// float * a = &mat1.content[0][0];
// 	// float * b = &mat2.content[0][0];
// 	// float * c = &result.content[0][0];

// 	// for (int i = 0; i < vectorizedSpan; i += _VECTOR_SIZE)
// 	// {
// 	// 	__m256 result = _mm256_add_ps(_mm256_loadu_ps(a), _mm256_loadu_ps(b));
// 	// 	_mm256_storeu_ps(c, result);
// 	// 	a += _VECTOR_SIZE;// * sizeof(double);
// 	// 	b += _VECTOR_SIZE;// * sizeof(double);
// 	// 	c += _VECTOR_SIZE;// * sizeof(double);
// 	// }
	
// 	// //handle remainders, if exist
// 	// for (int i = vectorizedSpan; i < size; i++)
// 	// {
// 	// 	*c = (*a) + (*b);
// 	// 	a++;
// 	// 	b++;
// 	// 	c++;
// 	// }

// 	return result;
// }

// Matrix_f32 Matrix_f32::VectorizedMultiplication(const Matrix_f32 &mat1, const Matrix_f32 & mat2)
// {
    
// }

template <typename T>
void Matrix_f32::CopyFromArray2D(Array2D<T> sourceArr)
{
    DeleteContent();

	rows = sourceArr.Rows();
	columns = sourceArr.Columns();

	if (sourceArr.IsEmpty())
	{
		content = NULL;
		return;
	}

	Alloc(rows, columns);

	for (_INDEX i = 0; i < rows; i++)
		for (_INDEX j = 0; j < columns; j++)
			content[i][j] = static_cast<float>(sourceArr.GetValue(i, j));
}

Matrix_f32 Matrix_f32::AddMAtrices(const Matrix_f32 & mat1, const Matrix_f32 & mat2)
{
    if (!AreOfSameSize(mat1, mat2))
	{
		std::cout << "ERROR! Attempting to add arrays of different sizes." << std::endl;
		return Matrix_f32();
	}
	
	Matrix_f32 result(mat1.Rows(), mat1.Columns());

	for (_INDEX i = 0; i < mat1.Rows(); i++)
	{
		for (_INDEX j = 0; j < mat1.Columns(); j++)
		{
			double value = mat1.GetValue(i, j) + mat2.GetValue(i, j);
			result.SetValue(i, j, value);
		}
	}
	
	return result;
}

Matrix_f32 Matrix_f32::SubtractMatrices(const Matrix_f32 & mat1, const Matrix_f32 & mat2)
{
    if (!AreOfSameSize(mat1, mat2))
	{
		std::cout << "ERROR! Attempting to add arrays of different sizes." << std::endl;
		return Matrix_f32();
	}
	
	Matrix_f32 result(mat1.Rows(), mat1.Columns());

	for (_INDEX i = 0; i < mat1.Rows(); i++)
	{
		for (_INDEX j = 0; j < mat1.Columns(); j++)
		{
			double value = mat1.GetValue(i, j) - mat2.GetValue(i, j);
			result.SetValue(i, j, value);
		}
	}
	
	return result;
}

Matrix_f32 Matrix_f32::MultiplyMatrices(const Matrix_f32 & mat1, const Matrix_f32 & mat2)
{
    if (!AreMultipliable(mat1, mat2))
	{
		std::cout << "ERROR! Attempting to multiply array of " << mat1.Columns() << "columns with an array of " << mat2.Rows() << " rows." << std::endl;
		return Matrix_f32();  //really need to figure out how to make this thing more gracefull.
	}
	
	Matrix_f32 result(mat1.Rows(), mat2.Columns());

	for (_INDEX i = 0; i < mat1.Rows(); i++)
	{
		for (_INDEX j = 0; j < mat2.Columns(); j++)
		{
			float cellValue = 0.0f;

			for (_INDEX k = 0; k < mat1.Columns(); k++)
				cellValue += mat1.GetValue(i, k) * mat2.GetValue(k, j);

			result.SetValue(i, j, cellValue);
		}
	}

	return result;
}

Matrix_f32 Matrix_f32::MultiplayMatrixWithScalar(const Matrix_f32 & mat1, const double scalar)
{
    Matrix_f32 result = mat1;

	for (_INDEX i = 0; i < result.Rows(); i++)
	{
		for (_INDEX j = 0; j < result.Columns(); j++)
		{
			result.SetValue(i, j, result.GetValue(i, j) * scalar);
		}
	}

	return result;
}

Matrix_f32 Matrix_f32::AddMatricesVectorized(const Matrix_f32 & mat1, const Matrix_f32 & mat2)
{
    if (!AreOfSameSize(mat1, mat2))
	{
		std::cout << "ERROR! Attempting to add arrays of different sizes." << std::endl;
		return Matrix_f32();
	}

	Matrix_f32 result(mat1.Rows(), mat1.Columns());

	size_t size = result.Rows() * result.Columns();
	
	int remainder = size % _VECTOR_SIZE;
	int vectorizedSpan = size - remainder;
	float * a = &mat1.content[0][0];
	float * b = &mat2.content[0][0];
	float * c = &result.content[0][0];

	for (int i = 0; i < vectorizedSpan; i += _VECTOR_SIZE)
	{
        #ifdef _USE_AVX512
            AddVectors512_f32(a, b, c);
        #elif defined(_USE_AVX256)
            AddVectors256_f32(a, b, c);
        #elif defined(_USE_SSE)
            AddVectors128_f32(a, b, c);
        #endif

		a += _VECTOR_SIZE;// * sizeof(double);
		b += _VECTOR_SIZE;// * sizeof(double);
		c += _VECTOR_SIZE;// * sizeof(double);
	}
	
	//handle remainders, if exist
	for (int i = vectorizedSpan; i < size; i++)
	{
		*c = (*a) + (*b);
		a++;
		b++;
		c++;
	}

	return result;
}

Matrix_f32 Matrix_f32::SubtractMatricesVectorized(const Matrix_f32 & mat1, const Matrix_f32 & mat2)
{
    if (!AreOfSameSize(mat1, mat2))
	{
		std::cout << "ERROR! Attempting to add arrays of different sizes." << std::endl;
		return Matrix_f32();
	}

	Matrix_f32 result(mat1.Rows(), mat1.Columns());

	size_t size = result.Rows() * result.Columns();
	
	int remainder = size % _VECTOR_SIZE;
	int vectorizedSpan = size - remainder;
	float * a = &mat1.content[0][0];
	float * b = &mat2.content[0][0];
	float * c = &result.content[0][0];

	for (int i = 0; i < vectorizedSpan; i += _VECTOR_SIZE)
	{
        #ifdef _USE_AVX512
            SubtractVectors512_f32(a, b, c);
        #elif defined(_USE_AVX256)
            SubtractVectors256_f32(a, b, c);
        #elif defined(_USE_SSE)
            SubtractVectors128_f32(a, b, c);
        #endif

		a += _VECTOR_SIZE;// * sizeof(double);
		b += _VECTOR_SIZE;// * sizeof(double);
		c += _VECTOR_SIZE;// * sizeof(double);
	}
	
	//handle remainders, if exist
	for (int i = vectorizedSpan; i < size; i++)
	{
		*c = (*a) + (*b);
		a++;
		b++;
		c++;
	}

	return result;
}

Matrix_f32 Matrix_f32::MultiplyMatricesVectorized(const Matrix_f32 & mat1, const Matrix_f32 & mat2)
{
    Matrix_f32 result(mat1.Rows(), mat1.Columns());
	

	return result;
}

Matrix_f32 Matrix_f32::MultiplayMatrixWithScalarVectorized(const Matrix_f32 & mat1, const double scalar)
{
    Matrix_f32 result = mat1;
	

	return result;
}

Matrix_f32 Matrix_f32::InvertMatrix(const Matrix_f32 & sourceMat, MatrixInversionMethod method)
{
    if (!IsInvertible(sourceMat))
	{
		std::cout << "ERROR! Attempting to invert a non-invertible array." << std::endl;
		return Matrix_f32();
	}
	
	switch (method)
	{
	case MatrixInversionMethod::Gauss_Jordan:
		return GausJordanElimination(sourceMat);
	default:
		std::cout << "ERROR! Received unsupported MatrixInversionMethod in InverArray()." << std::endl; //shouldn't happen (but still...)
		return Matrix_f32();
	}
}

double Matrix_f32::CalculateDeterminant(const Matrix_f32 & mat)
{
    //This recurive method calculate the determinant for a matrix of arbitrary dimensions. If the recieved matrix is less than 2x2, directly return the determinant, else will make use of
	//the method GetMinorSubMatrix() and work recuresively untill reaching the 2x2 matrices. 
	double result = 0.0f;

	if (!IsSquared(mat))
	{
		std::cout << "ERROR! Attempting to calculate the determinant of a non-squared matrix." << std::endl;
		return result;  //really need to figure out how to make this thing more gracefull.
	}
	
	if (mat.Rows() > 2)
	{
		for (size_t i = 0; i < mat.Rows(); i++)
		{
			result += pow(-1.0f, i) * mat.GetValue(0, i) * CalculateDeterminant(GetMinorSubMatrix(mat, 0, i));	//this is where the recurssion happens. The pow(-1.0f, i) term is used
		}																													//to flip the sign of the addition based on which column we're at.
	}
	else
	{
		result = (mat.GetValue(0, 0) * mat.GetValue(1, 1)) - (mat.GetValue(1, 0) * mat.GetValue(0, 1));
	}

	return result;
}

Matrix_f32 Matrix_f32::GausJordanElimination(const Array2D & sourceMat)
{
    Matrix_f32 result(sourceMat.Rows(), sourceMat.Columns());
	Matrix_f32 augmentedArr = MergeArrays(sourceMat, Identity(sourceMat.Rows())); //augmentedArr is the augment matrix, which is the original matrix with a identity matrix attached to its right.

	//If first pivot value is zero, must swap the row with another that has a non-zero value. If none exist, can't use this method.
	if (augmentedArr.GetValue(0, 0) == 0.0f)
	{
		//look for a row where first element is not zero
		for (int i = 0; i < augmentedArr.Rows(); i++)
		{
			if (augmentedArr.GetValue(i, 0) != 0.0f)
			{
				augmentedArr.SwapRows(0, i);
				break;
			}
			if (i == augmentedArr.Rows())
			{
				std::cout << "ERROR! Could not find a suitable pivot value for first row. Can't invert useing Gauss-Jordan Elimmination" << std::endl;
				return Matrix_f32();
			}
		}
	}

	//For an array of n*n size, n steps are needed to get the inverse.
	for (size_t step = 0; step < result.Rows(); step++) //the variable "step" here will be our pivot row, and by virtue of being a squared array, our pivot column as well.
	{
		double pivotValue = augmentedArr.GetValue(step, step);  //the value of th pivot cell is always on the diagonal.
		
		//mFactor is the value that, when multiplied with our pivot row and substracted from current row, should help reduce it towards zero (except diagonal value)
		//we extract mFactors before the loops because, their cell values will change inside the loops at first iterations, but we'll still be needing them for remaining iterations.
		double * mFactors = new double[result.Rows()]; 
		for (size_t i = 0; i < augmentedArr.Rows(); i++)
			mFactors[i] = augmentedArr.GetValue(i, step);

		for (size_t j = 0; j < augmentedArr.Columns(); j++)
		{
			double newColumnValueAtPivotRow = augmentedArr.GetValue(step, j) / pivotValue;
			augmentedArr.SetValue(step, j, newColumnValueAtPivotRow); //the pivot row value is adjusted before iterating over other rows.

			for (size_t i = 0; i < augmentedArr.Rows(); i++) //we then iterate over other rows, using the if-statement bellow to avoid modifying our pivot row.
			{
				if (i != step)
				{
					double newValueAtOtherRow = augmentedArr.GetValue(i, j) - (mFactors[i] * newColumnValueAtPivotRow);
					augmentedArr.SetValue(i, j, newValueAtOtherRow);
				}
				//augmentedArr.DisplayArrayInCLI();
			}
		}
		delete[] mFactors;
	}

	//extract result from augmentedArr
	result = augmentedArr.GetSubMatrix(0, result.Rows(), sourceMat.Columns(), result.Columns());

	return result;
}


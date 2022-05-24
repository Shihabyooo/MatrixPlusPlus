#pragma once
#include "ArrayPP.hpp"
#include "VectorMaths.hpp"

//instaniate necessary versions of the Array2D template.
template class Array2D<float>;
template class Array2D<int>;
template class Array2D<long>;

class Matrix_f32 : public Array2D<float> {
public:

    using Array2D<float>::Array2D; //To inherit base class constructors.

    Matrix_f32(Array2D<float> sourceArr);
    Matrix_f32(Array2D<int> sourceArr);
    Matrix_f32(Array2D<long> sourceArr);

    Matrix_f32 operator* (const Matrix_f32 & mat2);
    Matrix_f32 operator* (const double & scalar);
    Matrix_f32 operator+ (const Matrix_f32 & mat2);
	Matrix_f32 operator- (const Matrix_f32 & mat2);

    //todo add a special GetSubMAtrix class to return Matrix_f32.
    Matrix_f32 GetSubMatrix(_INDEX beginRow, _INDEX noOfRows, _INDEX beginColumn, _INDEX noOfColumns);
    Matrix_f32 Invert();
    double Determinant();
    void Overlay(const Matrix_f32 mat2, _INDEX rowOffset, _INDEX columnOffset); //Add another Array2D of non-equal size to this Array2D element by element. If the second Array2D is larger, elements outside the boundary will be clipped. rowOffset and columnOffset determine which elements of the first Array2D the first element of the second Array2D will be added to. //TODO vectorize this
    Matrix_f32 ** DecomposeLU(); //For squared matrices only. Decomposes the matrix into upper and lower triangular matrices, returns as an array of two [pointers to] Array2Ds, the first being the lower decomposition, and the second the upper. Returns NULL if not decomposable.
    Matrix_f32 ** DecomposeLUP(); //For squared matrices only. Decomposes the matrix into upper and lower triangular matrices with permuation, returns as an array of three [pointers to] Array2Ds, the first being the lower decomposition, and the second the upper, the third is the permuation matrix. Returns NULL if not decomposable.
    bool NearlyEquals(const Matrix_f32 &mat2, double tolerance);
    bool IsSymmetric(float tolerance);

    //Static methods
    static Matrix_f32 Identity(_INDEX dimension);	//simple function for quick construction of a identity matrix of size: dimension*dimension.
    static bool AreMultipliable(const Matrix_f32 &mat1, const Array2D &mat2);	//for m1*n1 and m2*n2 matrices, tests that n1 == m2.
    static bool IsInvertible(Matrix_f32 mat, bool checkSingular = false);	//Computing the determinant can take ages with current algorithm for large matrices, so the singularity check is optional.
    static bool IsSymmetric(const Matrix_f32 &mat, double tolerance); //For squared matrices only. Returns false for non-squared matrices. tolerance is the (absolute) allowed error in difference between counterpart values, bellow which they are considered equal.
    static bool AreNearlyEquall(const Matrix_f32 &mat1, const Array2D &mat2, double tolerance);
    static Matrix_f32 ** DecomposeLU(const Matrix_f32 &mat);  //For squared matrices only. Decomposes the matrix into upper and lower triangular matrices, returns as an array of two [pointers to] Array2Ds, the first being the lower decomposition, and the second the upper. Returns NULL if not decomposable. Deleting memory is responcibility of caller.
	static Matrix_f32 ** DecomposeLUP(const Matrix_f32 &mat); //For squared matrices only. Decomposes the matrix into upper and lower triangular matrices with permuation, returns as an array of three [pointers to] Array2Ds, the first being the lower decomposition, and the second the upper, the third is the permuation matrix. Returns NULL if not decomposable. Deleting memory is responcibility of caller.
	//static Matrix_f32 VectorizedAddition(const Matrix_f32 &mat1, const Matrix_f32 &mat2); //test
	//static Matrix_f32 VectorizedMultiplication(const Matrix_f32 &mat1, const Matrix_f32 &mat2); //test

protected:
    template <typename T>
    void CopyFromArray2D(Array2D<T> sourceArr);

    Matrix_f32 AddMAtrices(const Matrix_f32 &mat1, const Matrix_f32 &mat2);
    Matrix_f32 SubtractMatrices(const Matrix_f32 &mat1, const Matrix_f32 &mat2);
    Matrix_f32 MultiplyMatrices(const Matrix_f32 & mat1, const Matrix_f32 & mat2);
	Matrix_f32 MultiplayMatrixWithScalar(const Matrix_f32 &mat1, const double scalar);    

	Matrix_f32 AddMatricesVectorized(const Matrix_f32 &mat1, const Matrix_f32 &mat2);
    Matrix_f32 SubtractMatricesVectorized(const Matrix_f32 &mat1, const Matrix_f32 &mat2);
    Matrix_f32 MultiplyMatricesVectorized(const Matrix_f32 & mat1, const Matrix_f32 & mat2);
    Matrix_f32 MultiplayMatrixWithScalarVectorized(const Matrix_f32 &mat1, const double scalar);

	Matrix_f32 InvertMatrix(const Matrix_f32 & sourceMat, MatrixInversionMethod method = MatrixInversionMethod::Gauss_Jordan);	//switch-statement based on method to use appropriate implementation. Now only Gauss_Jordan, more in future.

    double CalculateDeterminant(const Matrix_f32 & mat);

    //matrix Inversion Methods
	Matrix_f32 GausJordanElimination(const Array2D & sourceMat);
};
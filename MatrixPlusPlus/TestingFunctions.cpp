#include <iostream>
//#include <math.h>
//#include <vector>
#include <chrono>

//#include "Array2D.hpp"
#include "MatricesPP.hpp"
#include "VectorMaths.hpp"

#define NOW std::chrono::high_resolution_clock::now()
#define TIME std::chrono::high_resolution_clock::time_point
#define DURATION(x) std::chrono::duration_cast<std::chrono::microseconds>(NOW - x).count()

int main(int argc, char ** argv)
{	
	//int size1 = 1000;
	//int size2 = 1000;
	//int trialsCount = 10;

	// for (int  trial = 0; trial < trialsCount; trial++)
	// {
	// 	Matrix_f32 matA = Matrix_f32(size1, size2);
	// 	Matrix_f32 matB = Matrix_f32(size2, size1);


	// 	srand(time(0));
	// 	for (int i = 0; i < size1; i++)
	// 	{
	// 		for (int j = 0; j < size2; j++)
	// 		{
	// 			matA[i][j] =  rand()%20;
	// 			matB[j][i] =  rand()%20;
	// 		}
	// 	}

	// 	std::chrono::high_resolution_clock::time_point start = NOW;
	// 	std::cout << "Mult1: " << std::endl;
	// 	Matrix_f32 multResultVector =  Matrix_f32::MultiplyMatricesVectorized_N(matA, matB);
	// 	long vecTime = std::chrono::duration_cast<std::chrono::microseconds>(NOW - start).count();

	// 	std::cout << "Mult2: " << std::endl;
	// 	start = NOW;
	// 	Matrix_f32 multResult_Normal = Matrix_f32::MultiplyMatrices(matA, matB);
	// 	long normTime = std::chrono::duration_cast<std::chrono::microseconds>(NOW - start).count();
		
	// 	//multResultVector.DisplayOnCLI();
	// 	//multResult_Normal.DisplayOnCLI();
	// 	//(multResultVector - multResult_Normal).DisplayOnCLI();

	// 	std::cout << "Normal time: " << normTime << std::endl;
	// 	std::cout << "Vector Time" << vecTime << std::endl;
	// 	std::cout << "gain: " << 100.0f * (normTime - vecTime) / vecTime << "%" << std::endl;

	// }
	// std::cout << "Press any key to continue." << std::endl;
	// std::cin.sync();
	// std::cin.get();

	// Matrix_f32 oneByOne = Matrix_f32(1, 1, 5.0f);
	// Matrix_f32 twoByTwo = Matrix_f32(2, 2);
	// Matrix_f32 threeByThree = Matrix_f32 (3, 3);

	// srand(time(0));
	// for (size_t i = 0; i < 2; i++)	
	// 	for (size_t j = 0; j < 2; j++)		
	// 		twoByTwo[i][j] = rand()%20;

	// for (size_t i = 0; i < 3; i++)	
	// 	for (size_t j = 0; j < 3; j++)		
	// 		threeByThree[i][j] = rand()%20;

	// std::cout << "1 by 1:\n";
	// (oneByOne * Matrix_f32::Invert1x1Matrix(oneByOne)).DisplayOnCLI(10);
	
	// std::cout << "2 by 2:\n";
	// (twoByTwo * Matrix_f32::Invert2x2Matrix(twoByTwo)).DisplayOnCLI(10);
	
	// std::cout << "3 by 3:\n";
	// (threeByThree * Matrix_f32::Invert3x3Matrix(threeByThree)).DisplayOnCLI(10);
	
	//int size = 10;
	//Matrix_f32 invTest(size, size);

	//srand(time(0));
	//for (size_t i = 0; i < size; i++)	
 //		for (size_t j = 0; j < size; j++)		
	//		invTest[i][j] = rand()%20;

	//invTest.DisplayOnCLI();
	//Matrix_f32::InvertMatrix(invTest, MatrixInversionMethod::Blockwise).DisplayOnCLI();
	//(invTest * Matrix_f32::InvertMatrix(invTest, MatrixInversionMethod::Blockwise)).DisplayOnCLI(7);
	//std::cout << "Determinant of result: " << (invTest *Matrix_f32::InvertMatrix(invTest, MatrixInversionMethod::Blockwise)).Determinant() << std::endl;

	/*Vector_f32 testVec(5, 1.0f);
	Vector_f32 testVec2(5, 2.0f);
	
	testVec[1] = testVec[2] = 99.0f;
	testVec = testVec + testVec2;

	testVec.DisplayOnCLI();*/

	/*srand(time(0));
	Matrix_f32 testMat(5, 5);
	for (int i = 0; i < 5; i++)
		for (int j = 0; j < 5; j++)
			testMat[i][j] = rand() % 10;*/

	/*testMat.DisplayOnCLI();
	float ** column0 = testMat.GetColumnPtr(0);
	float ** column1 = testMat.GetColumnPtr(1);
	float ** column2 = testMat.GetColumnPtr(2);
	
	for (int i = 0; i < 5; i++)
		std::cout << *(column0[i]) << std::endl;
	std::cout << "\n";
	for (int i = 0; i < 5; i++)
		std::cout << *(column1[i]) << std::endl;
	std::cout << "\n";
	for (int i = 0; i < 5; i++)
		std::cout << *(column2[i]) << std::endl;
	std::cout << "\n";

	auto columnCopy1 = testMat.GetColumn(1);
	auto rowCopy1 = testMat.GetRow(1);
	
	for (int i = 0; i < 5; i++)
		std::cout << columnCopy1[i] << std::endl;
	std::cout << "\n";

	for (int i = 0; i < 5; i++)
		std::cout << rowCopy1[i] << "\t";
	std::cout << "\n";*/
	/*Array2D<float> arrayF(10, 10, 2.0f);
	Matrix_f32 matrix(10, 10, 2.0f);
	Vector_f32 vector(10, 4.0f);*/

	//arrayF.DisplayOnCLI();

	/*Array2D<float> array2 = std::move(arrayF);
	Matrix_f32 matrix2 = std::move(matrix);
	Vector_f32 vector2 = std::move(vector);*/

	/*arrayF.DisplayOnCLI();
	array2.DisplayOnCLI();
	array2 = std::move(arrayF);
	array2.DisplayOnCLI();*/

	//array2 = matrix2 * matrix2;
	//array2 = std::move(matrix2);
	//array2 = std::move(vector2);

	/*Vector_f32 vector3 = std::move(matrix2);
	vector3.DisplayOnCLI();
	vector2.DisplayOnCLI();
	float result = vector3.DotProduct(vector2);
	std::cout << result << "\n";*/
	

	return 0;
}
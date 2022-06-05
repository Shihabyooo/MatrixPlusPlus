#include <iostream>
//#include <math.h>
//#include <vector>
#include <chrono>

//#include "Array2D.hpp"
#include "MatricesPP.hpp"
#include "VectorMaths.hpp"


int main(int argc, char ** argv)
{	
	int size1 = 1000;
	int size2 = 1000;
	int trialsCount = 10;
	for (int  trial = 0; trial < trialsCount; trial++)
	{
		Matrix_f32 matA = Matrix_f32(size1, size2);
		Matrix_f32 matB = Matrix_f32(size2, size1);


		srand(time(0));
		for (int i = 0; i < size1; i++)
		{
			for (int j = 0; j < size2; j++)
			{
				
				matA[i][j] =  rand()%20;
				matB[j][i] =  rand()%20;
			}
		}

		#define NOW std::chrono::high_resolution_clock::now()

		std::chrono::high_resolution_clock::time_point start = NOW;
		std::cout << "Mult1: " << std::endl;
		Matrix_f32 multResultVector =  Matrix_f32::MultiplyMatricesVectorized_N(matA, matB);
		long vecTime = std::chrono::duration_cast<std::chrono::microseconds>(NOW - start).count();

		std::cout << "Mult2: " << std::endl;
		start = NOW;
		Matrix_f32 multResult_Normal = Matrix_f32::MultiplyMatrices(matA, matB);
		long normTime = std::chrono::duration_cast<std::chrono::microseconds>(NOW - start).count();
		
		//multResultVector.DisplayArrayInCLI();
		//multResult_Normal.DisplayArrayInCLI();
		//(multResultVector - multResult_Normal).DisplayArrayInCLI();

		std::cout << "Normal time: " << normTime << std::endl;
		std::cout << "Vector Time" << vecTime << std::endl;
		std::cout << "gain: " << 100.0f * (normTime - vecTime) / vecTime << "%" << std::endl;

	}
	// std::cout << "Press any key to continue." << std::endl;
	// std::cin.sync();
	// std::cin.get();

	return 0;
}
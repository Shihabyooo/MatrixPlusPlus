#include <iostream>
//#include <math.h>
//#include <vector>
#include <chrono>

//#include "Array2D.hpp"
#include "MatricesPP.hpp"


int main(int argc, char ** argv)
{
	std::cout << "TESTTTT!" << std::endl;
	std::cout<< "Create 3x3 int array with init value of 777" << std::endl;
	Array2D<int> intArray1(3, 3, 777);
	intArray1.DisplayArrayInCLI();

	intArray1[0][0] = 9;
	intArray1[1][1] = 6;
	intArray1[2][2] = 'c';
	std::cout << " val : " << intArray1[0][0] << ", " << intArray1[1][1] <<", " << intArray1[2][2] <<std::endl;

	int * haha;
	*intArray1[1] = 2;

	intArray1.DisplayArrayInCLI();

	std::cout<< "Create 3x3 char array with init value of x" << std::endl;
	Array2D<char> charArray(3, 3, 'X');
	charArray.DisplayArrayInCLI();

	std::cout<< "Create 3x3 string array with init value of xxx" << std::endl;
	Array2D<std::string> stringArray(3, 3, "xxx");
	stringArray.DisplayArrayInCLI();

	stringArray(1,1) = 4;
	std::cout << " val : " << stringArray[1][1] << std::endl;


	std::cout << "Press any key to continue." << std::endl;
	std::cin.sync();
	std::cin.get();

	return 0;
}
#include "VectorPP_f32.hpp"

Vector_f32::Vector_f32()
{
	rows = 0;
	columns = 1;
	content = NULL;
}

Vector_f32::Vector_f32(_INDEX const size)
{
	Alloc(size, 1);

	columns = 1;
	rows = size;
}

Vector_f32::Vector_f32(_INDEX const size, float defaultValue)
{
	Alloc(size, 1);

	columns = 1;
	rows = size;
	SetEntireArrayToFixedValue(defaultValue);
}

Vector_f32::Vector_f32(Vector_f32 const & sourceVec)
{
	*this = sourceVec;
}

Vector_f32::Vector_f32(Matrix_f32 const & sourceMat, _INDEX column)
{
	VectorFromMatrix(sourceMat, column);
}

float & Vector_f32::operator[](const _INDEX row)
{
#ifdef _USE_BOUNDS_CHECK
	if (content == NULL			//Checking whether this object is empty. Making an assumption that initializing the first level of the content is automatically followed by init of sublevel.
		|| _row >= rows)		//Checking out-of-bound writes.
		throw std::out_of_range("ERROR! Row value out of range or content set to NULL");
#endif

	return content[row][0];
}

Vector_f32::~Vector_f32()
{
	DeleteContent();
}

Vector_f32 Vector_f32::operator+(Vector_f32 const & vec2) const
{
#ifdef _VECTORIZED_CODE
	return AddVectorsVectorized(*this, vec2);
#else
	return AddVectors(*this, vec2);
#endif
}

Vector_f32 Vector_f32::operator-(Vector_f32 const & vec2) const
{
#ifdef _VECTORIZED_CODE
	return SubtractVectorsVectorized(*this, vec2);
#else
	return SubtractVectors(*this, vec2);
#endif
}

double Vector_f32::Magnitude()
{
	double mag = 0;
	for (_INDEX i = 0; i < rows; i++)
		mag += pow(content[i][0], 2);
	return sqrt(mag);
}

void Vector_f32::VectorFromMatrix(Matrix_f32 const & sourceMat, _INDEX column)
{
	DeleteContent();
	rows = sourceMat.Rows();
	columns = 1;
	Alloc(rows, 1);

	for (int i = 0; i < rows; i++)
		content[i][0] = sourceMat.GetValue(i, column);
}

Vector_f32 Vector_f32::AddVectors(Vector_f32 const & vec1, Vector_f32 const & vec2)
{
#ifdef _USE_BOUNDS_CHECK
	if (!AreOfSameSize(vec1, vec1))
	{
		std::cout << "ERROR! Attempting to add vectors of different sizes." << std::endl;
		return Vector_f32();
	}
#endif

	Vector_f32 result(vec1.Rows(), 1);

	float * a = &vec1.content[0][0];
	float * b = &vec2.content[0][0];
	float * c = &result.content[0][0];

	for (size_t i = 0; i < result.rows; i++)
	{
		*c = *a + *b;
		a++;
		b++;
		c++;
	}

	return result;
}

Vector_f32 Vector_f32::SubtractVectors(Vector_f32 const & vec1, Vector_f32 const & vec2)
{
#ifdef _USE_BOUNDS_CHECK
	if (!AreOfSameSize(vec1, vec1))
	{
		std::cout << "ERROR! Attempting to add vectors of different sizes." << std::endl;
		return Vector_f32();
	}
#endif

	Vector_f32 result(vec1.Rows(), 1);

	float * a = &vec1.content[0][0];
	float * b = &vec2.content[0][0];
	float * c = &result.content[0][0];

	for (size_t i = 0; i < result.rows; i++)
	{
		*c = *a - *b;
		a++;
		b++;
		c++;
	}

	return result;
}


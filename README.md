#MatricesPlusPLus
Rework of Array2D project.

A set of C++ classes (.hpp + .cpps files) defining Matrix types to work-around C++'s weakness with dynamic matrices.
Comes with some maths operators overloaded for basic matrix operations plus methods for some advanced operations.

This project initially started as a quick class to use in other projetcs. Later it evolved into a project primarily for learning purposes.

Four matrix types are implemented (or planned to be):
Matrix_f32 - stores 32bit floating point values (floats), and supports all matrix mathematica ops (that were implemented).
Matrix_f64 - stores 64bit floating point values (doubles), ditto.
Matrix_i32 - stores 32 bit integer values (ints), supports a subset of math ops.
Matrix_i64 - stores 64 bit integer values (longs), ditto.

Additionally, an Array2D class is exposed that allows storing any type. Naturally, it doesn't allow for any mathematical ops.

Typically you would only need to include MatricesPP.hpp to access all those types.

Vectorized Arithematics:
This code tries to make use of the vectorization capablity of modern processors (i.e. SIMD instructions). Supported capablity depends on the processor. Most processors these days support SSE and AVX instruction sets, many AVX2 and FMA.

The actual vectorized capability the code will *try* to compile to are controlled using flags defined in VectorMaths.cpp (namely: "USE_SSE," "_USE_AVX256," "_USE_AVX512" and "_USE_FMA256"). The code will use the widest of the first three (SSE, AVX256, AVX512) set. i.e. if all were set, the code will try to use AVX512. 

Note that the compiler flags must be set independantly to support the selected vectorization capability.
128 bit functions require SSE/SSE2, 256bit floating point functions require AVX, 256bit integer function require AVX2, 512bit functions require AVX-512.
Also note that, in the current form, all 256bit operations require AVX2 support.

See the following link for more details:
https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html

The included Visual Studio Code tasks.json is configured for clang-12 with AVX, AVX2 and FMA support. SSE/SSE2 are implied when targeting x32-64

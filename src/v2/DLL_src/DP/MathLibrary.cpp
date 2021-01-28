// #include "pch.h"
#include "MathLibrary.h"


double add_c(int a, int b)
{
	double c = 0.0;

	c += a * b;

	for (int i = 0; i < 100000000; i++)
	{
		c += 0.001;
	}
	return c;
}

// #include "main_api.cpp"
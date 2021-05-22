#include<iostream>
#include"omp.h"

using namespace std;

void main()
{
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < 16; i++) 	{
		printf("i=%d, thread_id=%d\n", i, omp_get_thread_num());
	}
}
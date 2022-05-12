#ifndef KERNEL_H
#define KERNEL_H
#pragma once
#include<math.h>
#include <fstream>
const double M_PI = 3.14159265358979323846;
void calc_kernel_high_pass(int filter_length, double cut_off, double* fltr_kernel_arr_, bool );
#endif



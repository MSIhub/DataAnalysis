#pragma once
#include <iostream>
#include <math.h>

void lp_windowed_sinc_fltr(double* sig_src_arr, double* sig_dest_arr, double* fltr_kernel_dest_arr, double cutoff_freq, int filter_length, int input_sig_length);

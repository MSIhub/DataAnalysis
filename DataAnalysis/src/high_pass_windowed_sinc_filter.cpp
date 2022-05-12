#include "high_pass_windowed_sinc_filter.h"

void hp_windowed_sinc_fltr(
	double* fltr_kernel_dest_arr,
	double cutoff,
	int filter_length,
	double* sig_src_arr,
	double* sig_dest_arr,
	int input_sig_length)
{
	const double M_PI = 3.14159265358979323846;
	double threshold = 0;
	//Calculate the first low-pass filter kernel 
	for (int i = 0; i < filter_length; i++)
	{
		threshold = (i - (double)(filter_length / 2));
		if (threshold == 0)
		{
			fltr_kernel_dest_arr[i] = 2 * M_PI * cutoff;
		}
		if (threshold != 0)
		{
			fltr_kernel_dest_arr[i] = sin(2 * M_PI * cutoff * threshold) / threshold;
			fltr_kernel_dest_arr[i] *= (0.42 - 0.5 * cos(2 * M_PI * i / filter_length) + 0.08 * cos(4 * M_PI * i / filter_length));
		}
	}

	//Change low-pass filter to a high-pass filter using spectral inversion
	for (int i = 0; i < filter_length; i++)
	{
		fltr_kernel_dest_arr[i] = -fltr_kernel_dest_arr[i];
	}
	fltr_kernel_dest_arr[filter_length / 2] += 1;

	//Convolve the input signal and filter kernel 
	for (int j = filter_length; j < input_sig_length; j++)
	{
		sig_dest_arr[j] = 0;
		for (int i = 0; i < filter_length; i++)
		{
			sig_dest_arr[j] = sig_dest_arr[j] + sig_src_arr[j - i] * fltr_kernel_dest_arr[i];

		}

	}
}
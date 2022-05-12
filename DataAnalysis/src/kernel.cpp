#include "kernel.h"

void calc_kernel_high_pass(int filter_length, double cut_off, double* fltr_kernel_arr_, bool isLogging)
{
	double threshold = 0.0;
	//Calculate the first low-pass filter kernel 
	for (int i = 0; i < filter_length; i++)
	{
		threshold = (i - (double)(filter_length / 2));
		if (threshold == 0)
		{
			fltr_kernel_arr_[i] = 2 * M_PI * cut_off;
		}
		if (threshold != 0)
		{
			fltr_kernel_arr_[i] = sin(2 * M_PI * cut_off * threshold) / threshold;
			fltr_kernel_arr_[i] *= (0.42 - 0.5 * cos(2 * M_PI * i / filter_length) + 0.08 * cos(4 * M_PI * i / filter_length));
		}
	}

	//Change low-pass filter to a high-pass filter using spectral inversion
	for (int i = 0; i < filter_length; i++)
	{
		fltr_kernel_arr_[i] = -fltr_kernel_arr_[i];
	}
	fltr_kernel_arr_[filter_length / 2] += 1;

	if (isLogging)
	{
		std::fstream output_sig_rt_kernel_fptr;
		output_sig_rt_kernel_fptr.open("log/rt_ay/output_sig_rt_kernel.dat" , std::fstream::in | std::fstream::out | std::fstream::app);

		if (output_sig_rt_kernel_fptr.is_open())
		{
			for (int i = 0; i < filter_length; i++)
			{
				output_sig_rt_kernel_fptr << "\n";
				output_sig_rt_kernel_fptr << fltr_kernel_arr_[i];//column vector
			}
			output_sig_rt_kernel_fptr.close();
		}

	}

}
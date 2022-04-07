#include "band_pass_windowed_sinc_filter.h"
//10kHz-------> 0.2 : (10 / 24) * 0.5 = 0.2
void bp_windowed_sinc_fltr(
	double * lower_cutoff_state_buff,
	double * upper_cutoff_state_buff,
	double * fltr_kernel_dest_arr,
	double lower_cutoff,
	double upper_cutoff,
	int filter_length,
	double * sig_src_arr,
	double * sig_dest_arr,
	int input_sig_length)
{
	double M_PI = 3.14159265358979323846;
	double threshold = 0;
	//Calculate the first low-pass filter kernel 
	for (int i = 0; i < filter_length; i++)
	{
		threshold = (i - (double)(filter_length / 2));
		if (threshold == 0)
		{
			lower_cutoff_state_buff[i] = 2 * M_PI * lower_cutoff;
		}
		if (threshold != 0)
		{
			lower_cutoff_state_buff[i] = sin(2* M_PI * lower_cutoff * threshold)/ threshold;
			lower_cutoff_state_buff[i] *= (0.42 - 0.5 * cos(2 * M_PI * i / filter_length) + 0.08 * cos(4 * M_PI * i / filter_length));
		}
	}

	//Calculate the second low-pass filter kernel
	for (int i = 0; i < filter_length; i++)
	{
		threshold = (i - (double)(filter_length / 2));
		if (threshold == 0)
		{
			upper_cutoff_state_buff[i] = 2 * M_PI * upper_cutoff;
		}
		if (threshold != 0)
		{
			upper_cutoff_state_buff[i] = sin(2* M_PI * upper_cutoff * threshold)/ threshold;
			upper_cutoff_state_buff[i] *= (0.42 - 0.5 * cos(2 * M_PI * i / filter_length) + 0.08 * cos(4 * M_PI * i / filter_length));
		}
	}

	//Change low-pass filter to a high-pass filter using spectral inversion
	for (int i = 0; i < filter_length; i++)
	{
		upper_cutoff_state_buff[i] = -upper_cutoff_state_buff[i];
	}
	upper_cutoff_state_buff[filter_length / 2] += 1;


	//Add low-pass filter to high-pass filter kernel to form a band-reject filter kernel
	for (int i = 0; i < filter_length; i++)
	{
		fltr_kernel_dest_arr[i] = lower_cutoff_state_buff[i] + upper_cutoff_state_buff[i];
	}

	//Change band reject filter into band pass filter using spectral inversion
	for (int i = 0; i < filter_length; i++)
	{
		fltr_kernel_dest_arr[i] = -(fltr_kernel_dest_arr[i]);
	}
	fltr_kernel_dest_arr[filter_length/2] = fltr_kernel_dest_arr[filter_length / 2] + 1;


	//Convolve the input signal and the filter kernel 
	for (int j = filter_length; j < input_sig_length; j++)
	{
		sig_dest_arr[j] = 0; 
		for (int i = 0; i < filter_length; i++)
		{
			sig_dest_arr[j] += sig_src_arr[j - i] * fltr_kernel_dest_arr[i];
		}
	}

}
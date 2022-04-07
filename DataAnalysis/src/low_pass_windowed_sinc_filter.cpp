#include "low_pass_windowed_sinc_filter.h"

/**@note Cutoff frequency is normalized and must be between 0 and 0.5, where 0.5
represents the nyquist frequency. 
In this example, the signal was sampled at 48kHz, therfore the nyquist frequency is 24kHz.
i.e. 24kHz ------> 0.5, 
-- 10kHz  -------> 0.2 : (10/24) *0.5 = 0.2 
*/
void lp_windowed_sinc_fltr(double* sig_src_arr, double* sig_dest_arr, double* fltr_kernel_dest_arr, double cutoff_freq, int filter_length, int input_sig_length)
{
	double M_PI = 3.14159265358979323846;
	double threshold = 0;
	//Calculate the low pass filter kernel 
	for (int i = 0; i < filter_length; i++)
	{
		threshold = (i - (double)(filter_length / 2));
		if (threshold == 0)
		{
			fltr_kernel_dest_arr[i] = 2 * M_PI * cutoff_freq;
		}
		if (threshold != 0)
		{
			fltr_kernel_dest_arr[i] = sin(2*M_PI*cutoff_freq * threshold)/ threshold;
			fltr_kernel_dest_arr[i] = fltr_kernel_dest_arr[i] * (0.54 - 0.46 * cos(2 * M_PI * i / filter_length));// Hamming window
		}
	}
	
	//Convolve the input signal and filter kernel 
	for (int j = filter_length; j < input_sig_length; j++)
	{
		sig_dest_arr[j] = 0;
		for (int i = 0; i < filter_length; i++)
		{
			sig_dest_arr[j] = sig_dest_arr[j] + sig_src_arr[j-i] * fltr_kernel_dest_arr[i];
			
		}
		
	}
}
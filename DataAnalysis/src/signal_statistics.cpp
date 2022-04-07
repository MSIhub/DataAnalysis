#include "signal_statistics.h"

double calc_signal_mean(double* sig_src_arr, int sig_length)
{
	double _mean = 0.0;
	for (int i = 0; i < sig_length; i++)
	{
		_mean += sig_src_arr[i];
	}
	_mean /= (double)sig_length;
	return _mean;
}

double calc_signal_variance(double* sig_src_arr, double sig_mean, int sig_length)
{
	double _variance = 0.0;
	for (int i = 0; i < sig_length; i++)
	{
		_variance += pow((sig_src_arr[i] - sig_mean), 2);
	}
	_variance /= (double)(sig_length - 1.0);
	return _variance;
}

double calc_signal_std(double sig_variance)
{
	return sqrt(sig_variance);
}
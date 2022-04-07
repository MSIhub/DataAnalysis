void convolution(double* sig_src_arr, double* sig_dest_arr, double* imp_response_arr, int sig_src_length, int imp_response_length)
{
	int i, j;
	//Initialization of the dest/result array
	for (int i = 0; i < (sig_src_length+imp_response_length); i++)
	{
		sig_dest_arr[i] = 0.0;
	}

	for (int i = 0; i < sig_src_length; i++)
	{
		for (int j = 0; j < imp_response_length; j++)
		{
			sig_dest_arr[i + j] += sig_src_arr[i]*imp_response_arr[j];
		}
	}

}
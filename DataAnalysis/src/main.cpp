#include "main.h"


int main()
{
	//signal_statistics_test();
	//convolution_test();
	//running_sum_test();
	//dft_test();
	//idft_test();
	//ecg_dft_idft_test();
	//r2p_test();
	//complex_dft_test();
	//low_pass_windowed_sinc_filter_test();
	//band_pass_test();

	//Cueing_offline_test();

	//Translational
	Cueing_online_test();
	return 0;
}

#pragma region Cueing Online 
void Cueing_online_test()
{

	bool log_data = true;
	double scale_factor = 1;// 0.0005;
	int data_idx = 0;
	double hp_cutoff = 0.089; 
	calc_kernel_high_pass(RT_KER_LEN, hp_cutoff, &Output_signal_rt_kernel[0], log_data);

	//logging to file
	std::string log_fldr = "log/rt_ay/";
	std::string delimiter = "\t"; //tab limited text file with 8 point precision
	std::string log_filename = "log_data_";
	log_filename += std::to_string(std::time(nullptr));
	log_filename += ".dat";
	std::fstream log_fptr;
	log_fptr.open(log_fldr + log_filename, std::fstream::in | std::fstream::out | std::fstream::app);


	for (int m = 0; m < RT_TOTAL_SIG_LEN; m++)
	{
		SP7Pose* pos = new SP7Pose{0,0,0,0,0,0};// initialization to avoid sending garbage value
		SP7Vel* vel = new SP7Vel{0,0,0,0,0,0};
		cueing_acceleration_online(xplane_ay_test3[m], xplane_t_test3[m], &Output_signal_rt_kernel[0], scale_factor, 2, &(pos->y), &(vel->vy));

		
		if (!log_data) return;
		//Preparing the data stream
		std::stringstream ss;
		ss.precision(8);// max to micro meter
		ss << std::fixed << t << delimiter << pos->x << delimiter <<
			pos->y << delimiter <<
			pos->z << delimiter <<
			pos->roll << delimiter <<
			pos->pitch << delimiter <<
			pos->yaw << delimiter <<
			vel->vx << delimiter << 
			vel->vy << delimiter << 
			vel->vz << delimiter << 
			vel->vroll << delimiter << 
			vel->vpitch << delimiter << 
			vel->vyaw;		
		
		if (log_fptr.is_open())
		{
			log_fptr << "\n";
			log_fptr << ss.str();
			
		}

		delete pos, vel;
		
	}
	log_fptr.close();
	std::cout << "Log file location: " << log_fldr << log_filename ;
}

void cueing_acceleration_online(double sig_acc_input, double sig_time, double* kernel, double scale_factor, int data_index, double* out_pos_, double* out_vel_)
{
	/* filter => scale => integrate => integrate => update prev*/
	//Convolving with designed kernel will filter the signal
	//Offestting should not affect the integration value as the value will be modified  based on that

	t = sig_time;
	acc_fltrd = Convolve_rt(&kernel[0], RT_KER_LEN, sig_acc_input, &Input_Buff[0]);
	acc_fltrd_scaled = scale_factor * acc_fltrd; //scale 
	velocity = Intergration_Trapezoidal(acc_fltrd_scaled, acc_fltrd_scaled_prev, velocity_prev, t_prev, t); //Integration
	position = Intergration_Trapezoidal(velocity, velocity_prev, position_prev, t_prev, t); //Integration
	position += SP7_ZERO_POSE[data_index];//Offseting with respect to zeropose of SP7

	//Updating the output values
	*(out_pos_) = position;
	*(out_vel_) = velocity;

	//Updating the previous to current
	t_prev = t;
	acc_fltrd_scaled_prev = acc_fltrd_scaled;
	velocity_prev = velocity;
	position_prev = position;

	
}
#pragma endregion

#pragma region Helper functions

double Intergration_Trapezoidal(double input_curr, double input_prev, double output_prev, double t_prev, double t_curr)
{
	return output_prev + ((t_curr - t_prev) * ((input_curr + input_prev) / 2)); // Trapedzoidal intergral
}
//double Intergration_Trapezoidal(double acc_curr, double acc_prev, double vel_prev, double t_prev, double t_curr)

void convolve(double* sig1, int sig1_len, double* sig2, int sig2_len, double* convolved_sig_)
{
	//Convolution is commutative
	bool isSig2Large = (sig1_len < sig2_len);

	if (isSig2Large)
	{
		for (int j = sig1_len; j < sig2_len; j++)
		{
			convolved_sig_[j] = 0;
			for (int i = 0; i < sig1_len; i++)
			{
				convolved_sig_[j] = convolved_sig_[j] + sig2[j - i] * sig1[i];

			}
		}
	}
	else
	{
		for (int j = sig2_len; j < sig1_len; j++)
		{
			convolved_sig_[j] = 0;
			for (int i = 0; i < sig2_len; i++)
			{
				convolved_sig_[j] = convolved_sig_[j] + sig1[j - i] * sig2[i];

			}
		}
	}

}//O((sig1_len+sig2_len)^2) //could be improved

// Real time convolution for one input signal and filter kernel
// Covolution with MAC(multiply and accumulate) and circultion buffering (shifting)
double Convolve_rt(double* h, int h_size, double x_in , double* x)
{
	// executes only when run first and then treats as a global variable
	//static int idx = h_size - 1;
	static int i = h_size - 1;

	x[i] = x_in;

	//Covolution
	double y_out = 0.0;
	for (int k = 0; k < h_size; k++)
		y_out += h[k] * x[(k + i) % h_size]; //MAC

	//idx++;
	//Signal shift circularly through  array x in time-reversed order
	i += h_size - 1; 
	i %= h_size;
	return y_out;
}

void negate_scale_rmpadding(double* filtered_sig, int kernel_len, int sig_len, int scale_factor, double* sig_out_)
{
	//Scale value calculation, removing padding data and negating the effect of spectral inversion
	int j = 0;
	int total_length = kernel_len + sig_len;
	int mid_idx = (total_length / 2) - 1;
	if (total_length % 2 != 0)
		mid_idx++;

	for (int i = mid_idx; i < mid_idx + sig_len; i++)
	{
		sig_out_[j] = scale_factor * (-filtered_sig[i]);
		j++;
	}
}

#pragma endregion


#pragma region Cueing Offline 

void Cueing_offline_test()
{
	//Translational
	double hp_cutoff = 0.089; //  
	double scale_factor = 1; //0.0005
	std::string log_fldr = "log/ax/";
	int data_idx = 0;
	cueing_acceleration(&xplane_ax_test3[0], hp_cutoff, scale_factor, log_fldr, data_idx);

	hp_cutoff = 0.089; //  
	scale_factor = 1; //0.0005
	log_fldr = "log/ay/";
	data_idx = 1;
	cueing_acceleration(&xplane_ay_test3[0], hp_cutoff, scale_factor, log_fldr, data_idx);
	
	hp_cutoff = 0.089; //  
	scale_factor = 1; //0.0005
	log_fldr = "log/az/";
	data_idx = 2;
	cueing_acceleration(&xplane_az_test3[0], hp_cutoff, scale_factor, log_fldr, data_idx);


	//Rotational
	double hpv_cutoff = 0.089; //  
	double scale_factor_v = 1; //0.025;
	log_fldr = "log/vroll/";
	data_idx = 3;
	cueing_velocity(&xplane_vroll_test3[0], hpv_cutoff, scale_factor_v, log_fldr, data_idx);

	hpv_cutoff = 0.089; //  
	scale_factor_v = 1; //0.01;
	log_fldr = "log/vpitch/";
	data_idx = 4;
	cueing_velocity(&xplane_vpitch_test3[0], hpv_cutoff, scale_factor_v, log_fldr, data_idx);

	hpv_cutoff = 0.089; //  
	scale_factor_v = 1; //0.05;
	log_fldr = "log/vyaw/";
	data_idx = 5;
	cueing_velocity(&xplane_vyaw_test3[0], hpv_cutoff, scale_factor_v, log_fldr, data_idx);
}

void cueing_velocity(double* sig_acc_input, double hp_cutoff, double scale_factor, std::string log_folder, int data_index)
{
	// frequency analysis of the input signal 
	calc_sig_dft(sig_acc_input, &Output_signal_xp_rex[0], &Output_signal_xp_imx[0], XP_SIG_LENGTH);
	get_dft_output_mag(&Output_signal_xp_mag[0], &Output_signal_xp_rex[0], &Output_signal_xp_imx[0], XP_SIG_LENGTH);
	rect_to_polar_conversion(&Output_signal_xp_rex[0], &Output_signal_xp_imx[0], &Output_signal_xp_pol_mag[0], &Output_signal_xp_pol_phase[0], XP_SIG_LENGTH / 2);

	// filtering 
	hp_windowed_sinc_fltr(&Output_signal_xp_kernel[0], hp_cutoff, KERNEL_LENGTH, sig_acc_input, &Output_signal_xp_dest[0], XP_SIG_LENGTH);

	//freq analysis of the output signal
	calc_sig_dft(&Output_signal_xp_dest[0], &Output_signal2_xp_rex[0], &Output_signal2_xp_imx[0], XP_SIG_LENGTH);
	rect_to_polar_conversion(&Output_signal2_xp_rex[0], &Output_signal2_xp_imx[0], &Output_signal2_xp_pol_mag[0], &Output_signal2_xp_pol_phase[0], XP_SIG_LENGTH / 2);

	//Scale value calculation
	double scale_value = 0.1;
	std::pair<double*, double*> minmax = std::minmax_element(std::begin(Output_signal_xp_dest), std::end(Output_signal_xp_dest));

	if (std::abs(*(minmax.first)) > std::abs(*(minmax.second)))
	{
		scale_value = (1 * scale_factor) / std::abs(*(minmax.first));
	}
	else
	{
		scale_value = (1 * scale_factor) / std::abs(*(minmax.second));
	}

	//Scaling and negating the effect of spectral inversion
	int j = 0;
	for (int i = KERNEL_LENGTH / 2; i < XP_SIG_LENGTH - (KERNEL_LENGTH / 2); i++)
	{

		Output_signal_xp_dest2[j] = scale_value * (-Output_signal_xp_dest[i]);//negate the filtered to get proper signed signal. (effect of spectral inversion)
		j++;
	}


	//Integration
	Output_signal_xp_integral[0] = 0;
	

	for (int i = 1; i < XP_SIG_LENGTH; i++)
	{
		Output_signal_xp_integral[i] = Intergration_Trapezoidal(Output_signal_xp_dest2[i], Output_signal_xp_dest2[i - 1], Output_signal_xp_integral[i - 1], xplane_t_test3[i - 1], xplane_t_test3[i]);
	}

	//Offseting with respect to zeropose of SP7
	//THIS OFFSET SHOULD NOT AFFECT THE INTEGRATION VALUE AS THE VALUE WILL BE MODIFIED BASED ON THAT
	for (int i = 0; i < XP_SIG_LENGTH; i++)
	{
		Output_signal_xp_integral[i] += SP7_ZERO_POSE[data_index];
	}


	//Data log to file
	std::fstream input_sig_fptr, output_sig_rex_fptr, output_sig_imx_fptr, output_sig_mag_fptr, output_sig_pol_mag_fptr, output_sig_pol_phase_fptr, output_sig_hp_dest_fptr, output_sig_hp_kernel_fptr,
		output_sig2_rex_fptr, output_sig2_imx_fptr, output_sig2_pol_mag_fptr, output_sig2_pol_phase_fptr,
		output_sig_xp_integral_ftr, time_fptr;



	input_sig_fptr.open(log_folder + "input_signal.dat", std::ios::out);
	output_sig_rex_fptr.open(log_folder + "output_signal_dft_rex.dat", std::ios::out);
	output_sig_imx_fptr.open(log_folder + "output_signal_dft_imx.dat", std::ios::out);
	output_sig_mag_fptr.open(log_folder + "output_signal_dft_mag.dat", std::ios::out);
	output_sig_pol_mag_fptr.open(log_folder + "output_signal_pol_mag.dat", std::ios::out);
	output_sig_pol_phase_fptr.open(log_folder + "output_signal_pol_phase.dat", std::ios::out);
	output_sig_hp_dest_fptr.open(log_folder + "output_sig_hp_dest.dat", std::ios::out);
	output_sig_hp_kernel_fptr.open(log_folder + "output_sig_hp_kernel.dat", std::ios::out);
	//Output Freq analysis
	output_sig2_rex_fptr.open(log_folder + "output_signal2_dft_rex.dat", std::ios::out);
	output_sig2_imx_fptr.open(log_folder + "output_signal2_dft_imx.dat", std::ios::out);
	output_sig2_pol_mag_fptr.open(log_folder + "output_signal2_pol_mag.dat", std::ios::out);
	output_sig2_pol_phase_fptr.open(log_folder + "output_signal2_pol_phase.dat", std::ios::out);
	//Intregration
	output_sig_xp_integral_ftr.open(log_folder + "output_sig_xp_integral.dat", std::ios::out);
	

	time_fptr.open(log_folder + "time.dat", std::ios::out);

	if (input_sig_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH; i++)
		{
			input_sig_fptr << "\n";
			input_sig_fptr << *(sig_acc_input + i);//column vector
		}
		input_sig_fptr.close();
	}

	if (output_sig_rex_fptr.is_open())
	{
		for (int j = 0; j < XP_SIG_LENGTH / 2; j++)
		{
			output_sig_rex_fptr << "\n";
			output_sig_rex_fptr << Output_signal_xp_rex[j];//column vector
		}
		output_sig_rex_fptr.close();
	}

	if (output_sig_imx_fptr.is_open())
	{
		for (int k = 0; k < XP_SIG_LENGTH / 2; k++)
		{
			output_sig_imx_fptr << "\n";
			output_sig_imx_fptr << Output_signal_xp_imx[k];//column vector
		}
		output_sig_imx_fptr.close();
	}

	if (output_sig_mag_fptr.is_open())
	{
		for (int l = 0; l < XP_SIG_LENGTH / 2; l++)
		{
			output_sig_mag_fptr << "\n";
			output_sig_mag_fptr << Output_signal_xp_mag[l];//column vector
		}
		output_sig_mag_fptr.close();
	}

	if (output_sig_pol_mag_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH / 2; i++)
		{
			output_sig_pol_mag_fptr << "\n";
			output_sig_pol_mag_fptr << Output_signal_xp_pol_mag[i];//column vector
		}
		output_sig_pol_mag_fptr.close();
	}

	if (output_sig_pol_phase_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH / 2; i++)
		{
			output_sig_pol_phase_fptr << "\n";
			output_sig_pol_phase_fptr << Output_signal_xp_pol_phase[i];//column vector
		}
		output_sig_pol_phase_fptr.close();
	}


	if (output_sig_hp_dest_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH; i++)
		{
			output_sig_hp_dest_fptr << "\n";
			output_sig_hp_dest_fptr << Output_signal_xp_dest2[i];//column vector
		}
		output_sig_hp_dest_fptr.close();
	}


	if (output_sig_hp_kernel_fptr.is_open())
	{
		for (int i = 0; i < KERNEL_LENGTH; i++)
		{
			output_sig_hp_kernel_fptr << "\n";
			output_sig_hp_kernel_fptr << Output_signal_xp_kernel[i];//column vector
		}
		output_sig_hp_kernel_fptr.close();
	}


	//
	if (output_sig2_rex_fptr.is_open())
	{
		for (int j = 0; j < XP_SIG_LENGTH / 2; j++)
		{
			output_sig2_rex_fptr << "\n";
			output_sig2_rex_fptr << Output_signal2_xp_rex[j];//column vector
		}
		output_sig2_rex_fptr.close();
	}

	if (output_sig2_imx_fptr.is_open())
	{
		for (int k = 0; k < XP_SIG_LENGTH / 2; k++)
		{
			output_sig2_imx_fptr << "\n";
			output_sig2_imx_fptr << Output_signal2_xp_imx[k];//column vector
		}
		output_sig2_imx_fptr.close();
	}


	if (output_sig2_pol_mag_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH / 2; i++)
		{
			output_sig2_pol_mag_fptr << "\n";
			output_sig2_pol_mag_fptr << Output_signal2_xp_pol_mag[i];//column vector
		}
		output_sig2_pol_mag_fptr.close();
	}

	if (output_sig2_pol_phase_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH / 2; i++)
		{
			output_sig2_pol_phase_fptr << "\n";
			output_sig2_pol_phase_fptr << Output_signal2_xp_pol_phase[i];//column vector
		}
		output_sig2_pol_phase_fptr.close();
	}

	//
	if (output_sig_xp_integral_ftr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH; i++)
		{
			output_sig_xp_integral_ftr << "\n";
			output_sig_xp_integral_ftr << Output_signal_xp_integral[i];//column vector
		}
		output_sig_xp_integral_ftr.close();
	}

	if (time_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH; i++)
		{
			time_fptr << "\n";
			time_fptr << xplane_t_test3[i];//column vector
		}
		time_fptr.close();
	}

}
void cueing_acceleration(double * sig_acc_input, double hp_cutoff, double scale_factor, std::string log_folder, int data_index)
{
	// frequency analysis of the input signal 
	calc_sig_dft(sig_acc_input, &Output_signal_xp_rex[0], &Output_signal_xp_imx[0], XP_SIG_LENGTH);
	get_dft_output_mag(&Output_signal_xp_mag[0], &Output_signal_xp_rex[0], &Output_signal_xp_imx[0], XP_SIG_LENGTH);
	rect_to_polar_conversion(&Output_signal_xp_rex[0], &Output_signal_xp_imx[0], &Output_signal_xp_pol_mag[0], &Output_signal_xp_pol_phase[0], XP_SIG_LENGTH / 2);

	// filtering 
	hp_windowed_sinc_fltr(&Output_signal_xp_kernel[0], hp_cutoff, KERNEL_LENGTH, sig_acc_input, &Output_signal_xp_dest[0], XP_SIG_LENGTH);

	//freq analysis of the output signal
	calc_sig_dft(&Output_signal_xp_dest[0], &Output_signal2_xp_rex[0], &Output_signal2_xp_imx[0], XP_SIG_LENGTH);
	rect_to_polar_conversion(&Output_signal2_xp_rex[0], &Output_signal2_xp_imx[0], &Output_signal2_xp_pol_mag[0], &Output_signal2_xp_pol_phase[0], XP_SIG_LENGTH / 2);

	//Scale value calculation
	double scale_value = 0.1;
	std::pair<double*, double*> minmax = std::minmax_element(std::begin(Output_signal_xp_dest), std::end(Output_signal_xp_dest));

	if (std::abs(*(minmax.first)) > std::abs(*(minmax.second)))
	{
		scale_value = (1 * scale_factor) / std::abs(*(minmax.first));
	}
	else
	{
		scale_value = (1 * scale_factor) / std::abs(*(minmax.second));
	}
	
	//Scaling and negating the effect of spectral inversion
	int j = 0;
	for (int i = KERNEL_LENGTH / 2; i < XP_SIG_LENGTH - (KERNEL_LENGTH / 2); i++)
	{

		Output_signal_xp_dest2[j] = scale_value * (-Output_signal_xp_dest[i]);//negate the filtered to get proper signed signal. (effect of spectral inversion)
		j++;
	}


	//Integration
	Output_signal_xp_integral[0] = 0;
	Output_signal_xp_double_integral[0] = 0;

	for (int i = 1; i < XP_SIG_LENGTH; i++)
	{
		Output_signal_xp_integral[i] = Intergration_Trapezoidal(Output_signal_xp_dest2[i], Output_signal_xp_dest2[i - 1], Output_signal_xp_integral[i - 1], xplane_t_test3[i - 1], xplane_t_test3[i]);
	}


	for (int i = 1; i < XP_SIG_LENGTH; i++)
	{
		Output_signal_xp_double_integral[i] =  Intergration_Trapezoidal(Output_signal_xp_integral[i], Output_signal_xp_integral[i - 1], Output_signal_xp_double_integral[i - 1], xplane_t_test3[i - 1], xplane_t_test3[i]);
	}

	//Offseting with respect to zeropose of SP7
	//THIS OFFSET SHOULD NOT AFFECT THE INTEGRATION VALUE AS THE VALUE WILL BE MODIFIED BASED ON THAT
	for (int i = 0; i < XP_SIG_LENGTH; i++)
	{
		Output_signal_xp_double_integral[i] += SP7_ZERO_POSE[data_index];
	}


	//Data log to file
	std::fstream input_sig_fptr, output_sig_rex_fptr, output_sig_imx_fptr, output_sig_mag_fptr, output_sig_pol_mag_fptr, output_sig_pol_phase_fptr, output_sig_hp_dest_fptr, output_sig_hp_kernel_fptr,
		output_sig2_rex_fptr, output_sig2_imx_fptr, output_sig2_pol_mag_fptr, output_sig2_pol_phase_fptr,
		output_sig_xp_integral_ftr, output_sig_xp_integral2_ftr, time_fptr;

	

	input_sig_fptr.open(log_folder+"input_signal.dat", std::ios::out);
	output_sig_rex_fptr.open(log_folder + "output_signal_dft_rex.dat", std::ios::out);
	output_sig_imx_fptr.open(log_folder + "output_signal_dft_imx.dat", std::ios::out);
	output_sig_mag_fptr.open(log_folder + "output_signal_dft_mag.dat", std::ios::out);
	output_sig_pol_mag_fptr.open(log_folder + "output_signal_pol_mag.dat", std::ios::out);
	output_sig_pol_phase_fptr.open(log_folder + "output_signal_pol_phase.dat", std::ios::out);
	output_sig_hp_dest_fptr.open(log_folder + "output_sig_hp_dest.dat", std::ios::out);
	output_sig_hp_kernel_fptr.open(log_folder + "output_sig_hp_kernel.dat", std::ios::out);
	//Output Freq analysis
	output_sig2_rex_fptr.open(log_folder + "output_signal2_dft_rex.dat", std::ios::out);
	output_sig2_imx_fptr.open(log_folder + "output_signal2_dft_imx.dat", std::ios::out);
	output_sig2_pol_mag_fptr.open(log_folder + "output_signal2_pol_mag.dat", std::ios::out);
	output_sig2_pol_phase_fptr.open(log_folder + "output_signal2_pol_phase.dat", std::ios::out);
	//Intregration
	output_sig_xp_integral_ftr.open(log_folder + "output_sig_xp_integral.dat", std::ios::out);
	output_sig_xp_integral2_ftr.open(log_folder + "output_sig_xp_integral2.dat", std::ios::out);

	time_fptr.open(log_folder + "time.dat", std::ios::out);

	if (input_sig_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH; i++)
		{
			input_sig_fptr << "\n";
			input_sig_fptr << *(sig_acc_input+i);//column vector
		}
		input_sig_fptr.close();
	}

	if (output_sig_rex_fptr.is_open())
	{
		for (int j = 0; j < XP_SIG_LENGTH / 2; j++)
		{
			output_sig_rex_fptr << "\n";
			output_sig_rex_fptr << Output_signal_xp_rex[j];//column vector
		}
		output_sig_rex_fptr.close();
	}

	if (output_sig_imx_fptr.is_open())
	{
		for (int k = 0; k < XP_SIG_LENGTH / 2; k++)
		{
			output_sig_imx_fptr << "\n";
			output_sig_imx_fptr << Output_signal_xp_imx[k];//column vector
		}
		output_sig_imx_fptr.close();
	}

	if (output_sig_mag_fptr.is_open())
	{
		for (int l = 0; l < XP_SIG_LENGTH / 2; l++)
		{
			output_sig_mag_fptr << "\n";
			output_sig_mag_fptr << Output_signal_xp_mag[l];//column vector
		}
		output_sig_mag_fptr.close();
	}

	if (output_sig_pol_mag_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH / 2; i++)
		{
			output_sig_pol_mag_fptr << "\n";
			output_sig_pol_mag_fptr << Output_signal_xp_pol_mag[i];//column vector
		}
		output_sig_pol_mag_fptr.close();
	}

	if (output_sig_pol_phase_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH / 2; i++)
		{
			output_sig_pol_phase_fptr << "\n";
			output_sig_pol_phase_fptr << Output_signal_xp_pol_phase[i];//column vector
		}
		output_sig_pol_phase_fptr.close();
	}


	if (output_sig_hp_dest_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH; i++)
		{
			output_sig_hp_dest_fptr << "\n";
			output_sig_hp_dest_fptr << Output_signal_xp_dest2[i];//column vector
		}
		output_sig_hp_dest_fptr.close();
	}


	if (output_sig_hp_kernel_fptr.is_open())
	{
		for (int i = 0; i < KERNEL_LENGTH; i++)
		{
			output_sig_hp_kernel_fptr << "\n";
			output_sig_hp_kernel_fptr << Output_signal_xp_kernel[i];//column vector
		}
		output_sig_hp_kernel_fptr.close();
	}


	//
	if (output_sig2_rex_fptr.is_open())
	{
		for (int j = 0; j < XP_SIG_LENGTH / 2; j++)
		{
			output_sig2_rex_fptr << "\n";
			output_sig2_rex_fptr << Output_signal2_xp_rex[j];//column vector
		}
		output_sig2_rex_fptr.close();
	}

	if (output_sig2_imx_fptr.is_open())
	{
		for (int k = 0; k < XP_SIG_LENGTH / 2; k++)
		{
			output_sig2_imx_fptr << "\n";
			output_sig2_imx_fptr << Output_signal2_xp_imx[k];//column vector
		}
		output_sig2_imx_fptr.close();
	}


	if (output_sig2_pol_mag_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH / 2; i++)
		{
			output_sig2_pol_mag_fptr << "\n";
			output_sig2_pol_mag_fptr << Output_signal2_xp_pol_mag[i];//column vector
		}
		output_sig2_pol_mag_fptr.close();
	}

	if (output_sig2_pol_phase_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH / 2; i++)
		{
			output_sig2_pol_phase_fptr << "\n";
			output_sig2_pol_phase_fptr << Output_signal2_xp_pol_phase[i];//column vector
		}
		output_sig2_pol_phase_fptr.close();
	}

	//
	if (output_sig_xp_integral_ftr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH; i++)
		{
			output_sig_xp_integral_ftr << "\n";
			output_sig_xp_integral_ftr << Output_signal_xp_integral[i];//column vector
		}
		output_sig_xp_integral_ftr.close();
	}

	if (output_sig_xp_integral2_ftr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH; i++)
		{
			output_sig_xp_integral2_ftr << "\n";
			output_sig_xp_integral2_ftr << Output_signal_xp_double_integral[i];//column vector
		}
		output_sig_xp_integral2_ftr.close();
	}

	if (time_fptr.is_open())
	{
		for (int i = 0; i < XP_SIG_LENGTH; i++)
		{
			time_fptr << "\n";
			time_fptr << xplane_t_test3[i];//column vector
		}
		time_fptr.close();
	}

}
#pragma endregion 

#pragma region DSP Tutorial Methods
void band_pass_test()
{
	double lower_cutoff_freq = 0.21; // 10kHz  ////10kHz.002; // 0.1kHz
	double upper_cutoff_freq = 0.33; // 16kHz  ////0.11; //5.28kHz

	bp_windowed_sinc_fltr(
		&state_lw_cutoff_buff[0],
		&state_up_cutoff_buff[0],
		&Output_signal_bp_kernel[0],
		lower_cutoff_freq,
		upper_cutoff_freq,
		KERNEL_LENGTH,
		&InputSignal_f32_1kHz_15kHz[0],
		&Output_signal_bp_dest[0],
		SIG_LENGTH);


	//Prepare the data : Convert the double array to std::vector
	std::vector<double> d1(std::begin(InputSignal_f32_1kHz_15kHz), std::end(InputSignal_f32_1kHz_15kHz));
	std::vector<double> d2(std::begin(Output_signal_bp_dest), std::end(Output_signal_bp_dest));
	std::vector<double> d3(std::begin(Output_signal_bp_kernel), std::end(Output_signal_bp_kernel));
	std::vector<double> d4(std::begin(InputSignal_f32_1kHz_15kHz), std::end(InputSignal_f32_1kHz_15kHz));
	plot(d1, d2, d3, d4);
}


void low_pass_windowed_sinc_filter_test()
{
	double cutoff_freq = 0.2; //10kHz:refer the low_pass_windowed_sinc_filter.cpp comment about normalization

	lp_windowed_sinc_fltr(&InputSignal_f32_1kHz_15kHz[0], &Output_signal_lp_dest[0], &Output_signal_lp_kernel[0], cutoff_freq, KERNEL_LENGTH, SIG_LENGTH
	);

	//Prepare the data : Convert the double array to std::vector
	std::vector<double> d1(std::begin(InputSignal_f32_1kHz_15kHz), std::end(InputSignal_f32_1kHz_15kHz));
	std::vector<double> d2(std::begin(Output_signal_lp_dest), std::end(Output_signal_lp_dest));
	std::vector<double> d3(std::begin(Output_signal_lp_kernel), std::end(Output_signal_lp_kernel));
	std::vector<double> d4(std::begin(InputSignal_f32_1kHz_15kHz), std::end(InputSignal_f32_1kHz_15kHz));
	plot(d1, d2, d3, d4);
}

void complex_dft_test() 
{
	complex_dft(&_501pts_20Hz_sig_rex[0], &_501pts_20Hz_sig_imx[0], &Output_signal_cdft_freq_rex[0], &Output_signal_cdft_freq_imx[0], COMPLEX_SIG_LENGTH);

	//Prepare the data : Convert the double array to std::vector
	std::vector<double> d1(std::begin(_501pts_20Hz_sig_rex), std::end(_501pts_20Hz_sig_rex));
	std::vector<double> d2(std::begin(_501pts_20Hz_sig_imx), std::end(_501pts_20Hz_sig_imx));
	std::vector<double> d3(std::begin(Output_signal_cdft_freq_rex), std::end(Output_signal_cdft_freq_rex));
	std::vector<double> d4(std::begin(Output_signal_cdft_freq_imx), std::end(Output_signal_cdft_freq_imx));
	plot(d1, d2, d3, d4);
}

void r2p_test()
{
	rect_to_polar_conversion(&_320_pts_ecg_REX[0], &_320_pts_ecg_IMX[0], &Output_signal_pol_mag[0], &Output_signal_pol_phase[0], POL_SIG_LENGTH);

	//Prepare the data : Convert the double array to std::vector
	std::vector<double> d1(std::begin(_320_pts_ecg_REX), std::end(_320_pts_ecg_REX));
	std::vector<double> d2(std::begin(_320_pts_ecg_IMX), std::end(_320_pts_ecg_IMX));
	std::vector<double> d3(std::begin(Output_signal_pol_mag), std::end(Output_signal_pol_mag));
	std::vector<double> d4(std::begin(Output_signal_pol_phase), std::end(Output_signal_pol_phase));
	plot(d1, d2, d3, d4);
}

void ecg_dft_idft_test()
{
	calc_sig_dft(&_640_points_ecg_[0], &Output_signal_ecg_dft_rex[0], &Output_signal_ecg_dft_imx[0], ECG_SIG_LENGTH);
	calc_sig_idft((double*)&Output_signal_ecg_idft[0], (double*)&Output_signal_ecg_dft_rex[0], (double*)&Output_signal_ecg_dft_imx[0], ECG_SIG_LENGTH);

	//Prepare the data : Convert the double array to std::vector
	std::vector<double> d1(std::begin(_640_points_ecg_), std::end(_640_points_ecg_));
	std::vector<double> d2(std::begin(Output_signal_ecg_dft_rex), std::end(Output_signal_ecg_dft_rex));
	std::vector<double> d3(std::begin(Output_signal_ecg_dft_imx), std::end(Output_signal_ecg_dft_imx));
	std::vector<double> d4(std::begin(Output_signal_ecg_idft), std::end(Output_signal_ecg_idft));
	plot(d1, d2, d3, d4);
}

void idft_test()
{
	calc_sig_dft(&InputSignal_f32_1kHz_15kHz[0], &Output_signal_dft_rex[0], &Output_signal_dft_imx[0], SIG_LENGTH);
	calc_sig_idft((double*)&Output_signal_idft[0], (double*)&Output_signal_dft_rex[0], (double*)&Output_signal_dft_imx[0], SIG_LENGTH);

	//Prepare the data : Convert the double array to std::vector
	std::vector<double> d1(std::begin(InputSignal_f32_1kHz_15kHz), std::end(InputSignal_f32_1kHz_15kHz));
	std::vector<double> d2(std::begin(Output_signal_dft_rex), std::end(Output_signal_dft_rex));
	std::vector<double> d3(std::begin(Output_signal_dft_imx), std::end(Output_signal_dft_imx));
	std::vector<double> d4(std::begin(Output_signal_idft), std::end(Output_signal_idft));
	plot(d1, d2, d3, d4);
}

void plot(std::vector<double> _d1, std::vector<double> _d2, std::vector<double> _d3, std::vector<double> _d4)
{
	//cw d1, d2, d3, d4 : This function makes a multiplot with 4 graphs
	//ploting of gnuplot => sending data via temporary files
	Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");
	gp << "set size 1,1\n";
	gp << "set multiplot\n";
	gp << "set size 0.5,0.5\n";
	gp << "set origin 0,0\n";
	gp << "plot" << gp.file1d(_d1) << "with lines lc rgb 'black' title 'd1'\n";
	gp << "set origin 0,0.5\n";
	gp << "plot" << gp.file1d(_d2) << "with lines lc rgb 'red' title 'd2'" << std::endl;
	gp << "set origin 0.5,0.5\n";
	gp << "plot" << gp.file1d(_d3) << "with lines lc rgb 'blue' title 'd3'\n";
	gp << "set origin 0.5,0\n";
	gp << "plot" << gp.file1d(_d4) << "with lines lc rgb 'green' title 'd4'" << std::endl;
	//gp << "pause -1 " << std::endl;
		
#ifdef _WIN32
	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Press enter to exit." << std::endl;
	std::cin.get();
#endif
	gp.clearTmpfiles();
}

void dft_test()
{
	std::fstream input_sig_fptr, output_sig_rex_fptr, output_sig_imx_fptr, output_sig_mag_fptr;

	calc_sig_dft(&InputSignal_f32_1kHz_15kHz[0], &Output_signal_dft_rex[0], &Output_signal_dft_imx[0], SIG_LENGTH);

	get_dft_output_mag(&Output_signal_dft_mag[0], &Output_signal_dft_rex[0], &Output_signal_dft_imx[0], SIG_LENGTH);

	input_sig_fptr.open("input_signal.dat", std::ios::out);
	output_sig_rex_fptr.open("output_signal_dft_rex.dat", std::ios::out);
	output_sig_imx_fptr.open("output_signal_dft_imx.dat", std::ios::out);
	output_sig_mag_fptr.open("output_signal_dft_mag.dat", std::ios::out);

	if (input_sig_fptr.is_open())
	{
		for (int i = 0; i < SIG_LENGTH; i++)
		{
			input_sig_fptr << "\n";
			input_sig_fptr << InputSignal_f32_1kHz_15kHz[i];//column vector
		}
		input_sig_fptr.close();
	}

	if (output_sig_rex_fptr.is_open())
	{
		for (int j = 0; j < SIG_LENGTH / 2; j++)
		{
			output_sig_rex_fptr << "\n";
			output_sig_rex_fptr << Output_signal_dft_rex[j];//column vector
		}
		output_sig_rex_fptr.close();
	}

	if (output_sig_imx_fptr.is_open())
	{
		for (int k = 0; k < SIG_LENGTH / 2; k++)
		{
			output_sig_imx_fptr << "\n";
			output_sig_imx_fptr << Output_signal_dft_imx[k];//column vector
		}
		output_sig_imx_fptr.close();
	}

	if (output_sig_mag_fptr.is_open())
	{
		for (int l = 0; l < SIG_LENGTH / 2; l++)
		{
			output_sig_mag_fptr << "\n";
			output_sig_mag_fptr << Output_signal_dft_mag[l];//column vector
		}
		output_sig_mag_fptr.close();
	}

}

void running_sum_test()
{
	std::fstream input_sig_fptr, output_sig_fptr;
	calc_running_sum(&InputSignal_f32_1kHz_15kHz[0], &Output_signal_running_sum[0], SIG_LENGTH);
	input_sig_fptr.open("input_signal.dat", std::ios::out);
	output_sig_fptr.open("output_signal.dat", std::ios::out);

	if (input_sig_fptr.is_open())
	{
		for (int i = 0; i < SIG_LENGTH; i++)
		{
			input_sig_fptr << "\n";
			input_sig_fptr << InputSignal_f32_1kHz_15kHz[i];//column vector
		}
		input_sig_fptr.close();
	}


	if (output_sig_fptr.is_open())
	{
		for (int k = 0; k < SIG_LENGTH; k++)
		{
			output_sig_fptr << "\n";
			output_sig_fptr << Output_signal_running_sum[k];//column vector
		}
		output_sig_fptr.close();
	}
}

void calc_running_sum(double* sig_src_arr, double* sig_dest_arr, int sig_length)
{
	for (int i = 0; i < sig_length; i++)
	{
		sig_dest_arr[i] = sig_dest_arr[i - 1] + sig_src_arr[i];
	}
}

void signal_statistics_test()
{
	MEAN = calc_signal_mean(&InputSignal_f32_1kHz_15kHz[0], SIG_LENGTH);
	printf("\n\nMean = %f\n\n\n", MEAN);

	VARIANCE = calc_signal_variance(&InputSignal_f32_1kHz_15kHz[0], MEAN, SIG_LENGTH);
	printf("\n\nVariance = %f\n\n\n", VARIANCE);

	STD = calc_signal_std(VARIANCE);
	printf("\n\nStandard Deviation = %f\n\n\n", STD);
}

void convolution_test()
{
	std::fstream input_sig_fptr, imp_rsp_fptr, output_sig_fptr;

	convolution(&InputSignal_f32_1kHz_15kHz[0], &Output_signal[0], &Impulse_response[0], SIG_LENGTH, IMP_RSP_LENGTH);


	input_sig_fptr.open("input_signal.dat", std::ios::out);
	imp_rsp_fptr.open("impulse_response.dat", std::ios::out);
	output_sig_fptr.open("output_signal.dat", std::ios::out);

	if (input_sig_fptr.is_open())
	{
		for (int i = 0; i < SIG_LENGTH; i++)
		{
			input_sig_fptr << "\n";
			input_sig_fptr << InputSignal_f32_1kHz_15kHz[i];//column vector
		}
		input_sig_fptr.close();
	}

	if (imp_rsp_fptr.is_open())
	{
		for (int j = 0; j < IMP_RSP_LENGTH; j++)
		{
			imp_rsp_fptr << "\n";
			imp_rsp_fptr << Impulse_response[j];//column vector
		}
		imp_rsp_fptr.close();
	}

	if (output_sig_fptr.is_open())
	{
		for (int k = 0; k < SIG_LENGTH + IMP_RSP_LENGTH; k++)
		{
			output_sig_fptr << "\n";
			output_sig_fptr << Output_signal[k];//column vector
		}
		output_sig_fptr.close();
	}
}
#pragma endregion 



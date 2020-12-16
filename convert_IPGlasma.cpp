#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

const string file_mode = "Tmunu";

//int ix, iy;
//double T00, Txx, Tyy, Tzz, T0x, T0y, T0z, Txy, Tyz, Txz;

void do_Landau_matching(
		const double T00_in, const double T0x_in, const double T0y_in, const double T0z_in,
		const double Txx_in, const double Txy_in, const double Txz_in,
		const double Tyy_in, const double Tyz_in, const double Tzz_in );

int main(int argc, char *argv[])
{
	if (argc < 2) exit(8);

	// read path to input file from command line
	string path_to_file = string(argv[1]);
	
	// then read in file itself
	ifstream infile(path_to_file.c_str());
	if (infile.is_open())
	{
		string line;
		int lineCount = 0;
		if ( file_mode == "e_umu_pimunu" )
		while ( getline (infile, line) )
		{
			// skip header line
			if (lineCount++ == 1) continue;

			/*istringstream iss(line);
			iss >> eta >> x >> y >> e
				>> utau >> ux >> uy >> ueta
				>> resultpi00 >> resultpi0x >> resultpi0y >> resultpi0eta
				>> resultpixx >> resultpixy >> resultpixeta
				>> resultpiyy >> resultpiyeta >> resultpietaeta;*/
	
		}
		else if ( file_mode == "Tmunu" )
		while ( getline (infile, line) )
		{
			// skip header line
			//if (lineCount++ == 1) continue;

int ix, iy;
double T00, Txx, Tyy, Tzz, T0x, T0y, T0z, Txy, Tyz, Txz;

			istringstream iss(line);
			iss >> ix >> iy >> T00 >> Txx >> Tyy >> Tzz
				>> T0x >> T0y >> T0z >> Txy >> Tyz >> Txz;
			/*T00 = 268.131;
			Txx = 154.278;
			Tyy = 175.396;
			Tzz = -61.5434;
			T0x = 45.0853;
			T0y = -18.6764;
			T0z = 10.5892;
			Txy = 28.6759;
			Tyz = 22.1628;
			Txz = 0.947038;*/
			const double tau0 = 0.6;	//fm/c
			
			cout << "********************************************************************************" << endl;
			do_Landau_matching( T00, -T0x, -T0y, -T0z, Txx, -Txy, -Txz, Tyy, -Tyz, Tzz );
			cout << "********************************************************************************" << endl;
			do_Landau_matching( T00, -T0x, -T0y, 0.0, Txx, -Txy, 0.0, Tyy, 0.0, Tzz );
			cout << "********************************************************************************" << endl;


if (true) exit(8);
		}
		else
			exit(8);
	}
	
	infile.close();

	return (0);
}





void do_Landau_matching(
		const double T00_in, const double T0x_in, const double T0y_in, const double T0z_in,
		const double Txx_in, const double Txy_in, const double Txz_in,
		const double Tyy_in, const double Tyz_in, const double Tzz_in )
{
	/*double m[16] = { -T00_in, -T0x_in, -T0y_in, -T0z_in,
                      T0x_in,  Txx_in,  Txy_in,  Txz_in,
                      T0y_in,  Txy_in,  Tyy_in,  Tyz_in,
                      T0z_in,  Txz_in,  Tyz_in,  Tzz_in };
	double copy[16] = { -T00_in, -T0x_in, -T0y_in, -T0z_in,
                      T0x_in,  Txx_in,  Txy_in,  Txz_in,
                      T0y_in,  Txy_in,  Tyy_in,  Tyz_in,
                      T0z_in,  Txz_in,  Tyz_in,  Tzz_in };*/
	double m[16] = {  T00_in, -T0x_in, -T0y_in, -T0z_in,
                      T0x_in, -Txx_in, -Txy_in, -Txz_in,
                      T0y_in, -Txy_in, -Tyy_in, -Tyz_in,
                      T0z_in, -Txz_in, -Tyz_in, -Tzz_in };
	double copy[16] = {  T00_in, -T0x_in, -T0y_in, -T0z_in,
                         T0x_in, -Txx_in, -Txy_in, -Txz_in,
                         T0y_in, -Txy_in, -Tyy_in, -Tyz_in,
                         T0z_in, -Txz_in, -Tyz_in, -Tzz_in };


	double T_in[16] = { T00_in, T0x_in, T0y_in, T0z_in,
                        T0x_in, Txx_in, Txy_in, Txz_in,
                        T0y_in, Txy_in, Tyy_in, Tyz_in,
                        T0z_in, Txz_in, Tyz_in, Tzz_in };

	gsl_vector_complex *eval = gsl_vector_complex_alloc(4);
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc(4, 4);

	gsl_matrix_view mat = gsl_matrix_view_array(m, 4, 4);
	gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(4);
	int success = gsl_eigen_nonsymmv (&mat.matrix, eval, evec, w);
	gsl_eigen_nonsymmv_free(w);

	// sort by magnitude first
	gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

{
    int i, j;

    for (i = 0; i < 4; i++)
      {
        gsl_complex eval_i
           = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i
           = gsl_matrix_complex_column (evec, i);

		double ev[4];
		for (j = 0; j < 4; ++j)
		{
			gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
			//printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
			ev[j] = GSL_REAL(z);
		}
		if ( ev[0]*ev[0] - ev[1]*ev[1] - ev[2]*ev[2] - ev[3]*ev[3] > 0.0 )
		{
			/*printf ("eigenvalue = %g + %gi\n",
			GSL_REAL(eval_i), GSL_IMAG(eval_i));
			printf ("eigenvector = \n");*/
			double ev[4];
			for (j = 0; j < 4; ++j)
			{
				gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
				//printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
				ev[j] = GSL_REAL(z);
			}
			/*
			cout << "Test: " << ev[0]*ev[0] - ev[1]*ev[1] - ev[2]*ev[2] - ev[3]*ev[3] << endl << endl << endl;
			cout << "Check:" << endl;
			cout << ev[0]*copy[0]+ev[1]*copy[1]+ev[2]*copy[2]+ev[3]*copy[3] << " == " << GSL_REAL(eval_i) * ev[0] << endl;
			cout << ev[0]*copy[4]+ev[1]*copy[5]+ev[2]*copy[6]+ev[3]*copy[7] << " == " << GSL_REAL(eval_i) * ev[1] << endl;
			cout << ev[0]*copy[8]+ev[1]*copy[9]+ev[2]*copy[10]+ev[3]*copy[11] << " == " << GSL_REAL(eval_i) * ev[2] << endl;
			cout << ev[0]*copy[12]+ev[1]*copy[13]+ev[2]*copy[14]+ev[3]*copy[15] << " == " << GSL_REAL(eval_i) * ev[3] << endl;
			*/

			const double norm = sqrt(abs(ev[0]*ev[0] - ev[1]*ev[1] - ev[2]*ev[2] - ev[3]*ev[3]));
			for (int iev = 0; iev < 4; iev++) ev[iev] /= norm;

			if (ev[0] < 0.0) for (int iev = 0; iev < 4; iev++) ev[iev] *= -1.0;

			const double trace = T00_in - Txx_in - Tyy_in - Tzz_in;
			const double e = abs(GSL_REAL(eval_i));
			const double P = (e - trace)/3.0;
			/*cout << "Trace: " << trace << endl;
			cout << "e: " << e << endl;
			cout << "P: " << P << endl;
			cout << "u: " << endl;
			for (int ii = 0; ii < 4; ii++)
				cout << ev[ii] << "   ";
			cout << endl << endl;
			cout << "T:" << endl;
			for (int ii = 0; ii < 4; ii++)
			{
				for (int jj = 0; jj < 4; jj++)
					cout << T_in[ii*4+jj] << "   ";
				cout << endl;
			}			cout << endl << endl;
			cout << "pi:" << endl;
			for (int ii = 0; ii < 4; ii++)
			{
				for (int jj = 0; jj < 4; jj++)
				{
					int g = -static_cast<int>( ii == jj );
					if (ii==0 && jj==0) g = 1;
					cout << T_in[ii*4+jj] + P*g - (e+P)*ev[ii]*ev[jj] << "   ";
				}
				cout << endl;
			}*/
			cout << "e: " << e << endl;
			cout << "u: " << endl;
			for (int ii = 0; ii < 4; ii++)
			{
				double tau_factor = ( ii == 3 ) ? 1.0/0.6 : 1.0;
				cout << tau_factor*ev[ii] << "   ";
			}
			cout << endl << endl;
			cout << "pi:" << endl;
			for (int ii = 0; ii < 4; ii++)
			{
				for (int jj = ii; jj < 4; jj++)
				{
					double tau_factor = 1.0;
					if (ii==3) tau_factor /= 0.6;
					if (jj==3) tau_factor /= 0.6;
					int g = -static_cast<int>( ii == jj );
					if (ii==0 && jj==0) g = 1;
					cout << tau_factor*(T_in[ii*4+jj] + P*g - (e+P)*ev[ii]*ev[jj]) << "   ";
				}
				//cout << endl;
			}
			cout << endl;
		}
      }
  }
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);

//if (true) exit(8);

	return;
}



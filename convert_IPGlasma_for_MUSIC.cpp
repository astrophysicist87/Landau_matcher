#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

constexpr bool use_full_Tmunu = false;	// false means zero tau-eta, x-eta, and y-eta components
										// of T^{\mu\nu} before doing Landau matching
int input_mode;
//constexpr int input_mode = 1;			// 0 - assumes reading in T^{\mu\nu}
										// 1 - assumes reading in e, u^\mu, and \pi^{\mu\nu}

const double hbarc = 0.1973269718;
//const double tau0 = 0.6;
double tau0;

void do_Landau_matching(
		const double xc, const double yc,
		const double T00_in, const double T0x_in, const double T0y_in, const double T0z_in,
		const double Txx_in, const double Txy_in, const double Txz_in,
		const double Tyy_in, const double Tyz_in, const double Tzz_in );

int main(int argc, char *argv[])
{
	if (argc < 2) exit(8);

	// read path to input file from command line
	string path_to_file = string(argv[1]);

	input_mode = ( argc > 2 ) ? stoi( argv[2] ) : 0;	// input_mode == 0 is default
	tau0       = ( argc > 3 ) ? stod( argv[3] ) : 0.6;	// tau0 == 0.6 fm/c is default

cerr << "argc = " << argc << endl;
cerr << "input_mode = " << input_mode << endl;
cerr << "tau0 = " << tau0 << endl;

	int xsize = -1, ysize = -1;
	double dx = 0.0, dy = 0.0;
	
	// then read in file itself
	ifstream infile(path_to_file.c_str());
	if (infile.is_open())
	{
		string line;
		int lineCount = 0;
		while ( getline (infile, line) )
		{
			// skip header line; just process and print it back out
			if ( lineCount++ < 1 )
			{
				string dummy;
				istringstream issHeader(line);
				issHeader >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
                          >> xsize >> dummy >> ysize >> dummy >> dummy >> dummy
                          >> dx >> dummy >> dy;
				cout << line << endl;
				continue;
			}
			// skip empty lines; don't print them back out since they're not
			// included in IP-Glasma default output
			if ( line.empty() )
			{
				//cout << line << endl;
				continue;
			}

			int ix, iy;
			double T00, Txx, Tyy, Tzz, T0x, T0y, T0z, Txy, Tyz, Txz;

			istringstream iss(line);
			if ( input_mode == 0 )
			{
				iss >> ix >> iy >> T00 >> Txx >> Tyy >> Tzz
					>> T0x >> T0y >> T0z >> Txy >> Tyz >> Txz;

				T0x *= -1.0;
				T0y *= -1.0;
				T0z *= -1.0;
				Txy *= -1.0;
				Tyz *= -1.0;
				Txz *= -1.0;
			}
			else
			{
				double dummy;
				double e, utau, ux, uy, ueta;
				double pi00, pi0x, pi0y, pi0eta, pixx, pixy, pixeta, piyy, piyeta, pietaeta;
				double g00 = 1.0, gxx = -1.0, gyy = -1.0, gzz = -1.0;
				iss >> dummy >> ix >> iy >> e >> utau >> ux >> uy >> ueta
					>> pi00 >> pi0x >> pi0y >> pi0eta >> pixx >> pixy
					>> pixeta >> piyy >> piyeta >> pietaeta;

				double uz = tau0 * ueta;
				T00 = (4./3.)*e*utau*utau - e*g00/3.0 + hbarc*pi00;
				T0x = (4./3.)*e*utau*ux               + hbarc*pi0x;
				T0y = (4./3.)*e*utau*uy               + hbarc*pi0y;
				T0z = (4./3.)*e*utau*uz               + hbarc*tau0*pi0eta;
				Txx = (4./3.)*e*ux*ux     - e*gxx/3.0 + hbarc*pixx;
				Txy = (4./3.)*e*ux*uy                 + hbarc*pixy;
				Txz = (4./3.)*e*ux*uz                 + hbarc*tau0*pixeta;
				Tyy = (4./3.)*e*uy*uy     - e*gyy/3.0 + hbarc*piyy;
				Tyz = (4./3.)*e*uy*uz                 + hbarc*tau0*piyeta;
				Tzz = (4./3.)*e*uz*uz     - e*gzz/3.0 + hbarc*tau0*tau0*pietaeta;
			}
			
			double xc = dx*(ix-0.5*xsize);
			double yc = dy*(iy-0.5*ysize);

			if ( use_full_Tmunu )
				do_Landau_matching( xc, yc, T00, T0x, T0y, T0z, Txx, Txy, Txz, Tyy, Tyz, Tzz );
			else
				do_Landau_matching( xc, yc, T00, T0x, T0y, 0.0, Txx, Txy, 0.0, Tyy, 0.0, Tzz );

		}
	}
	
	infile.close();

	return (0);
}


void do_Landau_matching(
		const double xc, const double yc,
		const double T00_in, const double T0x_in, const double T0y_in, const double T0z_in,
		const double Txx_in, const double Txy_in, const double Txz_in,
		const double Tyy_in, const double Tyz_in, const double Tzz_in )
{
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
		double e, PPlusPi;
		double ev[4], u[4];
		double pi[16];
		
		for (i = 0; i < 4; i++)
		{
			gsl_complex eval_i = gsl_vector_complex_get (eval, i);
			gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
			
			for (j = 0; j < 4; ++j)
			{
				gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
				ev[j] = GSL_REAL(z);
			}
			if ( ev[0]*ev[0] - ev[1]*ev[1] - ev[2]*ev[2] - ev[3]*ev[3] > 0.0 )
			{
				const double norm = sqrt(abs(ev[0]*ev[0] - ev[1]*ev[1] - ev[2]*ev[2] - ev[3]*ev[3]));
				for (int iev = 0; iev < 4; iev++) ev[iev] /= norm;
				
				if (ev[0] < 0.0) for (int iev = 0; iev < 4; iev++) ev[iev] *= -1.0;
				
				const double trace = T00_in - Txx_in - Tyy_in - Tzz_in;
				e = abs(GSL_REAL(eval_i));
				PPlusPi = (e - trace)/3.0;
				//cout << "e: " << e << endl;
				//cout << "u: " << endl;
				for (int ii = 0; ii < 4; ii++)
				{
					//double tau_factor = ( ii == 3 ) ? 1.0/tau0 : 1.0;
					u[ii] = ev[ii];
					//cout << tau_factor*ev[ii] << "   ";
				}
				//cout << endl << endl;
				//cout << "pi:" << endl;
				for (int ii = 0; ii < 4; ii++)
				for (int jj = ii; jj < 4; jj++)
				{
					/*double tau_factor = 1.0;
					if (ii==3) tau_factor /= tau0;
					if (jj==3) tau_factor /= tau0;*/
					int g = -static_cast<int>( ii == jj );
					if (ii==0 && jj==0) g = 1;
					//cout << tau_factor*(T_in[ii*4+jj] + P*g - (e+P)*ev[ii]*ev[jj]) << "   ";
					pi[ii*4+jj] = T_in[ii*4+jj] + PPlusPi*g - (e+PPlusPi)*ev[ii]*ev[jj];
				}
				//cout << endl;

				break;
			}
		}

		// units are: x, y (fm), e (GeV/fm^3), pi (1/fm^4)
		cout << 0 << "   "
             << xc << "   "
			 << yc << "   "
			 << e << "   "
			 //<< PPlusPi << "   "
			 << u[0] << "   "
			 << u[1] << "   "
			 << u[2] << "   "
			 << u[3] << "   ";
		for (int ii = 0; ii < 4; ii++)
		for (int jj = ii; jj < 4; jj++)
		{
			double tau_factor = 1.0 / hbarc;
			if (ii==3) tau_factor /= tau0;
			if (jj==3) tau_factor /= tau0;
			cout << pi[ii*4+jj]*tau_factor << "   ";
		}
		cout << endl;
	}

	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);

	return;
}



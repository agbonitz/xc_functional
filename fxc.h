#include <math.h>


// Computes the reduced temperature theta = T_Ha / T_Fermi from an absolute temperature in Hartree
double theta( double T_Ha, double rs, double xi )
{
	double pi = 3.141592653589793238462643383279502884197169399375105820974;
	double n,n_up,kF_sq;
	
	// compute the total density
	n = 1.0 / ( 4.0 * pi / 3.0 *pow( rs, 3.0 ) );
	
	// compute the density of the spin-up electrons (n_up >= n_down by definition)
	n_up = 0.50 * n * (1.0 + xi);
	
	// compute the square of the Fermi wave vector of the spin-up electrons
	kF_sq = pow(6.0*pi*pi*n_up, 2.0/3.0);
	
	// return the reduced temperature
	return T_Ha * 2.0 / kF_sq;
}


double fxc( double rs, double theta, double xi )
{
	// Fit parameters for the unpolarized case:
	double pol_b1, pol_b2, pol_b3, pol_b4, pol_b5;
	double pol_c1, pol_c2;
	double pol_d1, pol_d2, pol_d3, pol_d4, pol_d5;
	double pol_e1, pol_e2, pol_e3, pol_e4, pol_e5;
	
	
	// Fit parameters for the unpolarized case:
	double upol_b1, upol_b2, upol_b3, upol_b4, upol_b5;
	double upol_c1, upol_c2;
	double upol_d1, upol_d2, upol_d3, upol_d4, upol_d5;
	double upol_e1, upol_e2, upol_e3, upol_e4, upol_e5;
	
	// Fit parameters for the spin-interpolation function 
	double h1,h2,lambda1,lambda2;
	
	
	double h_fun, lambda_fun, alpha_fun, phi_fun;
	double pol_b, upol_b, pol_d, upol_d, pol_e, upol_e, pol_c, upol_c, pol_a, upol_a;
	double fxc0, fxc1;
	
	
	// Reduced temperatures for the spin-polarized and unpolarized case
	double pol_theta, upol_theta;
	
	double lam, pi;
	pi = 3.141592653589793238462643383279502884197169399375105820974;
	lam = pow( 4.0/(pi*9.0), 1.0/3.0 );
	
	
	// spin-interpolation parameters (T=0)
	h1 = 3.187472580;
	h2 = 7.746628020;
	
	// spin-interpolation parameters (finite-T)
	lambda1 = 1.85909536;
	lambda2 = 0.0;
	
	
	// ground-state parameters (xi=0):
	upol_b1 = 0.34369020;
	upol_c1 = 0.87594420;
	upol_d1 = 0.727008760;
	upol_e1 = 0.253882140;
	
	
	// finite-T parameters (xi=0):
	upol_b2 = 7.821595313560;
	upol_b3 = 0.3004839866620;
	upol_b4 = 15.84434671250;
	upol_b5 = upol_b3*sqrt(1.50) / lam; // this relation ensures the Debye-Hueckel limit for large T
	upol_c2 = -0.2301308435510;
	upol_d2 = 2.382647341440;
	upol_d3 = 0.302212372510;
	upol_d4 = 4.393477183950;
	upol_d5 = 0.7299513398450;
	upol_e2 = 0.8157951385990;
	upol_e3 = 0.06468444104810;
	upol_e4 = 15.09846204770;
	upol_e5 = 0.2307613574740; 
	
	
	// ground-state parameters (xi=1):
	pol_b1 = 0.849877040;
	pol_c1 = 0.911268730;
	pol_d1 = 1.486587180;
	pol_e1 = 0.274540970;
	
	
	// finite-T parameters (xi=1):
	pol_b2 = 3.040330120730;
	pol_b3 = 0.07757301312480;
	pol_b4 = 7.577035924890;
	pol_b5 = pol_b3 * pow( 2.0, 1.0/3.0 ) * sqrt(1.50) / lam; // this relation ensures the Debye-Hueckel limit for large T
	pol_c2 = -0.0307957123308;
	pol_d2 = 4.926849055110;
	pol_d3 = 0.08493872251790;
	pol_d4 = 8.32698211880;
	pol_d5 = 0.2188649521260;
	pol_e2 = 0.4009948565550;
	pol_e3 = 2.887731949620;
	pol_e4 = 6.334992370920;
	pol_e5 = 24.8230087530;
	
	
	
	// Compute the reduced temperatures
	upol_theta = theta * pow( 1.0+xi, 2.0/3.0 ); 
	pol_theta = upol_theta * pow( 2.0, -2.0/3.0 );
	
	
	// Evaluate the spin-interpolation function
	h_fun = ( 2.0/3.0 + h1*rs ) / ( 1.0 + h2*rs );
	lambda_fun = lambda1 + lambda2*sqrt(rs)*upol_theta;
	alpha_fun = 2.0 - h_fun / exp( upol_theta*lambda_fun ); 
	phi_fun = ( pow( 1.0+xi, alpha_fun ) + pow( 1.0-xi, alpha_fun ) - 2.0 ) / ( pow( 2.0, alpha_fun ) - 2.0 );
	
	
	
	
	// Compute the boundary values, fxc0 and fxc1
	
	pol_a = 0.610887*tanh(1.0/pol_theta) * ( 0.75 + 3.04363*pow( pol_theta, 2 ) + 1.7035*pow( pol_theta, 4 )  - 0.09227*pow( pol_theta, 3 ) ) / (  1.0 + 8.31051*pow( pol_theta, 2 ) + 5.1105*pow( pol_theta, 4 ) );
	upol_a = 0.610887*tanh(1.0/upol_theta) * ( 0.75 + 3.04363*pow( upol_theta, 2 ) - 0.09227*pow( upol_theta, 3 ) + 1.7035*pow( upol_theta, 4 ) ) / (  1.0 + 8.31051*pow( upol_theta, 2 ) + 5.1105*pow( upol_theta, 4 ) );
	
	
	pol_b = tanh( 1.0 / sqrt( pol_theta ) ) * ( pol_b1 + pol_b2*pow( pol_theta, 2 ) + pol_b3*pow( pol_theta, 4 ) ) / ( 1.0 + pol_b4*pow( pol_theta, 2 ) + pol_b5*pow( pol_theta, 4 ) );
	upol_b = tanh( 1.0 / sqrt( upol_theta ) ) * ( upol_b1 + upol_b2*pow( upol_theta, 2 ) + upol_b3*pow( upol_theta, 4 ) ) / ( 1.0 + upol_b4*pow( upol_theta, 2 ) + upol_b5*pow( upol_theta, 4 ) );
	
	
	pol_d = tanh( 1.0 / sqrt( pol_theta ) ) * ( pol_d1 + pol_d2*pow( pol_theta, 2 ) + pol_d3*pow( pol_theta, 4 ) ) / ( 1.0 + pol_d4*pow( pol_theta, 2 ) + pol_d5*pow( pol_theta, 4 ) );
	upol_d = tanh( 1.0 / sqrt( upol_theta ) ) * ( upol_d1 + upol_d2*pow( upol_theta, 2 ) + upol_d3*pow( upol_theta, 4 ) ) / ( 1.0 + upol_d4*pow( upol_theta, 2 ) + upol_d5*pow( upol_theta, 4 ) );
	
	
	pol_e = tanh( 1.0 / pol_theta ) * ( pol_e1 + pol_e2*pow( pol_theta, 2 ) + pol_e3*pow( pol_theta, 4 ) ) / ( 1.0 + pol_e4*pow( pol_theta, 2 ) + pol_e5*pow( pol_theta, 4 ) );
	upol_e = tanh( 1.0 / upol_theta ) * ( upol_e1 + upol_e2*pow( upol_theta, 2 ) + upol_e3*pow( upol_theta, 4 ) ) / ( 1.0 + upol_e4*pow( upol_theta, 2 ) + upol_e5*pow( upol_theta, 4 ) );
	
	
	pol_c = pol_e * ( pol_c1 + pol_c2 / exp( 1.0 / pol_theta ) );
	upol_c = upol_e * ( upol_c1 + upol_c2 / exp( 1.0 / upol_theta ) );
	
	
	fxc0 = -1.0/rs * ( upol_a + sqrt(rs)*upol_b + rs*upol_c ) / ( 1.0 + sqrt(rs)*upol_d + rs*upol_e );
	fxc1 = -1.0/rs * ( pol_a*pow( 2.0, 1.0/3.0 ) + sqrt(rs)*pol_b + rs*pol_c ) / ( 1.0 + sqrt(rs)*pol_d + rs*pol_e );
	
	
	return fxc0 + phi_fun*( fxc1 - fxc0 );

}
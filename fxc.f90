! ####################################
! ### Parametrization of f_xc ########
! ####################################


! Computes the reduced temperature theta = T_Ha / T_Fermi from an absolute temperature in Hartree
DOUBLE PRECISION function theta(T_Ha,rs,xi) 
	IMPLICIT NONE 
	DOUBLE PRECISION, intent(in) :: T_Ha,rs,xi
	
	DOUBLE PRECISION :: pi
	DOUBLE PRECISION :: n, n_up, kF_sq
	
	pi = 3.141592653589793238462643383279502884197169399375105820974d0
	
	! compute the total density
	n = 1.0d0 / ( 4.0d0 * pi / 3.0d0 *rs**3.0d0 )
	
	! compute the density of the spin-up electrons (n_up >= n_down by definition)
	n_up = 0.50d0 * n * (1.0d0 + xi)
	
	! compute the square of the Fermi wave vector of the spin-up electrons
	kF_sq = (6.0d0*pi*pi*n_up)**(2.0d0/3.0d0)
	
	! return the reduced temperature
	theta = T_Ha * 2.0d0 / kF_sq
end function theta


! Return the XC free energy of the UEG
DOUBLE PRECISION function fxc(rs,theta,xi)
	IMPLICIT NONE
	DOUBLE PRECISION, intent(in) :: rs,theta,xi
	
	! Fit parameters for the polarized case:
	DOUBLE PRECISION :: pol_b1, pol_b2, pol_b3, pol_b4, pol_b5
	DOUBLE PRECISION :: pol_c1, pol_c2
	DOUBLE PRECISION :: pol_d1, pol_d2, pol_d3, pol_d4, pol_d5
	DOUBLE PRECISION :: pol_e1, pol_e2, pol_e3, pol_e4, pol_e5
	
	
	! Fit parameters for the unpolarized case:
	DOUBLE PRECISION :: upol_b1, upol_b2, upol_b3, upol_b4, upol_b5
	DOUBLE PRECISION :: upol_c1, upol_c2
	DOUBLE PRECISION :: upol_d1, upol_d2, upol_d3, upol_d4, upol_d5
	DOUBLE PRECISION :: upol_e1, upol_e2, upol_e3, upol_e4, upol_e5
	
	! Fit parameters for the spin-interpolation function 
	DOUBLE PRECISION :: h1,h2,lambda1,lambda2
	
	
	DOUBLE PRECISION :: h_fun, lambda_fun, alpha_fun, phi_fun
	DOUBLE PRECISION :: pol_b, upol_b, pol_d, upol_d, pol_e, upol_e, pol_c, upol_c, pol_a, upol_a
	DOUBLE PRECISION :: fxc0, fxc1
	
	
	! Reduced temperatures for the spin-polarized and unpolarized case
	DOUBLE PRECISION :: pol_theta, upol_theta
	
	DOUBLE PRECISION :: lam, pi 
	pi = 3.141592653589793238462643383279502884197169399375105820974d0
	lam = (4.0d0/(pi*9.0d0))**(1.0d0/3.0d0)
	
	
	! spin-interpolation parameters (T=0)
	h1 = 3.187472580d0
	h2 = 7.746628020d0
	
	! spin-interpolation parameters (finite-T)
	lambda1 = 1.85909536 
	lambda2 = 0.0d0
	
	! ground-state parameters (xi=0):
	upol_b1 = 0.34369020d0
	upol_c1 = 0.87594420d0
	upol_d1 = 0.727008760d0
	upol_e1 = 0.253882140d0
	
	
	! finite-T parameters (xi=0):
	upol_b2 = 7.821595313560d0
	upol_b3 = 0.3004839866620d0
	upol_b4 = 15.84434671250d0
	upol_b5 = upol_b3*sqrt(1.50d0) / lam ! this relation ensures the Debye-Hueckel limit for large T
	upol_c2 = -0.2301308435510d0
	upol_d2 = 2.382647341440d0
	upol_d3 = 0.302212372510d0
	upol_d4 = 4.393477183950d0
	upol_d5 = 0.7299513398450d0
	upol_e2 = 0.8157951385990d0
	upol_e3 = 0.06468444104810d0
	upol_e4 = 15.09846204770d0
	upol_e5 = 0.2307613574740d0 
		
	
	! ground-state parameters (xi=1):
	pol_b1 = 0.849877040d0
	pol_c1 = 0.911268730d0
	pol_d1 = 1.486587180d0
	pol_e1 = 0.274540970d0
	
	! finite-T parameters (xi=1):
	pol_b2 = 3.040330120730d0
	pol_b3 = 0.07757301312480d0
	pol_b4 = 7.577035924890d0
	pol_b5 = pol_b3 * 2.0d0**(1.0d0/3.0d0) * sqrt(1.50d0) / lam ! this relation ensures the Debye-Hueckel limit for large T
	pol_c2 = -0.03079571233080d0
	pol_d2 = 4.926849055110d0
	pol_d3 = 0.08493872251790d0
	pol_d4 = 8.32698211880d0
	pol_d5 = 0.2188649521260d0
	pol_e2 = 0.4009948565550d0
	pol_e3 = 2.887731949620d0
	pol_e4 = 6.334992370920d0
	pol_e5 = 24.8230087530d0
	
	
	
	! Compute the reduced temperatures
	upol_theta = theta * ( 1.0d0+xi )**(2.0d0/3.0d0) 
	pol_theta = upol_theta * 2.0d0**( -2.0d0/3.0d0 )
	
	
	! Evaluate the spin-interpolation function
	h_fun = ( 2.0d0/3.0d0 + h1*rs ) / ( 1.0d0 + h2*rs )
	lambda_fun = lambda1 + lambda2*sqrt(rs)*upol_theta 
	alpha_fun = 2.0d0 - h_fun / exp( upol_theta*lambda_fun ) 
	phi_fun = ( (1.0d0+xi)**alpha_fun + (1.0d0-xi)**alpha_fun - 2.0d0 ) / ( 2.0d0**alpha_fun - 2.0d0 )
	
	
	
	
	! Compute the boundary values, fxc0 and fxc1
	
	pol_a = 0.610887d0*tanh(1.0d0/pol_theta) * ( 0.75d0 + 3.04363d0*pol_theta**2 + 1.7035d0*pol_theta**4  &
     &  - 0.09227d0*pol_theta**3 ) &
     &  / (  1.0d0 + 8.31051d0*pol_theta**2 + 5.1105d0*pol_theta**4 )
	upol_a = 0.610887d0*tanh(1.0d0/upol_theta) * ( 0.75d0 + 3.04363d0*upol_theta**2 &
     &  - 0.09227d0*upol_theta**3 + 1.7035d0*upol_theta**4 ) &
     &  / (  1.0d0 + 8.31051d0*upol_theta**2 + 5.1105d0*upol_theta**4 )
     
     
	pol_b = tanh( 1.0d0 / sqrt( pol_theta ) ) * ( pol_b1 + pol_b2*pol_theta**2 + pol_b3*pol_theta**4 ) &
     &  / ( 1.0d0 + pol_b4*pol_theta**2 + pol_b5*pol_theta**4 )
        upol_b = tanh( 1.0d0 / sqrt( upol_theta ) ) * ( upol_b1 + upol_b2*upol_theta**2 + upol_b3*upol_theta**4 ) &
     &  / ( 1.0d0 + upol_b4*upol_theta**2 + upol_b5*upol_theta**4 )
     
     
        pol_d = tanh( 1.0d0 / sqrt( pol_theta ) ) * ( pol_d1 + pol_d2*pol_theta**2 + pol_d3*pol_theta**4 ) &
     &  / ( 1.0d0 + pol_d4*pol_theta**2 + pol_d5*pol_theta**4 )
	upol_d = tanh( 1.0d0 / sqrt( upol_theta ) ) * ( upol_d1 + upol_d2*upol_theta**2 + upol_d3*upol_theta**4 ) &
     &  / ( 1.0d0 + upol_d4*upol_theta**2 + upol_d5*upol_theta**4 )
     
     
	pol_e = tanh( 1.0d0 / pol_theta ) * ( pol_e1 + pol_e2*pol_theta**2 + pol_e3*pol_theta**4 ) &
     &  / ( 1.0d0 + pol_e4*pol_theta**2 + pol_e5*pol_theta**4 )
	upol_e = tanh( 1.0d0 / upol_theta ) * ( upol_e1 + upol_e2*upol_theta**2 + upol_e3*upol_theta**4 ) &
     &  / ( 1.0d0 + upol_e4*upol_theta**2 + upol_e5*upol_theta**4 )
! 	

	pol_c = pol_e * ( pol_c1 + pol_c2 / exp( 1.0d0 / pol_theta ) )
	upol_c = upol_e * ( upol_c1 + upol_c2 / exp( 1.0d0 / upol_theta ) )
	
	
	fxc0 = -1.0d0/rs * ( upol_a + sqrt(rs)*upol_b + rs*upol_c ) / ( 1.0d0 + sqrt(rs)*upol_d + rs*upol_e )
	fxc1 = -1.0d0/rs * ( pol_a*2.0d0**(1.0d0/3.0d0) + sqrt(rs)*pol_b + rs*pol_c ) / ( 1.0d0 + sqrt(rs)*pol_d + rs*pol_e )
	

	
	! Compute final result for fxc(rs,theta,xi):
	fxc = fxc0 + phi_fun*( fxc1 - fxc0 )
end function fxc
	
	
program fcn1
	implicit none
	DOUBLE PRECISION, external :: fxc
	DOUBLE PRECISION, external :: theta

	print *, "theta(2.25, 0.5, 0.3) = ", theta(2.250d0,0.50d0,0.30d0)
	print *, "fxc(0.2,0.3,0.3) = ", fxc(0.20d0,0.30d0,0.30d0)

end program fcn1


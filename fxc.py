import numpy as np

def a(t):
    return 0.610887*np.tanh(1.0/t) * (0.75 + 3.04363*t**2 - 0.09227*t**3 + 1.7035*t**4) / (1 + 8.31051*t**2 + 5.1105*t**4)

def b(t,xi,b1,b2,b3,b4): 
    if xi == 1:
        omega = 2**(1/3)
    else:
        omega = 1
    return np.tanh(1/np.sqrt(t)) * (b1 + b2*t**2 + b3*t**4) / (1 + b4*t**2 + b3*np.sqrt(3/2) *omega* (4/(9*np.pi))**(-1/3)*t**4)

def e(t,e1,e2,e3,e4,e5): 
    return np.tanh(1/t)* (e1 + e2*t**2 + e3*t**4) / (1 + e4*t**2 + e5*t**4)

def c(t,c1,c2,e1,e2,e3,e4,e5):
    return (c1 + c2*np.exp(-1/t)) * e(t,e1,e2,e3,e4,e5)

def d(t,d1,d2,d3,d4,d5):
    return np.tanh(1/np.sqrt(t)) * (d1 + d2*t**2 + d3*t**4) / (1 + d4*t**2 + d5*t**4)

def A(rs,t,xi,b1,b2,b3,b4,c1,c2,e1,e2,e3,e4,e5):
    if xi == 1:
        omega = 2**(1/3)
    else:
        omega = 1
    return omega*a(t) + b(t,xi,b1,b2,b3,b4)*rs**(1/2) + c(t,c1,c2,e1,e2,e3,e4,e5)*rs
    
def B(rs,t,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5):
    return -rs - d(t,d1,d2,d3,d4,d5)*rs**(3/2) - e(t,e1,e2,e3,e4,e5)*rs**2

def f_xc_xi(rs,t,xi,b1,b2,b3,b4,c1,c2,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5):
    return A(rs,t,xi,b1,b2,b3,b4,c1,c2,e1,e2,e3,e4,e5) / B(rs,t,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5)

def dA_drs(rs,t,xi,b1,b2,b3,b4,c1,c2,e1,e2,e3,e4,e5):
    return b(t,xi,b1,b2,b3,b4) / (2*rs**(1/2)) + c(t,c1,c2,e1,e2,e3,e4,e5)

def dB_drs(rs,t,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5):
    return -1 -3/2*d(t,d1,d2,d3,d4,d5)*rs**(1/2) - 2*e(t,e1,e2,e3,e4,e5)*rs

def dfxc_drs(rs,t,xi,b1,b2,b3,b4,c1,c2,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5):
    return 1/B(rs,t,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5) * dA_drs(rs,t,xi,b1,b2,b3,b4,c1,c2,e1,e2,e3,e4,e5) - A(rs,t,xi,b1,b2,b3,b4,c1,c2,e1,e2,e3,e4,e5) / (B(rs,t,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5) * B(rs,t,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5)) * dB_drs(rs,t,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5)

def Fit_fun_Vrs_xi(rst,xi,b1,b2,b3,b4,c1,c2,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5):
    rs,t = rst
    return 2*rs*f_xc_xi(rs,t,xi,b1,b2,b3,b4,c1,c2,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5) + rs**2*dfxc_drs(rs,t,xi,b1,b2,b3,b4,c1,c2,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5)

fitParZeroT_unpol=[ # b1,c1,d1,e1
0.3436902,0.8759442,0.72700876,0.25388214
]
b1_0,c1_0,d1_0,e1_0 = fitParZeroT_unpol

fitPar_unpol=[ #b2,b3,b4,c2,d2,d3,d4,d5,e2,e3,e4,e5
7.82159531356 ,
0.300483986662 ,
15.8443467125 ,
-0.230130843551 ,
2.38264734144 ,
0.30221237251 ,
4.39347718395 ,
0.729951339845 ,
0.815795138599 ,
0.0646844410481 ,
15.0984620477 ,
0.230761357474 ,
]
b2_0,b3_0,b4_0,c2_0,d2_0,d3_0,d4_0,d5_0,e2_0,e3_0,e4_0,e5_0 = fitPar_unpol

def f_xc_xi0(rs,t):
    return f_xc_xi(rs,t,0,b1_0,b2_0,b3_0,b4_0,c1_0,c2_0,d1_0,d2_0,d3_0,d4_0,d5_0,e1_0,e2_0,e3_0,e4_0,e5_0) 

def Fit_fun_Vrs_xi0(rs,t):
    rst = (rs,t)
    return Fit_fun_Vrs_xi(rst,0,b1_0,b2_0,b3_0,b4_0,c1_0,c2_0,d1_0,d2_0,d3_0,d4_0,d5_0,e1_0,e2_0,e3_0,e4_0,e5_0)

def dfxc_drs_xi0(rs,t):
    return dfxc_drs(rs,t,0,b1_0,b2_0,b3_0,b4_0,c1_0,c2_0,d1_0,d2_0,d3_0,d4_0,d5_0,e1_0,e2_0,e3_0,e4_0,e5_0) 


fitParZeroT_pol=[ # b1,c1,d1,e1
0.84987704,0.91126873,1.48658718,0.27454097
]
b1_1,c1_1,d1_1,e1_1 = fitParZeroT_pol


fitPar_pol=[ #b2,b3,b4,c2,d2,d3,d4,d5,e2,e3,e4,e5
3.04033012073 ,
0.0775730131248 ,
7.57703592489 ,
-0.0307957123308 ,
4.92684905511 ,
0.0849387225179 ,
8.3269821188 ,
0.218864952126 ,
0.400994856555 ,
2.88773194962 ,
6.33499237092 ,
24.823008753 ,
]
b2_1,b3_1,b4_1,c2_1,d2_1,d3_1,d4_1,d5_1,e2_1,e3_1,e4_1,e5_1 = fitPar_pol

def f_xc_xi1(rs,t):
    return f_xc_xi(rs,t,1,b1_1,b2_1,b3_1,b4_1,c1_1,c2_1,d1_1,d2_1,d3_1,d4_1,d5_1,e1_1,e2_1,e3_1,e4_1,e5_1) 

def Fit_fun_Vrs_xi1(rs,t):
    rst = (rs,t)
    return Fit_fun_Vrs_xi(rst,1,b1_1,b2_1,b3_1,b4_1,c1_1,c2_1,d1_1,d2_1,d3_1,d4_1,d5_1,e1_1,e2_1,e3_1,e4_1,e5_1)

def dfxc_drs_xi1(rs,t):
    return dfxc_drs(rs,t,1,b1_1,b2_1,b3_1,b4_1,c1_1,c2_1,d1_1,d2_1,d3_1,d4_1,d5_1,e1_1,e2_1,e3_1,e4_1,e5_1) 


### Spin interpolation
def g(rs,g1,g2,g3):
    ret = g1 + g2*rs
    ret /= 1 + g3*rs
    return ret 

def lam(rs,t,lam1,lam2):
    return lam1 + lam2*t*np.sqrt(rs)

def alpha(rs,t,g1,g2,g3,lam1,lam2):
    return 2.0 - g(rs,g1,g2,g3)*np.exp(-t*lam(rs,t,lam1,lam2))

def A_phi(rs,t,xi,g1,g2,g3,lam1,lam2):
    return (1+xi)**alpha(rs,t,g1,g2,g3,lam1,lam2) + (1-xi)**alpha(rs,t,g1,g2,g3,lam1,lam2) -2.0

def B_phi(rs,t,xi,g1,g2,g3,lam1,lam2):
    return 2.0**alpha(rs,t,g1,g2,g3,lam1,lam2) -2.0
    
def phi(rs,t,xi,g1,g2,g3,lam1,lam2):
    return A_phi(rs,t,xi,g1,g2,g3,lam1,lam2) / B_phi(rs,t,xi,g1,g2,g3,lam1,lam2)

def full_f_xc(rs,t,xi,g1,g2,g3,lam1,lam2):
    return f_xc_xi0(rs, t) + ( f_xc_xi1(rs,2.0**(-2.0/3.0)*t) - f_xc_xi0(rs, t) ) * phi(rs,t,xi,g1,g2,g3,lam1,lam2) 

def dA_phi_dalpha(rs,t,xi,g1,g2,g3,lam1,lam2):
    ret = (1+xi)**(alpha(rs,t,g1,g2,g3,lam1,lam2)) * np.log(1+xi)
    ret += (1-xi)**(alpha(rs,t,g1,g2,g3,lam1,lam2)) * np.log(1-xi+1e-15)
    return ret

def dB_phi_dalpha(rs,t,xi,g1,g2,g3,lam1,lam2):
    return 2.0**(alpha(rs,t,g1,g2,g3,lam1,lam2)) * np.log(2.0)

def dphi_dalpha(rs,t,xi,g1,g2,g3,lam1,lam2):
    ret = 1.0 / B_phi(rs,t,xi,g1,g2,g3,lam1,lam2) * dA_phi_dalpha(rs,t,xi,g1,g2,g3,lam1,lam2)
    ret -= A_phi(rs,t,xi,g1,g2,g3,lam1,lam2) * (B_phi(rs,t,xi,g1,g2,g3,lam1,lam2))**(-2.0) * dB_phi_dalpha(rs,t,xi,g1,g2,g3,lam1,lam2)
    return ret

def dg_drs(rs,t,g1,g2,g3):
    return g2 / (1.0+g3*rs) - g(rs,g1,g2,g3) / (1+g3*rs) * g3 

def dalhpa_drs(rs,t,xi,g1,g2,g3,lam1,lam2):
    ret = -dg_drs(rs,t,g1,g2,g3)*np.exp(-t*lam(rs,t,lam1,lam2)) 
    ret += g(rs,g1,g2,g3)*t*lam2*t*0.5*rs**(-1.0/2.0)*np.exp(-t*lam(rs,t,lam1,lam2))
    return ret

def dphi_drs(rs,t,xi,g1,g2,g3,lam1,lam2):
    return dphi_dalpha(rs,t,xi,g1,g2,g3,lam1,lam2) * dalhpa_drs(rs,t,xi,g1,g2,g3,lam1,lam2) 

def full_dfxc_drs(rs,t,xi,g1,g2,g3,lam1,lam2):
    ret = dfxc_drs_xi0(rs,t)
    ret += ( dfxc_drs_xi1(rs,2.0**(-2.0/3.0)*t) - dfxc_drs_xi0(rs,t) ) * phi(rs,t,xi,g1,g2,g3,lam1,lam2) 
    ret += ( f_xc_xi1(rs,2.0**(-2.0/3.0)*t) - f_xc_xi0(rs, t) ) * dphi_drs(rs,t,xi,g1,g2,g3,lam1,lam2)
    return ret

g1,g2,g3 = [2/3,3.18747258,7.74662802] 


def Fit_fun_Vrs_tot(rstxi,lam1):
    rs,t,xi = rstxi
    return 2*rs*full_f_xc(rs,t*(1+xi)**(2/3),xi,g1,g2,g3,lam1,0) + rs**2*full_dfxc_drs(rs,t*(1+xi)**(2/3),xi,g1,g2,g3,lam1,0)

fitPar_Phi_finiteT=[1.85909536]

def Vrs_tot(rs,t,xi):
    return Fit_fun_Vrs_tot((rs,t,xi),fitPar_Phi_finiteT[0])


# this function returns the xc free energy for given density parameter rs, 
# reduced spin up temperatuer t=T/T_{F\uparrow} and spin polarization \xi=(n_up-n_down)/ (n_up+n_down)
def fxc_tot(rs,t,xi):
    return full_f_xc(rs,t*(1+xi)**(2/3),xi,g1,g2,g3,fitPar_Phi_finiteT[0],0)

def n(rs):
	return 1 / (4*np.pi / 3 * rs**3);

def nup(xi,n):
	return (1 + xi)/2 * n;
    
def t(T,nup):
	return T / ((6*np.pi*np.pi*nup)**(2/3)) * 2.0;

# this function returns the xc free energy per particle for given density parameter rs, 
# reduced temperatuer T in Hartree and spin polarization \xi=(n_up-n_down)/ (n_up+n_down)
def fxc_tot_T(rs,T,xi):
    return fxc_tot(rs, t(T,nup(xi,n(rs))), xi)
    
#print(fxc_tot_T(0.5, 2.25, 0.3))
#print(fxc_tot(1.5,1,0),'blub')


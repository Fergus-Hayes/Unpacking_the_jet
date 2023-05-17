import numpy as np

c, H_0 = 3.e5, 67.8
Omega_M, Omega_D = 0.308, 0.692

z_max = 2.6

gamma = 0.55

R_g = 1.e-6
R_GW = 0.9
eps_NSBH = 0.2
eps_BNS = 1.

R_BNS, R_NSBH = 1.e-6, 1.e-7

T_0 = 9.8
N_z, N = 35., 107

f_r, f_FOV = N_z/N, 0.1

l10_L_0 = 51.5
alpha, beta = 0.6, 2.4
l10_min,l10_max = 50.5, 56.#l10_L_0 - 3., l10_L_0 + 3.

l10_L_02 = np.log10(8.55e+51)
gamma2 = 2.4

A = 1.
E_p = 800.
E_min, E_max = 1., 10000.
E_alpha, E_beta = -1., -2.25

e_min, e_max = 15., 150.

F_min = 2.8e-8#5.e-9#2.8e-8

SF_BAT = (0.47, -0.05, 1.46, 1.45, 0.6e-7, 5.5e-9)

model = 'PL'
theta_c =10.3*np.pi/180.
theta_j = 30.*np.pi/180.
s = 8.2

fr_a, fr_b, fr_c, fr_d = 0.017, 0.13, 3.3, 5.3

t_m=20.e6*(60*60*24*365)/3.086e19

cosmo = (c,H_0,Omega_M,Omega_D)
lf_paras = (np.power(10.,l10_L_0), alpha, beta,
           np.power(10.,l10_min), np.power(10.,l10_max))
slf_paras = (np.power(10.,l10_L_0), gamma, np.power(10.,l10_min))
spec_paras = (A, E_p, E_min, E_max, E_alpha, E_beta)
fr_paras = (fr_a, fr_b, fr_c, fr_d)
det_paras = (e_min, e_max)

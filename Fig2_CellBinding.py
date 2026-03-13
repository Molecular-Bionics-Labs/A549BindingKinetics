import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import scipy.optimize as sci
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from scipy.optimize import curve_fit

import _pat_pmpc as pat 
import constant as const

"""
This script has the fitting for the PMPC Psome binidng with a549 cells.
binding_odeint() is a function of bining. 
get_odeint_solution() is a function of resolution we can get form the binding. 
You need to input the experiment time (t_exp), the data to fit (y_exp). 
Then, use curve_fit(fun_name, t_exp, y_exp, bounds=boundary) to fit the data. 
You can also choose other algrithm to fit your parameters. 
c_mic = 2.35*10**-7 ; c_disk = 5.14*10**-7
"""

# 1. Loading data. 
fp_read = "./Data/A549_15Cell.xlsx"
NA = 6.022*10**23  
A549_CellType = np.array([61.2, 0.1, 334.1])

dsh = const.dsh_A549


def binding_odeint(ft, t, nu, kon, koff, ke):
    """
    Here is the binding function which can be solved by the 
    scipy.integrate.odeint. 

    Parameters: 
    ----------
    ft: solution for the function. 

    t: time. 
    
    nu: valence. 

    kon: on rate. 

    koff: off rate.

    ke: endocytosis rate. 
    """
    drhoPdt = -kon*ft[0]*ft[1]**nu + koff * ft[2]
    drhoRdt = -kon*ft[0]*ft[1]**nu + koff * ft[2]
    drhoPRdt = kon*ft[0]*ft[1]**nu - koff * ft[2] - ke*ft[2]

    drhoEdt = ke*ft[2]

    return [drhoPdt, drhoRdt, drhoPRdt, drhoEdt]

def binding_odeint_tube(ft, t, nu, kon, koff, ke):
    """
    Here is the binding function which can be solved by the 
    scipy.integrate.odeint. 

    Parameters: 
    ----------
    ft: solution for the function. 

    t: time. 

    kon: on rate. 

    koff: off rate.

    ke: endocytosis rate. 
    """
    kon = kon/rho0_tube**nu

    drhoPdt = -kon*ft[0]*ft[1]**nu + koff * ft[2]
    drhoRdt = -kon*ft[0]*ft[1]**nu + koff * ft[2]
    drhoPRdt = kon*ft[0]*ft[1]**nu - koff * ft[2] - ke*ft[2]

    drhoEdt = ke*ft[2]
    
    return [drhoPdt, drhoRdt, drhoPRdt, drhoEdt]

def get_odeint_solution(t, rhoR, nu, kon1, koff1, ke1, kon2, ke2, koff2):
    """
    Here is the fiting function. 

    Parameters:
    ----------
    fun_name: function name. 

    t_exp: experiment time point. 

    kon: on rate. 

    koff: off rate 

    ke: endocytosis rate. 
    """
    
    """Psome.""" 
    rhoPR0, rhoE0 = 0, 0 
    rhoR_Psome = rhoR/nu
    y_init = [rho0_Psome, rhoR_Psome, rhoPR0, rhoE0]
    args1 = (nu, kon1/rhoR_Psome**nu, koff1, ke1)
    # 2. Solve function and get the result. 
    sol1 = odeint(binding_odeint, y_init, t_exp_Psome, args=args1)
    y1 = sol1[:, 2]/np.max(sol1[:2])

    """Disk."""
    nu2 = nu/(2*45*Rg / 15**2)
    rhoR_disk = rhoR / nu2 
    y_init = [rho0_disk, rhoR_disk, rhoPR0, rhoE0]
    args2 = (nu2, kon2/rhoR_disk**nu2, koff2, ke2)
    # 2. Solve function and get the result. 
    sol2 = odeint(binding_odeint, y_init, t_exp_disk, args=args2)
    y2 = sol2[:, 2]/np.max(sol2[:, 2])
    y = np.hstack((y1, y2))
    
    return y

def get_odeint_solution_tube(t, nu, kon, koff, ke):
    # 1. Initial values. 
    rhoPR0, rhoE0 = 0, 0 
    y_init = [rho0_tube, rhoR_tube, rhoPR0, rhoE0]
    args1 = (nu, kon, koff, ke)

    # 2. Solve function and get the result. 
    sol1 = odeint(binding_odeint_tube, y_init, t_exp_tube_fit, args=args1)
    y = sol1[:, 2]/np.max(sol1[:2]) * 0.5

    return y


def pltfig(texp, yexp, yerror, tsim, ysim, colori, par_type, dashx=None, 
           dashy=None):
    plt.figure(figsize=(10,8))
    plt.errorbar(texp, yexp, yerr=yerror,     
                 fmt='o', # 'o' for scatter, no lines
                 markersize=8,  # Adjust scatter point size
                 capsize=5,     # Caps on error bars
                 capthick=2,    # Thickness of cap lines
                 elinewidth=2,  # Thickness of error bar lines
                 color=colori,   # Set color
                 alpha=0.8  # Transparency
                 )
    plt.grid(axis='y', linestyle='-', linewidth=0.7, alpha=0.7)
    plt.plot(tsim, ysim, color=colori)
    # Micelle 
    # plt.plot(t_sim[:4000], y_sim[:4000], color=color_i)
    # plt.plot(t_sim[4000:6000], y_sim[4000:6000], color=color_i, linestyle="--")
    fp_save = "./Figure/Fig2{0}.svg".format(par_type)
    plt.xlim([0, 60])
    plt.ylim([0, 1])
    if dashx is None:
        pass 
    else: 
        plt.plot(dashx, dashy, linestyle="--")
    plt.savefig(fp_save, dpi=800)

ves = 0
Rg = 5.5 # The length of the PMPC ligand. 
# Density concentraiton of the NPs
rhoR0 = 0
c_Psome = 1.43*10**-8
c_disk = 5.14*10**-7
c_tube = 1.63*10**-9    
rho0_Psome = c_Psome * NA * 10 **-9 / 15  # turn to NP per cell. 
rho0_tube = c_tube * NA * 10 **-9 / 15  # turn to NP per cell. 
rho0_disk = c_disk * NA * 10 **-9 / 15  # turn to NP per cell. 
# Experiment data and data for fitting.  
data_Psome = pd.read_excel(fp_read, sheet_name="Psome")
data_disk = pd.read_excel(fp_read, sheet_name="Disk")
data_tube = pd.read_excel(fp_read, sheet_name="Tube")

t_exp_Psome = data_Psome.values[:, 0]
y_Psome = data_Psome.values[:, 1]
error_Psome = data_Psome.values[:, 2]
y_max_Psome = np.max(y_Psome)
theta_Psome = y_Psome/y_max_Psome 

t_exp_disk = data_disk.values[:30, 0]

y_disk = data_disk.values[:30, 1]
error_disk = data_disk.values[:30, 2]
y_max_disk = np.max(y_disk)
theta_disk = y_disk/y_max_disk

t_exp_tube_fit = data_tube.values[:-14, 0]
y_tube_fit = data_tube.values[:-14, 1]
y_max_tube = np.max(y_tube_fit)
theta_tube_fit = y_tube_fit/y_max_tube

"""Fitting"""
fitting_function = "curve_fit"
t_exp = np.hstack((t_exp_Psome, t_exp_disk))
theta = np.hstack((theta_Psome, theta_disk))
colors = [[0.8, 0.1, 0.1], [0, 0.392, 0], [0.804, 0.522, 0.247]]
# rho_R, nu_psome, kon_psome, koff_psome, ke_psome, kon_disk, kedis,koff_disk, 
boundary = np.array([[1.25*10**4, 2.5,  0.05*10**-1, 10**-3, 8*10**-2, 2, 10**-3, 1.2*10**-2],
                     [25*10**4, 3,  0.1*10**-1, 10**-1, 9*10**-2,   2.1, 10**-1, 1.3*10**-2]]) 

if fitting_function=="curve_fit":

    popt, *pcov = curve_fit(get_odeint_solution, t_exp, theta, bounds=boundary)
    print("Orignal fitting results.", "\n", popt)
    
    R_Psome = 40
    R_disk = 16
    """
    Posome fit:
    """
    nu_Psome = popt[1]
    rhoR0 = popt[0]
    rhoR_Psome = popt[0]/nu_Psome
    popt_Psome = (popt[1], popt[2]/rhoR_Psome**nu_Psome, popt[3],popt[4])

    y_init = [rho0_Psome, rhoR_Psome, 0, 0]
    t_sim = np.linspace(0, 60, 6001)
    y_sol = odeint(binding_odeint, y_init, t_sim, args=popt_Psome)
    print("Psome Occupied Receptor", np.max(y_sol[:, 2])/rhoR_Psome)

    y_sol = y_sol[:, 2]/ np.max(y_sol[:,2]) * np.max(y_Psome)
    pltfig(t_exp_Psome, y_Psome, error_Psome, t_sim, y_sol, colors[0], "Psome")

    # Use PAT to modify the koff.<----------------------------------(pat).

    Psome_rec = nu_Psome*A549_CellType/ np.sum(A549_CellType)
    ka_Psome = pat.pmpc_ka(Psome_rec, dsh, R_Psome, 25, phi=1)
    kon_Psome = popt[2]/(rhoR_Psome*15/(10**-9*NA))**nu_Psome
    koff_Psome = kon_Psome / ka_Psome
    print("nu_Psomes=", "{:.3e}".format(popt[1]), 
          "; ka_Psomes=", "{:.3e}".format(ka_Psome), 
          "; kon_Psome=", "{:.3e}".format(kon_Psome/60), 
          "; koff_Psome=", "{:.3e}".format(koff_Psome/60),
          "; ke_Psome=", "{:.3e}".format(popt[4]/60),
          "; residence time", "{:.3e}".format(60/koff_Psome),
          "receptor density", "{:.3e}".format(rhoR0*nu_Psome),
          "Nanoparticle left", "{:.3e}".format(rho0_Psome/rhoR0)
          )

    """
    Disk.
    """
    t_exp_disk = data_disk.values[:, 0]
    y_disk = data_disk.values[:, 1]
    error_disk = data_disk.values[:, 2]
    y_max_disk = np.max(y_disk)
    theta_disk = y_disk/y_max_disk
    nu_disk = nu_Psome/((2*R_Psome*Rg- Rg**2)/ (R_disk)**2) # The geometry. 
    rhoR_disk = popt[0]/nu_disk
    popt_disk = (nu_disk, popt[5]/rho0_disk**nu_disk, popt[6], popt[7])

    y_init = [rho0_disk, rhoR_disk, 0, 0]
    y_sol = odeint(binding_odeint, y_init, t_sim, args=popt_disk)
    print("Disk Occupied Receptor", np.max(y_sol[:,2])/rhoR_disk)

    y_sol = y_sol[:, 2]/ np.max(y_sol[:,2]) * 0.84
    pltfig(t_exp_disk, y_disk, error_disk, t_sim[0:4000], y_sol[0:4000], 
           colors[1], "Disk", dashx=t_sim[4000:6001], dashy=y_sol[4000:6001])
    
    disk_rec = nu_disk*A549_CellType/ np.sum(A549_CellType)
    ka_disk = pat.pmpc_ka(disk_rec, dsh, R_disk, 25, phi=1, L=10, ParType="Disk")

    kon_disk = popt[5]/(rhoR_disk*R_disk/(10**-9*NA)) ** nu_disk
    koff_disk = kon_disk / ka_disk
    print("nu_disk=","{:.3e}".format(nu_disk), 
          "; ka_disk=","{:.3e}".format(ka_disk), 
          "; kon_disk=","{:.3e}".format(kon_disk/60), 
          "; koff_disk=","{:.3e}".format(koff_disk/60),
          "; ke_disk=","{:.3e}".format(popt[6]/60),
          "; residence time", "{:.3e}".format(60/koff_disk),
          "Nanoparticle in the liquid", "{:.3e}".format(rho0_Psome-rhoR0)
          )
    
    """Tube"""
    length = 300
    R_tube = 40
    c_R = nu_Psome / (2*R_Psome*Rg*np.pi) 
    nu_tube = (2*np.sqrt(2*R_tube*Rg-Rg**2))*length*c_R *0.5
    print("nu_tube----->", nu_tube)
    tube_list = nu_tube*A549_CellType/ np.sum(A549_CellType)
    ka_mic = pat.pmpc_ka(tube_list, dsh, 15, 25, phi=1)
    rhoR_tube = popt[0]/nu_tube 
 
    y_init = [rho0_tube, rhoR_tube, 0, 0]

    boundary = np.array([[nu_tube,     5900,   10**-9, 10**-10], 
                        [nu_tube +0.01,   6000, 10**-2, 10**-9]])

    popt, *pcov = curve_fit(get_odeint_solution_tube, t_exp_tube_fit, theta_tube_fit, bounds=boundary)
    print("tube results----->", popt)
    t_sim = np.linspace(0, 60, 6001)
    args_tube = (popt[0], popt[1], popt[2], popt[3])

    y_sol = odeint(binding_odeint_tube, y_init, t_sim, args=args_tube)
    print("Tube Occupied Receptor", np.max(y_sol[:,2])/rhoR_tube)
    y_sol = y_sol[:,2]/np.max(y_sol[:, 2]) * 0.5 
    t_exp_tube = data_tube.values[:, 0]
    y_tube = data_tube.values[:, 1]
    error_tube = data_tube.values[:, 2]
    y_max_tube = np.max(y_tube)
    theta_tube = y_tube/y_max_tube

    pltfig(t_exp_tube, y_tube, error_tube, t_sim[0:4000], y_sol[0:4000], colors[2], 
           "Tube", dashx=t_sim[4000:6001], dashy=y_sol[4000:6001])
    print(popt)
    tube_list = popt[0]*A549_CellType/ np.sum(A549_CellType)
    ka_tube = pat.pmpc_ka(tube_list, dsh, 40, 25, phi=1, ParType="Tube", L=length)
    kon_tube = popt[1]/(rhoR_tube*15/(10**-9*NA))**popt[0]
    koff_tube = kon_tube / ka_tube

    print("nu_tube_fit=","{:.3e}".format(nu_tube), 
          "; ka_tube=","{:.3e}".format(ka_tube), 
          "; kon_tube=","{:.3e}".format(kon_tube/60), 
          "; koff_tube=","{:.3e}".format(koff_tube/60),
          "; Residence time=","{:.3e}".format(60/koff_tube)
          )


    print("nu_posome", nu_Psome)
    print(nu_disk)
    print(nu_tube)



    
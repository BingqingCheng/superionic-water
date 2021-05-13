#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy import interpolate
from scipy import integrate
from scipy.interpolate import griddata

# constants
kbK_2_ev = 0.00008661733
evperh2o_2_kjperg = 5.3602963
A2perps_2_cm2pers = 0.0001
A3GPa2ev = 0.0062415091 # A^3*GPa in eV

# Thermodynamic Intergation along isotherm
# 
# ###    $G{(P_1,T)} = G{(P_0,T)}
#     +\int_{P_0}^{P_1} \left<V\right>_{P,T} dP$
#     
# ###    $G{(P_1,T)}/T = G{(P_0,T)}/T
#     +\dfrac{1}{T}\int_{P_0}^{P_1} \left<V\right>_{P,T} dP$
#     
# ###    $G{(P_1,T)}/T = G{(P_0,T)}/T
#     +\dfrac{1}{T}\int_{\log(P_0)}^{\log(P_1)} P\left<V\right>_{P,T} d \log(P)$

def get_dmu_TdP(PTVdata, T, P0, P):
    """
    effectively TI along isotherm, as a function of P. i.e. Compute mu(P,T)/T-mu_0(P_0,T)/T
    PTVdata: an array of [P[Gpa], T[K], molar_volume[A^3/f.u.]] 
    """   
    PV = np.zeros(len(PTVdata))
    for i in range(len(PTVdata)):
        PV[i] = A3GPa2ev*PTVdata[i,2]*PTVdata[i,0]
    
    p_int = np.arange(10, 800, 5) 
    grid_z1 = griddata(PTVdata[:,[0,1]], PV, (p_int,T), method='cubic').T
    V_P = interpolate.interp1d(np.log(p_int), grid_z1)
    return (1./T)*integrate.quad(V_P, np.log(P0), np.log(P))[0]


# ### Thermodynamic Intergation along isobar
# 
# ###    $\dfrac{G{(P,T_1)}}{T_1} = \dfrac{G{(P,T_0)}}{T_0} 
#     -\int_{T_0}^{T_1} \dfrac{\left<H\right>_{P,T}}{T^2} dT$
#     
# ### $y=\ln(T/T_0)$
# 
# ### $-\int_{T_0}^{T_1} \dfrac{\left<H\right>_{P,T}}{T^2} dT=-\int_0^{\ln(T/T_0)} \left<H\right>_{P,T}/T dy$

def get_dmu_TdT(PTVEdata, P, T0, T):
    """
    effectively TI along isobar, as a function of T. i.e. Compute mu(P,T)/T-mu_0(P,T_0)/T_0
    PTHdata: an array of [P[Gpa], T[K], molar_volume[A^3/f.u.], molar_PE[eV/f.u.]] 
    """  
    T_start = max(min(T0-200,T-200),100)
    T_finish = max(T+200,T0+200)
    T_step = min((np.abs(T-T0)+50)/20.,100)
    T_int = np.arange(T_start, T_finish, T_step)
    
    HdTT = np.zeros(len(PTVEdata))
    for i in range(len(PTVEdata)):
        HdTT[i] = A3GPa2ev*PTVEdata[i,2]*PTVEdata[i,0] + PTVEdata[i,3]
        HdTT[i] /= PTVEdata[i,1]
    grid_z1 = griddata(PTVEdata[:,[0,1]], HdTT, (P,T_int), method='cubic')
    #plt.plot(np.log(T_int/T0),grid_z1)
    #print(grid_z1)
    H_T = interpolate.interp1d(np.log(T_int/T0), grid_z1)
    return -1.0*integrate.quad(H_T, 0, np.log(T/T0))[0]


# # Thermodynamic integration along both P, T, in a straight line in log(P),Log(T)
# ###    $\dfrac{G{(P_1,T)}}{T} - \dfrac{G{(P_0,T)}}{T}
#    =\dfrac{1}{T}\int_{\log(P_0)}^{\log(P_1)} P\left<V\right>_{P,T} d \log(P)$
#     
# ###    $\dfrac{G{(P,T_1)}}{T_1} - \dfrac{G{(P,T_0)}}{T_0} 
#     =-\int_0^{\ln(T/T_0)} \left<H\right>_{P,T}/T dy$


def get_dmu_TdPT(PTVEdata, P0, T0, P, T, nbin=20):
    """
    effectively TI along a straight line from [log(T0), log(P0)] to [log(T), log(P)]
    PTVEdata: an array of [P[Gpa], T[K], molar_volume[A^3/f.u.],molar_PE[eV/f.u.]] 
    """   
    
    dlogT = np.log(T)-np.log(T0)
    T_int = np.linspace(np.log(T0), np.log(T), num=nbin, endpoint=True)
    #print(np.exp(T_int))

    dlogP = np.log(P)-np.log(P0)
    P_int = np.linspace(np.log(P0), np.log(P), num=nbin, endpoint=True)
    #print(np.exp(P_int))

    x_int = np.linspace(0, 1, num=nbin, endpoint=True)
    
    PV = np.zeros(len(PTVEdata))
    for i in range(len(PTVEdata)):
        PV[i] = A3GPa2ev*PTVEdata[i,2]*PTVEdata[i,0]
        PV[i] /= PTVEdata[i,1]
        
    HdTT = np.zeros(len(PTVEdata))
    for i in range(len(PTVEdata)):
        HdTT[i] = A3GPa2ev*PTVEdata[i,2]*PTVEdata[i,0] + PTVEdata[i,3]
        HdTT[i] /= PTVEdata[i,1]
        HdTT[i] *= -1.
        
    grid_z1 = griddata(PTVEdata[:,[0,1]], HdTT*dlogT+PV*dlogP, 
                       np.vstack((np.exp(P_int),np.exp(T_int))).T, method='cubic')


    dmu_PT = interpolate.interp1d(x_int, grid_z1)
    return integrate.quad(dmu_PT, 0, 1)[0]


def get_dmu_strict(PTVEdata, P0, T0, P, T):
    """get mu(P,T) - mu(P0, T)"""
    
    try:
        # first along isobar P0
        dmu_isobar_1 = get_dmu_TdT(PTVEdata, P0, T0, T)
        # then along isotherm T
        dmu_isotherm_1 = get_dmu_TdP(PTVEdata[:,[0,1,2]], T, P0, P)

        # first along isotherm T0
        dmu_isotherm_2 = get_dmu_TdP(PTVEdata[:,[0,1,2]], T0, P0, P)
        # then along isobar P
        dmu_isobar_2 = get_dmu_TdT(PTVEdata, P, T0, T)
    except:
        return float("NaN")

    if T*abs((dmu_isobar_1+dmu_isotherm_1)- (dmu_isobar_2+dmu_isotherm_2)) < 0.002:
        return T*((dmu_isobar_1+dmu_isotherm_1) + (dmu_isobar_2+dmu_isotherm_2))/2.
    else:
        return float("NaN")



def get_dmu(PTVEdata, P0, T0, P, T, use_type=0):
    """get mu(P,T) - mu(P0, T0)"""
    
    if use_type == 0:
        try:
            dmu = T*get_dmu_TdPT(PTVEdata, P0, T0, P, T)
        except:
            dmu = float("NaN")
        return dmu
    
    elif use_type == 2:
        try:
            # first along isotherm T0
            dmu_isotherm_2 = get_dmu_TdP(PTVEdata[:,[0,1,2]], T0, P0, P)
            # then along isobar P
            dmu_isobar_2 = get_dmu_TdT(PTVEdata, P, T0, T)
            dmu_r2 = T*(dmu_isobar_2+dmu_isotherm_2)
            #print( P, T, "route2", dmu_r2)
        except:
            dmu_r2 = float("NaN")
        return dmu_r2
    
    elif use_type == 1:   
        try:
            # first along isobar P0
            dmu_isobar_1 = get_dmu_TdT(PTVEdata, P0, T0, T)
            # then along isotherm T
            dmu_isotherm_1 = get_dmu_TdP(PTVEdata[:,[0,1,2]], T, P0, P)
            dmu_r1 = T*(dmu_isobar_1+dmu_isotherm_1)
            #print( P, T, "route1", dmu_r1)
        except:
            dmu_r1 = float("NaN")
        return dmu_r1


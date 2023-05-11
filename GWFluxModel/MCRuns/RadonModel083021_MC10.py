#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 12:34:07 2022

@author: keirajohnson
"""

##for 08-30-2022##

from StreamTran import *

from pylab import *

import numpy as np

import cfc_tools as cfc

#from get_atm_conc import convert2aqueous

import get_atm_conc as atm

import stable_isotope_tools as iso

import pdb

import datetime as dt

import noble_gas_tools as ng

import matplotlib as mpl

import pickle as p

import os as os

import pandas as ps






mpl.rcParams['figure.figsize']=[10.0,8.0];

mpl.rcParams['text.usetex'] = True;

mpl.rcParams['font.size'] = 15.0;

###for error around best 10 AIC###

os.chdir(path="/Users/johnkeir/Box/Keira_Johnson/Coal_Creek/Summer")

MC_runs=ps.read_csv("MC_10_df.csv")

MC_runs_date=MC_runs[MC_runs.date=="MC_output_all_0830.csv"]

print(MC_runs_date)

Rn_list=MC_runs_date.Radon

Gas_list=MC_runs_date.Gas

print(Rn_list)











#######################################################################

#######################################################################

#inputs

#######################################################################

#######################################################################

for i in range(40,50):

    ##################################
    
    #initialize the stream transport simulation
    
    sTran = StreamTranSim('CoalCreek_072021')
    
    FieldData = {}
    
    run_tran = 1 #estimate groundwater discharge
    
    calc_age = 0 #estimate groundwater age
    
    
    #stream temperature
    
    T = 14. #max water temperature for samples
    
    T_gw = 5. #mean recharge temperature guess for now
    
    
    
    #stream elevation - didnt change this, its minimal difference
    
    E = 2979. #maximum sampling elevation in meters
    
    E_gw = 3030. #mean sampling elevation for now in meters
    
    
    
    ##################################
    
    #initialize tracer input time series for age tracers only.
    
    ##################################
    
    
    
    #make tracer input time series (cfc,sf6,He)
    
    #df_cfc_atm = cfc.get_gas_conc()
    
    #P = cfc.lapse_rate(E)
    
    #P_gw = cfc.lapse_rate(E_gw)
    
    #df_cfc = atm.convert2aqueous(df_cfc_atm,T,P=P,addHe=True,addAr39=False,addKr81=False) #cfc should pmol/kg after this un
    
    #df_cfc_gw = atm.convert2aqueous(df_cfc_atm,T_gw,P=P_gw,addHe=True,addAr39=False,addKr81=False) #cfc should pmol/kg after this un
    
    #
    
    
    
    #Ottowa_GNIP = iso.GNIP('Ottowa')
    
    #Ottowa_GNIP.read_gnip_csv('./GNIP_Ottowa.csv')
    
    #df_T = Ottowa_GNIP.data
    
    
    
    #df_tracer_gw = atm.concat_resampled_time_series(df_cfc_gw,df_T)
    
    
    
    
    
    ##sampling time - don't worry about this right now.
    
    #t_i = df_cfc.index[-1]
    
    
    
    ##################################
    
    #groundwater age data - one age for all for now
    
    #xd = np.array([40.])*1000.
    
    #yd = np.array([30.*365])  #mean age of groundwater wells sampled
    
    #C_tau = (xd,yd)
    
    #sTran.C_tau = C_tau
    
    
    
    ##################################
    
    #Q-discharge field data - gauged data
    
    ##################################
    
    
    
    xd = np.array([0.,2163.,2687.])  #distance downstream in m
    
    yd = np.array([0.02,0.03,0.08]) #discharge in cms
    
    DischData=(xd,yd)
    
    
    
    #######################################################################
    
    #tracer info
    
    #######################################################################
    
    
    
    ##################################
    
    #tracer - 1
    
    ##################################
    
    tracer = 'Rn'
    
    k_exch = 2. # Gas Exchange Coefficient m/d  #THIS IS ESTIMATED BETTER AFTER GEOMETERY IS ADDED
    
    k_exch = k_exch/60./60./24. #m/s
    
    lamma = 3.8235/60./60./24. #Zero order decay coefficients-1
    
    C_atm = 0. #Atmospheric Concentration
    
    error_perc = 0.014 #Average uncertainty as a standard error - stdev/mean in pC/L
    
    
    
    ##################################
    
    #tracer field data
    
    ##################################
    
    
    
    #river data
    
    xd = np.array([0., 640., 1376., 2163., 2574., 2687.]) #distance down stream
    
    yd = np.array([3.7,2.4,4.0,2.7,12.5,9.0]) #concentration in pC/L
    
    FieldData[tracer]=(xd,yd)
    
    
    
    #groundwater data
    
    xd = np.array([581.]) #approximate distance along the downstream direction
    
    #yd = np.array([316.]) #groundwater concentration from Meadow Spring
    
    yd = np.array([Rn_list[i]])
    
    C_gw = (xd,yd)
    
    
    
    #add tracer to simulation
    
    sTran.add_tracer(tracer)
    
    sTran.Tracers[tracer].C_atm = C_atm
    
    sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
    
    sTran.Tracers[tracer].C_gw = C_gw
    
    sTran.Tracers[tracer].k = k_exch
    
    sTran.Tracers[tracer].error_perc = error_perc
    
    
    ##################################
    
    #tracer - 2
    
    ##################################
    #=============================================================================
    
    tracer = 'SO4'
    
    k_exch = 0. # gas exchange coefficient m/s - zero for a non-gaseous tracer
    
    lamma = 0. #s-1
    
    C_atm = 0.
    
    error_perc = 0.05 #unknown for Cl
    
    
    
    
    
    ##################################
    
    #tracer field data
    
    ##################################
    
    
    
    #river data
    
    xd = np.array([0., 640., 1376., 2163., 2574., 2687.]) #distance down stream
    
    yd = np.array([88.21, 85.74, 84.27, 83.9, 72.44, 162.47]) #UM
    
    FieldData[tracer]=(xd,yd)
    
    
    
    #groundwater data
    
    xd = np.array([581.])
    
    yd = np.array([366.2]) #366.2
    
    C_gw = (xd,yd)
    
    
    
    #add tracer to simulation
    
    sTran.add_tracer(tracer)
    
    sTran.Tracers[tracer].C_atm = C_atm
    
    sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
    
    sTran.Tracers[tracer].C_gw = C_gw
    
    sTran.Tracers[tracer].k = k_exch
    
    sTran.Tracers[tracer].error_perc = error_perc
    # =============================================================================
    
    
    
    
    
    ##################################
    
    #tracer - 3
    
    #add tracers at will by copying and pasting and changing the inputs appropriately
    
    tracer = 'd18O'
    
    k_exch = 0. # gas exchange coefficient m/s - zero for a non-gaseous tracer
    
    lamma = 0. #s-1
    
    C_atm = 0.
    
    error_perc = 0.1 #unknown for Cl
    
    
    
    
    
    ##################################
    
    #tracer field data
    
    ##################################
    
    
    
    #river data
    
    xd = np.array([0., 640., 1376., 2163., 2574., 2687.]) #distance down stream
    
    yd = np.array([-14.6, -14.55, -14.59, -14.47, -14.96, -15.38]) #UM
    
    FieldData[tracer]=(xd,yd)
    
    
    
    #groundwater data
    
    xd = np.array([581.])
    
    yd = np.array([-16.9]) #-16.9
    
    C_gw = (xd,yd)
    
    
    
    #add tracer to simulation
    
    sTran.add_tracer(tracer)
    
    sTran.Tracers[tracer].C_atm = C_atm
    
    sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
    
    sTran.Tracers[tracer].C_gw = C_gw
    
    sTran.Tracers[tracer].k = k_exch
    
    sTran.Tracers[tracer].error_perc = error_perc
    
    
    
    ##################################
    
    #tracer - 4
    
    ##################################
    #=============================================================================
    
    tracer = 'Ca'
    
    k_exch = 0. # gas exchange coefficient m/s - zero for a non-gaseous tracer
    
    lamma = 0. #s-1
    
    C_atm = 0.
    
    error_perc = 0.05 #unknown for Cl
    
    
    
    
    
    ##################################
    
    #tracer field data
    
    ##################################
    
    
    
    #river data
    
    xd = np.array([0., 640., 1376., 2163., 2574., 2687.]) #distance down stream
    
    yd = np.array([11.91, 11.92, 12.06, 12.01, 11.32, 15.55]) #UM
    
    FieldData[tracer]=(xd,yd)
    
    
    
    #groundwater data
    
    xd = np.array([581.])
    
    yd = np.array([23.2]) #23.2
    
    C_gw = (xd,yd)
    
    
    
    #add tracer to simulation
    
    sTran.add_tracer(tracer)
    
    sTran.Tracers[tracer].C_atm = C_atm
    
    sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
    
    sTran.Tracers[tracer].C_gw = C_gw
    
    sTran.Tracers[tracer].k = k_exch
    
    sTran.Tracers[tracer].error_perc = error_perc
    # =============================================================================
    
    ##################################
    
    #tracer - 5
    
    ##################################
    #=============================================================================
    
    tracer = 'Na'
    
    k_exch = 0. # gas exchange coefficient m/s - zero for a non-gaseous tracer
    
    lamma = 0. #s-1
    
    C_atm = 0.
    
    error_perc = 0.05 #unknown for Cl
    
    
    
    
    
    ##################################
    
    #tracer field data
    
    ##################################
    
    
    
    #river data
    
    xd = np.array([0., 640., 1376., 2163., 2574., 2687.]) #distance down stream
    
    yd = np.array([1.47, 1.54, 1.66, 1.78, 1.96, 2.01]) #UM
    
    FieldData[tracer]=(xd,yd)
    
    
    
    #groundwater data
    
    xd = np.array([581.])
    
    yd = np.array([2.2])
    
    C_gw = (xd,yd)
    
    
    
    #add tracer to simulation
    
    sTran.add_tracer(tracer)
    
    sTran.Tracers[tracer].C_atm = C_atm
    
    sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
    
    sTran.Tracers[tracer].C_gw = C_gw
    
    sTran.Tracers[tracer].k = k_exch
    
    sTran.Tracers[tracer].error_perc = error_perc
    # =============================================================================
    
    ##################################
    
    #tracer - 6
    
    ##################################
    #=============================================================================
    
    tracer = 'Si'
    
    k_exch = 0. # gas exchange coefficient m/s - zero for a non-gaseous tracer
    
    lamma = 0. #s-1
    
    C_atm = 0.
    
    error_perc = 0.05 #unknown for Cl
    
    
    
    
    
    ##################################
    
    #tracer field data
    
    ##################################
    
    
    
    #river data
    
    xd = np.array([0., 640., 1376., 2163., 2574., 2687.]) #distance down stream
    
    yd = np.array([1.23, 1.31, 1.36, 1.40, 1.97, 2.4]) #UM
    
    FieldData[tracer]=(xd,yd)
    
    
    
    #groundwater data
    
    xd = np.array([581.])
    
    yd = np.array([4.]) #4
    
    C_gw = (xd,yd)
    
    
    
    #add tracer to simulation
    
    sTran.add_tracer(tracer)
    
    sTran.Tracers[tracer].C_atm = C_atm
    
    sTran.Tracers[tracer].C_us = FieldData[tracer][1][0]
    
    sTran.Tracers[tracer].C_gw = C_gw
    
    sTran.Tracers[tracer].k = k_exch
    
    sTran.Tracers[tracer].error_perc = error_perc
    # =============================================================================
    
    
    #######################################################################
    
    #model geometry
    
    #######################################################################
    
    
    
    L = 2687.  #Modeled Total Lengthm
    
    nx = int(1.e4)   #number of cells - can keep this constant for now
    
    Evap = 4.9e-3/60./60./24. #Evaporation Rate
    
    Q_us = 0.02 #measured
    
    
    
    ##################################
    
    # geometery visually estimated at sampling locations
    
    ##################################
    
    #widths
    
    xd = np.array([0.,2163.,2687.]) #distance downstream
    
    yd = np.array([1.87, 1.51, 2.68]) #average width
    
    w = (xd,yd)
    
    
    
    #depths
    
    xd = np.array([0.,2163.,2687.]) #distance downstream
    
    yd = np.array([0.19, 0.11, 0.22]) #average depth
    
    d = (xd,yd)  #([dist.]),([depths])
    
    A = (w[0],w[1]*d[1])
    
    
    
    
    
    
    
    ##################################
    
    #Estimate Gas Exchange Coefficient from Geometry
    
    ##################################
    
    #Calculate the gas exchange coefficients from Raymond 2012 - equation 7.
    
    E_us = 2979. #elevation upstream
    
    E_ds = 2905. #elevation downstream
    
    
    
    #no need to change below
    
    S = (E_us-E_ds)/d[0][-1]
    
    Q_av = Q_us
    
    A_av = A[1].mean()
    
    V_av = Q_av/A_av
    
    D_av = d[1].mean()
    
    k600 = ng.k600(Q_av,V_av,D_av,S,eqn=7)
    
    Sch_Rn = ng.schmidt('Rn',T)
    
    k_Rn = (Sch_Rn/600.)**-0.5*k600
    
    #k_Rn = 97
    
    k_Rn=Gas_list[i]
    
    #k_Rn = 3.85 #lower MCMC 95%
    
    #k_Rn = 5.42 #lower MCMC 95%
    
    sTran.Tracers['Rn'].k = k_Rn/24./60./60.
    
    print('Rn Gas Exchange = ', sTran.Tracers['Rn'].k*60*60*24, ' m/d')
    
    
    
    
    
    ##################################
    
    #Tributary info
    
    #tributary (np.array([x]),[disch],{tracer:[C]})
    
    Q_trib = (np.array([2599.]),[0.01],{'Rn':[2.1],'SO4':[364.08],'d18O':[-16.54], 'Ca':[24.603], 'Na':[2.109], 'Si':[3.513]})
    
    
    
    
    
    ##################################
    
    #Groundwater inflows
    
    #inflows
    
    #just does even groundwater inflow steps, with the number of steps equal to the number of sampling points
    
    ql = 1.e-7*np.ones(len(FieldData['Rn'][0]))
    
    
    
    
    
    ##################################
    
    #end model input!!!
    
    ##################################
    
    
    
    #RUNTIME
    
    #######################################################################
    
    #Parameterize Simulation
    
    #######################################################################
    
    sTran.L = L
    
    sTran.nx = nx
    
    sTran.wi = w
    
    sTran.Ai = A
    
    sTran.q_lin = ql
    
    sTran.Q_us = Q_us
    
    sTran.Q_trib = Q_trib
    
    #sTran.C_t = df_cfc
    
    #sTran.C_t = df_tracer_gw
    
    sTran.tau_flag = 1
    
    #sTran.t_i = t_i
    
    sTran.age_dist = 'dispersion'
    
    
    
    
    
    
    
    #######################################################################
    
    
    
    #Run Simulation
    
    
    
    #######################################################################
    
    
    
    #sTran.CalcTran()
    
    
    
    print('Calculating groundwater inflow')
    
    
    
    tracer_list=['Rn']
    
    if run_tran:
    
        sTran.CalcGwInflow(tracer_list,FieldData,DischData,DischargeError=0.01)
    
        f = open('CoalCreekResults_v1.pkl','wb')
    
        p.dump(sTran,f)
    
        f.close()
    
    
    
    else:
    
        f = open('CoalCreekResults_v1.pkl','rb')
    
        sTran = p.load(f)
    
        f.close()
    
    
    final_data = None
    with open('CoalCreekResults_v1.pkl', 'rb') as f:
        data = p.load(f)
        sTran = data
        #Figure 1
        ##################################
        
        radon = sTran.SimRes.tracers['Rn']
        #chloride = sTran.SimRes.tracers['Cl']
        
    
        def extract_stuff(obj):
            keys = []
            extracted = []
            for k,v in obj.__dict__.items():
                if k=='name':
                    continue
                keys.append(k)
                extracted.append(v.value)
                
            return (keys, extracted)
    
        
        discharge_keys, discharge_data = (["Discharge"], [sTran.SimRes.Q.value])
        radon_keys, radon_data = extract_stuff(radon)
        #chloride_keys, chloride_data = extract_stuff(chloride)
        gwDischarge_keys, gwDischarge_data = (["GWDischarge"], [sTran.SimRes.q_lin.value])
        
        def export_to_csv(keys, data, name):
            with open(f"{name}_exported_083021_best10AIC.csv", "w+") as datafile:
                datafile.write(",".join(keys) + '\n')
                numpy_stuff = np.asarray(data).T
                for row in numpy_stuff:
                    datafile.write(','.join([str(x) for x in row]) + '\n')
                    
                csv_string = ','.join(radon_keys) + '\n'
                
                #for row in numpy_stuff:
                    #datafile.write(','.join([str(x) for x in row]) + '\n')
                    
                    
        export_to_csv(radon_keys, radon_data, "Radon_"+str(Rn_list[i])+"_"+str(Gas_list[i]))
        #export_to_csv(chloride_keys, chloride_data, "Chloride")
        export_to_csv(discharge_keys, discharge_data, "Discharge_"+str(Rn_list[i])+"_"+str(Gas_list[i]))
        export_to_csv(gwDischarge_keys, gwDischarge_data, "GWDischarge_"+str(Rn_list[i])+"_"+str(Gas_list[i]))
    
    



# #######################################################################

# #visulization the solution

# #######################################################################

    


# #######################################################################

# #Figure 1

# fig1 = figure(figsize=(4,4))

# ##################################

# #discharge

# error_perc=0.05

# x = np.linspace(0,L,num=nx)

# subplot(3,1,1)

# plot(x/1000.,sTran.SimRes.Q.value,'k',lw=2,alpha=0.6,label='modeled')





# #field data

# errorbar(DischData[0]/1000.,DischData[1],yerr=error_perc*DischData[1],fmt='ks',mec='k',label='Data')

# grid(b=True,which='major',linestyle='--')

# ylabel('discharge (m3/s)')



# ##################################

# #Rn

# tracer = 'Rn'

# error_perc = 0.15

# subplot(3,1,2)

# plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'k',alpha=0.6,lw=2,label='modeled')

# ylabel("Rn (piC/L)")





# #plot field data

# errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='ks',mec='k',label='Data')

# grid(b=True,which='major',linestyle='--')



# ##################################

# #=============================================================================
# #SO4

# #tracer = 'SO4'

# #error_perc = 0.15

# #subplot(8,1,3)

# #plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'k',alpha=0.6,lw=2,label='modeled')

# #ylabel(tracer)







# #plot field data

# #errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='ks',mec='k',label='Data')

# #grid(b=True,which='major',linestyle='--')
# #=============================================================================


# ##################################

# #d18O

# #tracer = 'd18O'

# #error_perc = 0.05

# #subplot(8,1,4)

# #plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'k',alpha=0.6,lw=2,label='modeled')

# #ylabel(tracer)




# #plot field data

# #errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='ks',mec='k',label='Data')

# #grid(b=True,which='major',linestyle='--')

# ##################################

# #=============================================================================
# #Ca

# #tracer = 'Ca'

# #error_perc = 0.15

# #subplot(8,1,5)

# #plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'k',alpha=0.6,lw=2,label='modeled')

# #ylabel(tracer)







# #plot field data

# #errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='ks',mec='k',label='Data')

# #grid(b=True,which='major',linestyle='--')

# ##################################

# #=============================================================================
# #Na

# #tracer = 'Na'

# #error_perc = 0.15

# #subplot(8,1,6)

# #plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'k',alpha=0.6,lw=2,label='modeled')

# #ylabel(tracer)







# #plot field data

# #errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='ks',mec='k',label='Data')

# #grid(b=True,which='major',linestyle='--')

# ##################################

# #=============================================================================
# #Si

# #tracer = 'Si'

# #error_perc = 0.15

# #subplot(8,1,7)

# #plot(x/1000.,sTran.SimRes.tracers[tracer].C.value,'k',alpha=0.6,lw=2,label='modeled')

# #ylabel(tracer)







# #plot field data

# #errorbar(FieldData[tracer][0]/1000.,FieldData[tracer][1],yerr=error_perc*FieldData[tracer][1], fmt='ks',mec='k',label='Data')

# #grid(b=True,which='major',linestyle='--')

# ##################################

# #Groundwater inflow

# subplot(3,1,3)

# plot(x/1000.,sTran.SimRes.q_lin.value,'k',lw=2,alpha=0.6,label='Groundwater Disch.')

# ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# grid(b=True,which='major',linestyle='--')

# xlabel('distance (km)')

# ylabel('lateral discharge (m/s)')

# #tight_layout()

# fig1.set_size_inches(8,7)

# fig1.savefig('PhysicalParamsFigure_083021_medMC.png',dpi=200,bbox_inches='tight')











########################################################################






#show()

#sTran.SimRes.dump2xl()




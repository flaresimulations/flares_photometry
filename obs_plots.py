import numpy as np
import pandas as pd
from uncertainties import unumpy

observation_path = './Obs_data1/'

RR_2014 = pd.read_table(observation_path+'Remy_Ruyer_2014_KINGFISH_z0.txt',sep='\t',comment='#',names=['DM', 'DM_err', 'SM', 'HI', 'HI_err', 'Metals', 'H2_mw', 'H2_z'])
RR_2014['DM_err_actual'] = RR_2014['DM']*(RR_2014['DM_err']/100.0)
RR_2014['DTG'] = RR_2014['DM'] / (RR_2014['HI'] + RR_2014['H2_mw'])
RR_2014['DTM'] = RR_2014['DM'] / (((10**(RR_2014['Metals'] - 8.69))*0.0134)*(RR_2014['HI'] + RR_2014['H2_mw']))

#RR2015 z=0
RR_2015 = pd.read_table(observation_path+'RR2015.txt',sep='\t',comment='#',index_col=False,names=['SM','SM_err','DM_1','DM_1_up','DM_1_down','DM_2', 'DM_2_up', 'DM_2_down', 'HI', 'HI_err', 'Oxygen', 'H2_mw','H2_z'])

RR_2015['DTG_1A'] = RR_2015['DM_1'] - np.log10(RR_2015['HI'] + RR_2015['H2_z'])
RR_2015['DTG_1B'] = RR_2015['DM_1'] - np.log10(RR_2015['HI'] + RR_2015['H2_mw'])
RR_2015['DTG_2A'] = RR_2015['DM_2'] - np.log10(RR_2015['HI'] + RR_2015['H2_z'])
RR_2015['DTG_2B'] = RR_2015['DM_2'] - np.log10(RR_2015['HI'] + RR_2015['H2_mw'])

RR_2015['DTM_1A'] = RR_2015['DM_1'] - np.log10((10**(RR_2015['Oxygen'] - 8.69 - 1.87289520164))*(RR_2015['HI'] + RR_2015['H2_z']))
RR_2015['DTM_1B'] = RR_2015['DM_1'] - np.log10((10**(RR_2015['Oxygen'] - 8.69 - 1.87289520164))*(RR_2015['HI'] + RR_2015['H2_mw']))
RR_2015['DTM_2A'] = RR_2015['DM_2'] - np.log10((10**(RR_2015['Oxygen'] - 8.69 - 1.87289520164))*(RR_2015['HI'] + RR_2015['H2_z']))
RR_2015['DTM_2B'] = RR_2015['DM_2'] - np.log10((10**(RR_2015['Oxygen'] - 8.69 - 1.87289520164))*(RR_2015['HI'] + RR_2015['H2_mw']))

#Santini 2014 z=0,1,2
Santini_2014_z0 = pd.read_table(observation_path+'Santini_2014_z0.txt',sep='\t',comment='#',index_col=False,names=['SM_z0', 'DM_z0', 'DM_up_err_z0', 'DM_down_err_z0'])
Santini_2014_z1 = pd.read_table(observation_path+'Santini_2014_z1.txt',sep='\t',comment='#',index_col=False,names=['SM_z1', 'DM_z1', 'DM_up_err_z1', 'DM_down_err_z1'])
Santini_2014_z2 = pd.read_table(observation_path+'Santini_2014_z2.txt',sep='\t',comment='#',index_col=False,names=['SM_z2', 'DM_z2', 'DM_up_err_z2', 'DM_down_err_z2'])

Santini_2014 = pd.concat((Santini_2014_z0,Santini_2014_z1,Santini_2014_z2), axis=1)

#daCunha 2015
daCunha_z2 = pd.read_table(observation_path+'daCunha_2015_z_2.txt',sep='\t',comment='#',index_col=False,names=['z_z2','z_up_err_z2','z_down_err_z2','SM_z2','SM_up_err_z2','SM_down_err_z2','DM_z2','DM_up_err_z2','DM_down_err_z2'])
daCunha_z3 = pd.read_table(observation_path+'daCunha_2015_z_3.txt',sep='\t',comment='#',index_col=False,names=['z_z3','z_up_err_z3','z_down_err_z3','SM_z3','SM_up_err_z3','SM_down_err_z3','DM_z3','DM_up_err_z3','DM_down_err_z3'])
daCunha_z4 = pd.read_table(observation_path+'daCunha_2015_z_4.txt',sep='\t',comment='#',index_col=False,names=['z_z4','z_up_err_z4','z_down_err_z4','SM_z4','SM_up_err_z4','SM_down_err_z4','DM_z4','DM_up_err_z4','DM_down_err_z4'])
daCunha_z5 = pd.read_table(observation_path+'daCunha_2015_z_5.txt',sep='\t',comment='#',index_col=False,names=['z_z5','z_up_err_z5','z_down_err_z5','SM_z5','SM_up_err_z5','SM_down_err_z5','DM_z5','DM_up_err_z5','DM_down_err_z5'])
daCunha_z6 = pd.read_table(observation_path+'daCunha_2015_z_6.txt',sep='\t',comment='#',index_col=False,names=['z_z6','z_up_err_z6','z_down_err_z6','SM_z6','SM_up_err_z6','SM_down_err_z6','DM_z6','DM_up_err_z6','DM_down_err_z6'])

daCunha_2015 = pd.concat((daCunha_z2,daCunha_z3,daCunha_z4,daCunha_z5,daCunha_z6), axis=1)

#Mancini2015
Mancini_2015 = pd.read_table(observation_path+'Mancini_2015_z6_z7.txt',sep=',',comment='#',index_col=False,names=['z', 'SM', 'SM_err', 'DM', 'DM_err', 'uplims'])

#Bourne2012
Bourne_2012 = pd.read_table(observation_path+'Bourne2012_z0.txt',sep='\t',comment='#',index_col=False,names=['MEDZ','MEDM','MEDC','MLIMITS_DOWN','MLIMITS_UP','MEDMSTAR','MEDMSTARERR','MEDMDUST','MEDMDUSTERR','MEDDMS','MEDDMSERR','NBIN'])

#Ciesla2014
Ciesla_2014 = pd.read_table(observation_path+'Ciesla_2014_z0.txt',sep='\t',comment='#',index_col=False,names=['ID1','DM','DM_err','ID2','SM'])

#-------------Wiseman2016
Wiseman_z2 = pd.read_table(observation_path+'wiseman2016_z2.txt',sep='\t',comment='#',index_col=False,names=['z_z2','SM_z2','SM_uperr_z2','SM_downerr_z2','SFR_z2','SFR_err_z2','Metals_z2','Metals_err_z2','DTM_z2','DTM_err_z2'])
Wiseman_z3 = pd.read_table(observation_path+'wiseman2016_z3.txt',sep='\t',comment='#',index_col=False,names=['z_z3','SM_z3','SM_uperr_z3','SM_downerr_z3','SFR_z3','SFR_err_z3','Metals_z3','Metals_err_z3','DTM_z3','DTM_err_z3'])
Wiseman_z4 = pd.read_table(observation_path+'wiseman2016_z4.txt',sep='\t',comment='#',index_col=False,names=['z_z4','SM_z4','SM_uperr_z4','SM_downerr_z4','SFR_z4','SFR_err_z4','Metals_z4','Metals_err_z4','DTM_z4','DTM_err_z4'])

Wiseman_2017 = pd.concat((Wiseman_z2,Wiseman_z3,Wiseman_z4),axis=1)

#-------------Clemens2013
Clemens_2013 = pd.read_table(observation_path+'Clemens2013.txt',sep='\t',comment='#',index_col=False,names=['DM', 'Phi', 'Phi_up_err', 'Phi_down_err'])

#-------------VlahakisA_2005
VlahakisA_2005 = pd.read_table(observation_path+'Vlahakis2005A.txt',sep='\t',comment='#',index_col=False,names=['DM', 'Phi', 'Phi_up_err', 'Phi_down_err'])

#-------------VlahakisA_2005
VlahakisB_2005 = pd.read_table(observation_path+'Vlahakis2005B.txt',sep='\t',comment='#',index_col=False,names=['DM', 'Phi', 'Phi_up_err', 'Phi_down_err'])

def DM_obs(axs, i):
    
    if i == 0:
        import simondata as hst
        import dustpedia_data as dp_dat
        hubbletype = dp_dat.hubbletype
        ETG = np.where(hubbletype < 0)
        LTG = np.where(hubbletype >= 0)
        axs.errorbar(unumpy.nominal_values(dp_dat.Mstar_all[ETG]), unumpy.nominal_values(dp_dat.Mdust_all[ETG]), xerr = unumpy.std_devs(dp_dat.Mstar_all[ETG]), yerr = unumpy.std_devs(dp_dat.Mdust_all[ETG]), color='pink',label=r'$\mathrm{DustPedia}$ $\mathrm{(ETG)}$',fmt='.', alpha = 0.6, elinewidth=0.05)
        axs.errorbar(unumpy.nominal_values(dp_dat.Mstar_all[LTG]), unumpy.nominal_values(dp_dat.Mdust_all[LTG]), xerr = unumpy.std_devs(dp_dat.Mstar_all[LTG]), yerr = unumpy.std_devs(dp_dat.Mdust_all[LTG]), color='violet',label=r'$\mathrm{DustPedia}$ $\mathrm{(LTG)}$',fmt='.', alpha = 0.8,elinewidth=0.05)
        axs.errorbar(RR_2015['SM'], RR_2015['DM_1'], xerr = RR_2015['SM_err'], yerr = (RR_2015['DM_1_down'], RR_2015['DM_1_up']),color='green',label=r'$\mathrm{Remy-Ruyer}$ $\mathrm{2015}$',fmt='.', alpha = 1)
        #axs[i].errorbar(np.log10(Bourne_2012['MEDMSTAR']), np.log10(Bourne_2012['MEDMDUST']), yerr = np.log10(Bourne_2012['MEDMDUSTERR']/Bourne_2012['MEDMDUST']) , color='orange',label=r'$\mathrm{Bourne}$ $\mathrm{2012}$',fmt='.')
        axs.errorbar(Ciesla_2014['SM'], Ciesla_2014['DM'], yerr = Ciesla_2014['DM_err'] , color='red',label=r'$\mathrm{Ciesla}$ $\mathrm{2014}$',fmt='.', alpha = 1)
        axs.errorbar(Santini_2014['SM_z0'], Santini_2014['DM_z0'], yerr = (Santini_2014['DM_down_err_z0'], Santini_2014['DM_up_err_z0']),label=r'$\mathrm{Santini}$ $\mathrm{2014}$',color='blue',fmt='.', alpha = 1)
    if i == 1:
        axs.errorbar(Santini_2014['SM_z1'], Santini_2014['DM_z1'], yerr = (Santini_2014['DM_down_err_z1'], Santini_2014['DM_up_err_z1']), color='blue',label=r'$\mathrm{Santini}$ $\mathrm{2014}$',fmt='.', alpha = 1)
        
        #ok = np.where(hst.z == 1)
        #axs.errorbar(hst.Mstar[ok], hst.Mdust[ok], xerr = [hst.Mstar_l[ok], hst.Mstar_u[ok]], yerr = [hst.Mdust_l[ok], hst.Mdust_u[ok]], color='indigo',label=r'$\mathrm{S.}$ $\mathrm{P.}$ $\mathrm{Driver}$ $\mathrm{2017}$',fmt='.', alpha = 0.7, elinewidth=0.05)
    if i == 2:
        axs.errorbar(Santini_2014['SM_z2'], Santini_2014['DM_z2'], yerr = (Santini_2014['DM_down_err_z2'], Santini_2014['DM_up_err_z2']), color='blue',label=r'$\mathrm{Santini}$ $\mathrm{2014}$',fmt='.', alpha = 1)
        axs.errorbar(daCunha_2015['SM_z2'], daCunha_2015['DM_z2'], xerr = (daCunha_2015['SM_down_err_z2'], daCunha_2015['SM_up_err_z2']), yerr = (daCunha_2015['DM_down_err_z2'], daCunha_2015['DM_up_err_z2']), color='crimson',label=r'$\mathrm{daCunha}$ $\mathrm{2015}$',fmt='.', alpha = 1)
        #ok = np.where(hst.z == 2)
        #axs.errorbar(hst.Mstar[ok], hst.Mdust[ok], xerr = [hst.Mstar_l[ok], hst.Mstar_u[ok]], yerr = [hst.Mdust_l[ok], hst.Mdust_u[ok]], color='indigo',label=r'$\mathrm{S.}$ $\mathrm{P.}$ $\mathrm{Driver}$ $\mathrm{2017}$',fmt='.', alpha = 0.7, elinewidth=0.05)
    if i == 3:
        axs.errorbar(daCunha_2015['SM_z3'], daCunha_2015['DM_z3'], xerr = (daCunha_2015['SM_down_err_z3'], daCunha_2015['SM_up_err_z3']), yerr = (daCunha_2015['DM_down_err_z3'], daCunha_2015['DM_up_err_z3']), color='crimson',label=r'$\mathrm{daCunha}$ $\mathrm{2015}$',fmt='.', alpha = 1)
        
        #ok = np.where(hst.z == 3)
        #axs.errorbar(hst.Mstar[ok], hst.Mdust[ok], xerr = [hst.Mstar_l[ok], hst.Mstar_u[ok]], yerr = [hst.Mdust_l[ok], hst.Mdust_u[ok]], color='indigo',label=r'$\mathrm{S.}$ $\mathrm{P.}$ $\mathrm{Driver}$ $\mathrm{2017}$',fmt='.', alpha = 0.7, elinewidth=0.05)
    if i == 4:
        axs.errorbar(daCunha_2015['SM_z4'], daCunha_2015['DM_z4'], xerr = (daCunha_2015['SM_down_err_z4'], daCunha_2015['SM_up_err_z4']), yerr = (daCunha_2015['DM_down_err_z4'], daCunha_2015['DM_up_err_z4']), color='crimson',label=r'$\mathrm{daCunha}$ $\mathrm{2015}$',fmt='.', alpha = 1)
        
        #ok = np.where(hst.z == 4)
        #axs.errorbar(hst.Mstar[ok], hst.Mdust[ok], xerr = [hst.Mstar_l[ok], hst.Mstar_u[ok]], yerr = [hst.Mdust_l[ok], hst.Mdust_u[ok]], color='indigo',label=r'$\mathrm{S.}$ $\mathrm{P.}$ $\mathrm{Driver}$ $\mathrm{2017}$',fmt='.', alpha = 0.7, elinewidth=0.05)
    if i == 5:
        axs.errorbar(daCunha_2015['SM_z5'], daCunha_2015['DM_z5'], xerr = (daCunha_2015['SM_down_err_z5'], daCunha_2015['SM_up_err_z5']), yerr = (daCunha_2015['DM_down_err_z5'], daCunha_2015['DM_up_err_z5']), color='crimson',label=r'$\mathrm{daCunha}$ $\mathrm{2015}$',fmt='.', alpha = 1)
        
        #ok = np.where(hst.z == 5)
        #axs.errorbar(hst.Mstar[ok], hst.Mdust[ok], xerr = [hst.Mstar_l[ok], hst.Mstar_u[ok]], yerr = [hst.Mdust_l[ok], hst.Mdust_u[ok]], color='indigo',label=r'$\mathrm{S.}$ $\mathrm{P.}$ $\mathrm{Driver}$ $\mathrm{2017}$',fmt='.', alpha = 0.7, elinewidth=0.05)
    if i == 6:
        axs.errorbar(daCunha_2015['SM_z6'], daCunha_2015['DM_z6'], xerr = (daCunha_2015['SM_down_err_z6'], daCunha_2015['SM_up_err_z6']), yerr = (daCunha_2015['DM_down_err_z6'], daCunha_2015['DM_up_err_z6']), color='crimson',label=r'$\mathrm{daCunha}$ $\mathrm{2015}$',fmt='.', alpha = 1)
        
        #ok = np.where(hst.z == 6)
        #axs.errorbar(hst.Mstar[ok], hst.Mdust[ok], xerr = [hst.Mstar_l[ok], hst.Mstar_u[ok]], yerr = [hst.Mdust_l[ok], hst.Mdust_u[ok]], color='indigo',label=r'$\mathrm{S.}$ $\mathrm{P.}$ $\mathrm{Driver}$ $\mathrm{2017}$',fmt='.', alpha = 0.7, elinewidth=0.05)
    if i == 7:
        axs.errorbar(Mancini_2015['SM'], Mancini_2015['DM'], yerr = Mancini_2015['DM_err'], xerr = Mancini_2015['SM_err'], uplims = Mancini_2015['uplims'], color='brown',label=r'$\mathrm{Mancini}$ $\mathrm{2015}$',fmt='.', alpha = 1, marker = ".", markersize = 6)
        
        
def D_Met_obs(axs, i):

    if(i == 0):
        axs.errorbar(RR_2015['Oxygen'], RR_2015['DM_1'], yerr = (RR_2015['DM_1_down'], RR_2015['DM_1_up']),color='g',label=r'$\mathrm{Remy-Ruyer}$ $\mathrm{2015}$',fmt='.')
        
def DG_Mstar_obs(axs, i):

    if(i == 0):
        RR_DTG1B_err = 10**RR_2015['DTG_1B'] * np.sqrt( (RR_2015['DM_1_up']/RR_2015['DM_1'])**2 + (RR_2015['HI_err']/RR_2015['HI'])**2 )
        log_RR_DTG1B_err =0.434* (np.log10(RR_DTG1B_err)/RR_2015['DTG_1B'])
        axs.errorbar(RR_2015['SM'], RR_2015['DTG_1B'], yerr=(log_RR_DTG1B_err), color='g',label=r'$\mathrm{Remy-Ruyer}$ $\mathrm{2015}$',fmt='.') 
 
def DG_met_obs(axs, i):        
    
    if(i == 0):
        RR_DTG1B_err = 10**RR_2015['DTG_1B'] * np.sqrt( (RR_2015['DM_1_up']/RR_2015['DM_1'])**2 + (RR_2015['HI_err']/RR_2015['HI'])**2 )
        log_RR_DTG1B_err =0.434* (np.log10(RR_DTG1B_err)/RR_2015['DTG_1B'])
        axs.errorbar(RR_2015['Oxygen'], RR_2015['DTG_1B'], yerr=(log_RR_DTG1B_err), color='g',label=r'$\mathrm{Remy-Ruyer}$ $\mathrm{2015}$',fmt='.') 

def DMF(axs, i):

    if(i == 0):
        axs.errorbar(VlahakisA_2005['DM'], np.log10(VlahakisA_2005['Phi']), yerr=[np.log10((VlahakisA_2005['Phi']+VlahakisA_2005['Phi_down_err'])/VlahakisA_2005['Phi']),np.log10((VlahakisA_2005['Phi']+VlahakisA_2005['Phi_up_err'])/VlahakisA_2005['Phi'])] ,color='g',label=r'$\mathrm{VlahakisA}$ $\mathrm{2005}$',fmt='.')
        axs.errorbar(VlahakisB_2005['DM'] ,np.log10(VlahakisB_2005['Phi']) ,yerr=[np.log10((VlahakisB_2005['Phi']+VlahakisB_2005['Phi_down_err'])/VlahakisB_2005['Phi']),np.log10((VlahakisB_2005['Phi']+VlahakisB_2005['Phi_up_err'])/VlahakisB_2005['Phi'])] ,color='red',label=r'$\mathrm{VlahakisB}$ $\mathrm{2005}$',fmt='.')
        axs.errorbar(Clemens_2013['DM'], Clemens_2013['Phi'], xerr = np.ones(len(Clemens_2013['DM']))*0.3, yerr = [Clemens_2013['Phi_down_err'],Clemens_2013['Phi_up_err']],color='b',label=r'$\mathrm{Clemens}$ $\mathrm{2013}$',fmt='.')
    if (i == 1):
        data = np.genfromtxt('./Obs_data/Eales2009.txt', skip_header=2, delimiter = ',')
        x, y, y_up, y_low = np.log10(data[:,0]), np.log10(data[:,1]), np.log10(data[:,2])-np.log10(data[:,1]), np.log10(data[:,1])-np.log10(data[:,3])
        axs.errorbar(x, y, yerr=[y_low, y_up] ,color='g',label=r'$\mathrm{Eales}$ $\mathrm{2009}$',fmt='.')
    
def DTM_oxy(axs, i):
    
    import Wiseman_data_OH
    Wiseman_OH = Wiseman_data_OH.OH_abund
    Wiseman_z = Wiseman_data_OH.z
    Wiseman_DTM = Wiseman_data_OH.DTM_all
    
    import De_Cia_data 
    Cia_OH = De_Cia_data.OH_abund
    Cia_z = De_Cia_data.z
    Cia_DTM = De_Cia_data.DTM_all
    Cia_xuplims = De_Cia_data.xuplims
    Cia_yuplims = De_Cia_data.yuplims
    
    if(i == 0):
        x = RR_2015['Oxygen']
        y = RR_2015['DTM_1B']
        yerror = np.abs((0.75*RR_2015['DTM_1B']))
        yall = unumpy.uarray(y, yerror)
        ynew = yall - unumpy.log10(1. + 10**yall)
        axs.errorbar(x, unumpy.nominal_values(ynew), yerr=unumpy.std_devs(ynew), color='g',label=r'$\mathrm{Remy-Ruyer}$ $\mathrm{2015}$',fmt='o') 
    if(i == 2):
        ok = np.where(Wiseman_z == 2)
        x = Wiseman_OH[ok]
        y = Wiseman_DTM[ok]
        axs.errorbar(unumpy.nominal_values(x), unumpy.nominal_values(y), xerr =  unumpy.std_devs(x), yerr = unumpy.std_devs(y), color='r',label=r'$\mathrm{Wiseman}$ $\mathrm{2017}$',fmt='o')
        
        ok = np.where(Cia_z == 2)
        x = Cia_OH[ok]
        y = Cia_DTM[ok]
        xuplims = Cia_xuplims[ok]
        yuplims = Cia_yuplims[ok]
        up = np.logical_and(xuplims == 0, yuplims == 0)
        axs.errorbar(unumpy.nominal_values(x[up]), unumpy.nominal_values(y[up]), xerr =  unumpy.std_devs(x[up]), yerr = unumpy.std_devs(y[up]), color='b',label=r'$\mathrm{De}$ $\mathrm{Cia}$ $\mathrm{2016}$',fmt='o')
        axs.errorbar(unumpy.nominal_values(x[~up]), unumpy.nominal_values(y[~up]), xerr =  unumpy.std_devs(x[~up]), yerr = unumpy.std_devs(y[~up]), xuplims = xuplims[~up], uplims = yuplims[~up], color='b', markerfacecolor='None', markeredgecolor = 'b', fmt='o', alpha = 0.8, elinewidth=0.5)
    if(i == 3):
        ok = np.where(Wiseman_z == 3)
        x = Wiseman_OH[ok]
        y = Wiseman_DTM[ok]
        axs.errorbar(unumpy.nominal_values(x), unumpy.nominal_values(y), xerr =  unumpy.std_devs(x), yerr = unumpy.std_devs(y), label=r'$\mathrm{Wiseman}$ $\mathrm{2017}$', color='r',fmt='o')
        
        ok = np.where(Cia_z == 3)
        x = Cia_OH[ok]
        y = Cia_DTM[ok]
        xuplims = Cia_xuplims[ok]
        yuplims = Cia_yuplims[ok]
        up = np.logical_and(xuplims == 0, yuplims == 0)
        axs.errorbar(unumpy.nominal_values(x[up]), unumpy.nominal_values(y[up]), xerr =  unumpy.std_devs(x[up]), yerr = unumpy.std_devs(y[up]), color='b',label=r'$\mathrm{De}$ $\mathrm{Cia}$ $\mathrm{2016}$',fmt='o')
        axs.errorbar(unumpy.nominal_values(x[~up]), unumpy.nominal_values(y[~up]), xerr =  unumpy.std_devs(x[~up]), yerr = unumpy.std_devs(y[~up]), xuplims = xuplims[~up], uplims = yuplims[~up], color='b', markerfacecolor='None', markeredgecolor = 'b', fmt='o', alpha = 0.8, elinewidth=0.5)
    if(i == 4):
        ok = np.where(Wiseman_z == 4)
        x = Wiseman_OH[ok]
        y = Wiseman_DTM[ok]
        axs.errorbar(unumpy.nominal_values(x), unumpy.nominal_values(y), xerr =  unumpy.std_devs(x), yerr = unumpy.std_devs(y), label=r'$\mathrm{Wiseman}$ $\mathrm{2017}$', color='r',fmt='o')
    if(i == 5):
        ok = np.where(Wiseman_z == 5)
        x = Wiseman_OH[ok]
        y = Wiseman_DTM[ok]
        axs.errorbar(unumpy.nominal_values(x), unumpy.nominal_values(y), xerr =  unumpy.std_devs(x), yerr = unumpy.std_devs(y), label=r'$\mathrm{Wiseman}$ $\mathrm{2017}$', color='r',fmt='o')
		
def DTM_stell(axs, i):

    if(i == 0):
        x = RR_2015['SM']
        y = RR_2015['DTM_1B']
        yerror = np.abs((0.75*RR_2015['DTM_1B']))
        yall = unumpy.uarray(y, yerror)
        ynew = yall - unumpy.log10(1. + 10**yall)
        axs.errorbar(x, unumpy.nominal_values(ynew), yerr=unumpy.std_devs(ynew), color='g',label=r'$\mathrm{Remy-Ruyer}$ $\mathrm{2015}$',fmt='o') 

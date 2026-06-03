import os,sys,h5py
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import LineCollection
from scipy.constants import Boltzmann, Avogadro, R

### constants corner ###
RT          = R*298/1000
dx          = 0.1
step        = dx * RT
grid_points = np.arange(0,RT*5+step,step)
########################
dx          = 0.25
step        = dx
grid_points = np.arange(0,12+step,step)
########################


def Evaluate(FEL_data):
    thermal_limit = []
    for thermal_threshold in grid_points:
        mask              = (FEL_data <= thermal_threshold)
        low_energy_grid   = np.where(mask, FEL_data, np.nan)
        total_grid_points = low_energy_grid.shape[0] * low_energy_grid.shape[1]
        not_nan_mask      = ~np.isnan(low_energy_grid)
        count             = (np.count_nonzero(not_nan_mask) / total_grid_points) * 100
        thermal_limit.append(count)
    return thermal_limit

def Plot_Thermal_Units(ax,pos):
        skips = [0,6]
        for level in range(7):
            pos1  = [RT*level,RT*level]
            y     = [-1,101]
            ax[pos].plot(pos1,y,linestyle='-',color='grey',alpha=0.6,linewidth=0.25)
            ##################### Add text label
            if level not in skips:
                x     = pos1[0]+0.2
                y_fix = 95
                ax[pos].text(x,y_fix,f'{level} RT',color='grey')
            #####################

def Manage_Analysis(sequences,waters,systems,runs,colors_map,time='',label='',expt_label=''):

        df_columns=['Sequences','Run1','Run2','Run3','Mean','Std']

        if label == 'DYF':
            empty_array = np.zeros((len(sequences),len(df_columns)), dtype=object)

        for water in waters:

            j = 0
            fig,ax = plt.subplots(1,3,figsize=(12,3))
            Plot_Thermal_Units(ax,0)
            Plot_Thermal_Units(ax,1)
            Plot_Thermal_Units(ax,2)

            for k, mol in enumerate(sequences):

                for system in systems:

                    system_label     = f'{mol}{system}'
                    experiment_areas = [system_label]
                    print(system_label)
                    for run in runs:

                        if label == 'DYF' and system == '':
                            filepath = f'path/to/trajectory/folder/{seq}1_C36m-{water}/Run{run}/1{mol}_plots/FEL.csv'
                            
                        if os.path.exists(filepath):
                            FEL_data = np.genfromtxt(filepath,delimiter=',')
                            thermal_limit = Evaluate(FEL_data)
                            function_area = np.column_stack((grid_points,thermal_limit))   

                            ############################ integrate 
                            x_array = function_area[:,0]
                            y_array = function_area[:,1]
                            area_under_curve = round(np.trapz(y_array,x_array,dx=dx),0)
                            # print(area_under_curve,run,mol,system)
                            ##################################################
                            experiment_areas.append(area_under_curve)
                            title = f'{mol}_{system}_run{run}'
                            if run == '1':
                                ax[0].scatter(x_array,y_array,label=title,s=12,edgecolor='None')
                            elif run == '2':
                                ax[1].scatter(x_array,y_array,label=title,s=12,edgecolor='None')
                            elif run == '3':
                                ax[2].scatter(x_array,y_array,label=title,s=12,edgecolor='None')
                        else:

                            # print('No such file!')
                            experiment_areas.append(np.nan)

                    #########################################################
                    experiment_areas.append(round(np.nanmean(experiment_areas[1:]),0)) ## acounts for presence of np.nan in pandas df
                    experiment_areas.append(round(np.nanstd(experiment_areas[1:]),0))

                if label == 'DYF':
                    empty_array[k,:] = experiment_areas

            ax[0].legend(frameon=False,fontsize=8,loc='center left',handletextpad=0.05)
            ax[1].legend(frameon=False,fontsize=8,loc='center left',handletextpad=0.05)
            ax[2].legend(frameon=False,fontsize=8,loc='center left',handletextpad=0.05)

            ax[0].set_ylabel("τ",va='center',rotation='horizontal',fontsize=16)
            ax[1].set_xlabel('Energy (kJ/mol)',fontsize=16)
            plt.tight_layout()
            plt.subplots_adjust(left=0.05)
            plt.savefig(f'AUC_traces_{label}_{water}',dpi=800)
            plt.show()

            ##############################################################################
            df         = pd.DataFrame(empty_array,columns=df_columns)
            labels     = df['Sequences'].to_list()
            heights    = df['Mean']
            y_errors   = df['Std']
            colors     = [colors_map[label] for label in labels]
            x          = np.arange(1,len(labels)+1,1)
            plt_length = len(labels) * 2
            
            print(df)

            
            fig,ax = plt.subplots(figsize=(plt_length,3))

            if label == 'DYF':
                plt.close()
                fig,ax = plt.subplots(figsize=(plt_length,3))

            bars = plt.bar(x, heights, color=colors, yerr=y_errors, capsize=5,edgecolor='black',alpha=0.5)

            ax.set_ylabel(r'AUC ($\% \cdot \mathrm{kJ/mol}$)',fontsize=16)
            ax.set_xlabel('Sequences',fontsize=16)
            ax.set_xticks(x, labels)
            ax.set_ylim(0,500)
            ax.tick_params(axis='both',labelsize=12)
            plt.tight_layout()

            if label == 'DYF':
                output_name = f'DYF_{water}_AUC'
                df.to_csv(output_name+'.csv')
                plt.savefig(output_name,dpi=800,)

            plt.show()
            plt.close()
            ##############################################################################

colors_map = {'DFY':'cornflowerblue',
              'DYF':'royalblue',
              'FDY':'violet',
              'YDF':'orchid',
              'FYD':'turquoise',
              'YFD':'aquamarine'}

mols    = ['FDY']
waters  = ['sTIP3P']
systems = [''] 
runs    = ['1','2','3']

experiment_labels = ['D0']
water_models      = {'D0': ['sTIP3P','mTIP3P']} 

for expt_label in experiment_labels:

    water_models.get(expt_label)
    Manage_Analysis(mols,waters,systems,runs,colors_map,time='500ns',label='DYF',expt_label=expt_label)

sys.exit()

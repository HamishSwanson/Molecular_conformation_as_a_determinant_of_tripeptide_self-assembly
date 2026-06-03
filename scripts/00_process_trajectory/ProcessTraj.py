import os,sys,math,glob
import pandas as pd
import numpy as np
import shutil 
from pathlib import Path

single = True

if single:

    mols   = ['DYF','DFY','FDY','YDF','YFD','FYD']
    waters = ['mTIP3P','sTIP3P']
    concs  = ['1']
    runs   = ['1','2','3']
    time   = ['250ns','500ns']
    test   = True

    if test:
        mols   = ['FDY']        
        waters = ['sTIP3P']
        concs  = ['1']
        runs   = ['1','2','3']
        time   = ['500ns']

    for seq in mols:
        for water in waters:
            for conc in concs:
                for subdir in time:
                    for run in runs:
                        folder       = f'path/to/trajectory/folder/{seq}1_C36m-{water}/Run{run}/'
                        folder_above = f'path/to/trajectory/folder/{seq}1_C36m-{water}/'
                        simulation_output = f'{seq}{conc}_C36m-{water}_run{run}'
                        xtc               = f'{folder}{seq}{conc}_C36m-{water}_run{run}.xtc'
                        gro               = f'{folder}{seq}{conc}_C36m-{water}_run{run}.gro'
                        outfile           = f'{folder}{conc}_{seq}_C36m-{water}_run{run}_noPBC'
                        print(folder,simulation_output,'single')
                        
                        if os.path.exists(xtc) and os.path.exists(gro):
                            
                            print(f'XTC and GRO found!')

                            if not os.path.exists(f'{outfile}.xtc'):

                                cmd1 = f'echo protein protein | gmx trjconv -f {xtc} -s {folder_above}/*tpr -pbc mol -center -o {outfile}.xtc' 
                                cmd2 = f'echo protein protein | gmx trjconv -f {gro} -s {folder_above}/*tpr -pbc mol -center -o {outfile}.gro' 
                                os.system(cmd1)
                                os.system(cmd2)
                            
                            else:

                                print(f'\tDone >> {outfile}.xtc !')

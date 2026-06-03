import os,sys,math,re,h5py

default_expts   = True

if default_expts:

    mols   = ['DFY','DYF','FDY','FYD','YDF','YFD']
    waters = ['sTIP3P']
    runs   = ['1','2','3']
    time   = '500ns'
    test   = True

    if test:
        mols   = ['FDY']        
        waters = ['sTIP3P']
        runs   = ['1','2','3']
        time   = '500ns'

    for seq in mols:
        for water in waters:
            for run in runs:
                cmd1 = f'python SolveFEL.py {default_expts} {seq} {water} {run} {time}'
                print('Input>> ',cmd1)
                os.system(cmd1)

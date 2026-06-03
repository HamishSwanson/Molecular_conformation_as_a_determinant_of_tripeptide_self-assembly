import os,math,sys,h5py,re
import numpy as np
import pandas as pd
import MDAnalysis as mda
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from sklearn.metrics.pairwise import euclidean_distances
from scipy.constants import Boltzmann, Avogadro, R
from scipy.stats import gaussian_kde
from scipy import signal
from scipy import ndimage
import warnings
import config as global_config
from pathlib import Path
from multiprocessing import Pool, cpu_count

def generate_circular_ticks(start=180, step=90, max_angle=360):
    ticks = []
    current_angle = start
    while current_angle < max_angle:
        ticks.append(current_angle)
        current_angle += step
    if current_angle >= max_angle:
        for angle in range(0, 121, step):
            ticks.append(angle)
    return ticks
circular_ticks = generate_circular_ticks()

class LoadData:
    def __init__(self,data_file_path):
        self.file_path      = os.getcwd()
        self.data_file_path = data_file_path

    def read_hdf5_as_array(self):
        with h5py.File(self.data_file_path,'r') as hdf5_file:
            data = hdf5_file['dihedrals'][:]
        dihedrals = np.where(data<120,data+360,data) ## modify edge shifting!
        return dihedrals

class DataSaver:
    def __init__(self,conf_data,analysis_file_name):
        self.conf_data          = conf_data
        self.analysis_file_name = analysis_file_name

    def save_to_hdf5(self):
        file_info = f'{self.analysis_file_name}'
        with h5py.File(file_info, 'w') as hf:
            hf.create_dataset('EnergyLevels', data=self.conf_data)

class PopulationAnalysis:
    def __init__(self,data_file_path,h5_data,molNum,species,peptide_length,water):
        self.molNum          = molNum
        self.species         = species
        self.peptide_length  = peptide_length
        self.file_path       = os.getcwd()
        self.data_file_path  = data_file_path
        self.test_method     = global_config.TestMethod
        self.SampleWindow    = global_config.SampleWindow
        self.FirstDimension  = global_config.FirstDimension
        self.SecondDimension = global_config.SecondDimension
        self.num_processors  = global_config.num_processors
        self.minima_depth    = global_config.minima_depth
        self.grid_size       = 1
        self.torsional_dim   = np.arange(120,480,self.grid_size)
        self.h5_data         = h5_data
        self.water           = water

    def PlotTest1D(self,x,y,mol):
        torsions = [x,y]
        fig, axes = plt.subplots(nrows=2, ncols=1)
        for i,values in enumerate(torsions):
            axes[i].hist(values, bins=360, range=(120, 480), density=True,color='turquoise', edgecolor='black', alpha=0.7)
            axes[i].set_ylabel("Density", fontsize=8)
            axes[i].tick_params(axis='both', labelsize=8)
        axes[i].set_xlabel("Dihedral Angle (°)", fontsize=10)
        plt.tight_layout()
        plt.savefig(f'{self.h5_data}/{self.molNum}{self.species}_plots/{self.species}_{self.water}_X')
        plt.close()

    def PlotTest2D(self,x,y,mol):
        colors = {'DFY':'cornflowerblue','DYF':'royalblue','FDY':'violet','YDF':'orchid','FYD':'turquoise','YFD':'aquamarine'}
        fig    = plt.figure(figsize=(2,2))
        ###################################################
        tick_positions = [tick if tick >= 120 else tick + 360 for tick in circular_ticks]
        tick_labels    = [f'{tick}' for tick in circular_ticks]
        plt.xticks(tick_positions, tick_labels,fontsize=8)
        plt.yticks(tick_positions, tick_labels,fontsize=8)
        ###################################################
        color = colors.get(self.species)
        print(color)
        plt.scatter(x,y,color=color,s=0.5)
        plt.xlim(120,480)
        plt.ylim(120,480)
        plt.tight_layout()
        plt.savefig(f'{self.h5_data}/{self.molNum}{self.species}_plots/{self.species}_{self.water}_Y',dpi=800,transparent=True)
        plt.close()

    def PlotTest3D(self,x_grid,y_grid,z_grid,mol):
        fig           = plt.figure(figsize=(4,4))
        ax            = fig.add_subplot(111, projection='3d')
        surf          = ax.plot_surface(x_grid, y_grid, z_grid, cmap='Reds', edgecolor='k', alpha=0.8)
        ax.grid(True, color='black', linestyle='-', linewidth=0.5)
        label_library = {0:"Torsion 1-2 (°)",1:"Torsion 2-3 (°)",2:"Torsion 1-3 (°)"}
        ax.set_xlabel(label_library.get(self.FirstDimension), fontsize=10, labelpad=5)
        ax.set_ylabel(label_library.get(self.SecondDimension), fontsize=10, labelpad=5)
        ax.set_zlabel("Energy (kJ/mol)", fontsize=10,rotation=0)
        ax.view_init(elev=45, azim=-130) 
        ax.set_xlim(120,480)
        ax.set_ylim(120,480)
        ax.set_zlim(0,60)
        ax.set_zticks(np.arange(0,80,20))
        ###################################################
        tick_positions = [tick if tick >= 120 else tick + 360 for tick in circular_ticks]
        tick_labels    = [f'{tick}' for tick in circular_ticks]
        plt.xticks(tick_positions, tick_labels,fontsize=8)
        plt.yticks(tick_positions, tick_labels,fontsize=8)
        ###################################################
        # plt.suptitle(f'{self.molNum}{self.species}',fontsize=14)
        # plt.tight_layout()
        plt.subplots_adjust(bottom=0.02,top=0.98,left=0.12,right=0.98)
        plt.savefig(f'{self.h5_data}/{self.molNum}{self.species}_plots/{self.species}_{self.water}_Z')
        plt.close()

    def make_KDE_model(self,x,y,mol):
        combined_rows    = np.row_stack((x,y))
        self.PlotTest1D(x,y,mol)
        self.PlotTest2D(x,y,mol)
        grid_size        = 1 
        x_grid, y_grid   = np.meshgrid(self.torsional_dim,self.torsional_dim) 
        grid             = np.vstack([x_grid.ravel(), y_grid.ravel()])
        ################################################
        test_bandwidths  = {'250ns':0.327684,'500ns':0.291935}
        test_bandwidth   = test_bandwidths.get(time) #0.291935
        kde_model        = gaussian_kde(combined_rows,bw_method=test_bandwidth)
        ################################################ 
        pdf              = kde_model(grid)
        epsilon          = 1e-12                            ## fudge to remove where pdf = 0 giving np.log(pdf) = inf
        energy           = -R*298*np.log(pdf+epsilon)/1000
        energy           = energy - np.min(energy)
        z_grid           = energy.reshape(x_grid.shape)
        ################################################ 
        # max_energy = -R*298*np.log(epsilon)/1000
        # print(max_energy)
        ################################################ 
        print(f'\nSum of PDF: {np.sum(pdf)} > {mol}/{self.molNum}')
        v_unsampled   = 40 #np.max(z_grid)  # or use mode if there’s noise
        mask_sampled  = z_grid < v_unsampled
        z_grid        = np.where(mask_sampled, z_grid, np.nan)
        np.savetxt(f'{self.h5_data}/{self.molNum}{self.species}_plots/probabilities.csv',pdf,delimiter=',')
        np.savetxt(f'{self.h5_data}/{self.molNum}{self.species}_plots/FEL.csv',z_grid,delimiter=',')
        ################################################ 
        print('TestMethod!')
        self.PlotTest3D(x_grid, y_grid, z_grid,mol)    
        # if self.test_method and mol == 0:
        #     print('TestMethod!')
        #     self.PlotTest1D(x,y)
        #     self.PlotTest2D(x,y)
        #     self.PlotTest3D(x_grid, y_grid, z_grid)
        ################################################    
        local_minima               = ndimage.minimum_filter(z_grid,size=3,mode='wrap') == z_grid ### returns array where min in square grid of 9 are returned for that neighbourhood | then boolean to find where this == its z_grid value
        labeled_minima, num_minima = ndimage.label(local_minima)                     ### returns labels for each minima (1,...,n) and number of minima
        sorted_minima              = np.column_stack(np.where(local_minima))         ### location of minima within the z_grid stacked into an array
        ################################################ 
        coord_x = []
        coord_y = []
        energy  = []
        for x_value,y_value in sorted_minima:
            coord_x.append(x_grid[x_value, y_value])
            coord_y.append(y_grid[x_value, y_value])
            energy.append(z_grid[x_value, y_value])
        df = pd.DataFrame({
            'Coord_X': coord_x,
            'Coord_Y': coord_y,
            'Energy (kJ/mol)': energy})
        df_sorted = df.sort_values(by='Energy (kJ/mol)', ascending=True)
        top_n     = df_sorted.head(self.minima_depth).to_numpy()
        num_rows  = top_n.shape[0]
        num_cols  = top_n.shape[1]
        if num_rows < self.minima_depth:
            padded_array = np.full((self.minima_depth, num_cols), np.nan)
            padded_array[:num_rows, :] = top_n
        else:
            padded_array = top_n
        return padded_array

    def IterateMolecules(self):
        print(self.data_file_path)
        self.torsional_data = LoadData(self.data_file_path).read_hdf5_as_array()
        # results_array       = np.full((self.minima_depth,3,self.molNum), np.nan)
        results_array       = np.full((self.minima_depth,3), np.nan)
        for mol in range(self.molNum):
            x = self.torsional_data[self.SampleWindow:,self.FirstDimension,mol]
            y = self.torsional_data[self.SampleWindow:,self.SecondDimension,mol]
            FEL_results      = self.make_KDE_model(x,y,mol)
            FEL_results[:,0] = FEL_results[:,0]
            FEL_results[:,1] = FEL_results[:,1]
            np.set_printoptions(suppress=True, precision=5)
            # results_array[:,:,mol] = FEL_results
            results_array[:,:] = FEL_results
        print(results_array)
        return results_array
        
class FELAnalysis:
    def __init__(self,path,water):
        self.file_path        = os.getcwd()
        self.SampleWindow     = global_config.SampleWindow
        self.FirstDimension   = global_config.FirstDimension
        self.SecondDimension  = global_config.SecondDimension
        self.num_processors   = global_config.num_processors
        self.minima_depth     = global_config.minima_depth
        self.SampleWindowTime = global_config.SampleWindowTime
        self.species          = ""    ## default is unknown and inferred in process_file_pair() - will fail if not defined correctly
        self.molNum           = ""    ## default is unknown and inferred in process_file_pair() - will fail if not defined correctly
        self.peptide_length   = ""    ## default is unknown and inferred in process_file_pair() - will fail if not defined correctly
        self.path             = path
        self.water            = water
        # self.data_loader     = LoadData(data_source)

    def init_data(self):
        h5_data      = list(self.path.rglob('data*.h5'))
        self.h5_data = sorted(h5_data)

    def extract_simulation_info(self,data_file):
        pattern = r"data_(\d+)_(\w+)\.h5"
        print(data_file)
        match   = re.fullmatch(pattern,data_file)
        if match:
            molNum  = int(match.group(1))
            species = match.group(2)
        else:
            print('Discontinuity!')
        simulation_info = np.array((molNum,species,len(species)),dtype=object)
        return simulation_info

    def process_file_pair(self,data_file):
        h5_data             = data_file.parent
        h5_name             = data_file.name
        data_file_path      = str(data_file)
        simulation_info     = self.extract_simulation_info(h5_name)
        self.molNum         = simulation_info[0]
        self.species        = simulation_info[1].replace('p','') ## added 12/18/2025 for FFDp tripeptide 
        self.peptide_length = simulation_info[2]
        analysis_file_name  = f'{h5_data}/{self.molNum}{self.species}_{self.FirstDimension}-{self.SecondDimension}_{self.SampleWindowTime}ns_{self.minima_depth}levels.h5'
        if not os.path.exists(f'{h5_data}/{self.molNum}{self.species}_plots'):
            os.system(f'mkdir {h5_data}/{self.molNum}{self.species}_plots')
        if not os.path.exists(analysis_file_name+'test'):
            traj_analysis    = PopulationAnalysis(data_file_path,h5_data,self.molNum,self.species,self.peptide_length,self.water)
            analysis_results = traj_analysis.IterateMolecules()
            saver            = DataSaver(analysis_results,analysis_file_name) #+'.h5')
            hdf5_file        = saver.save_to_hdf5()
            np.savetxt(analysis_file_name.replace('.h5','.csv'), analysis_results, delimiter=',')
        else:
            print(f'Already Done: {h5_name}')
            pass

    def parallel_run(self): #,file_pair):
        data_file = list(self.h5_data)
        with Pool(processes=self.num_processors) as pool:
            pool.map(self.process_file_pair,data_file)

    def main(self):
        self.init_data()
        self.parallel_run()

class PlotTest:
    def __init__(self):
        self.file_path = os.getcwd()

    def generate_circular_ticks(start=180, step=90, max_angle=360):
        ticks = []
        current_angle = start
        while current_angle < max_angle:
            ticks.append(current_angle)
            current_angle += step
        if current_angle >= max_angle:
            for angle in range(0, 121, step):
                ticks.append(angle)
        return ticks

if __name__ == "__main__":
    TestMethod2      = True
    ###############################
    # SampleWindow     = 0  -- analyze the full trajectory    
    # FirstDimension   = 0  -- analyze 1-2 torsion
    # SecondDimension  = 1  -- analyze 2-3 torsion 
    # num_processors   = 1  -- use one CPU
    # minima_depth     = 10 -- save first 10 minimas found in FEL
    ###############################

    default_expts   = sys.argv[1] == "True"

    if default_expts:
        print('Default Experiments Selected')
        sequence  = sys.argv[2]
        water     = sys.argv[3]
        run       = sys.argv[4]
        time      = sys.argv[5]
        path      = f'path/to/trajectory/folder/{seq}1_C36m-{water}/'

    if not os.path.exists(path):
        print(f'{path} does not exist!')
        sys.exit()
    else:
        path_object = Path(path)
        analysis    = FELAnalysis(path_object,water)
        analysis.main()
        sys.exit()
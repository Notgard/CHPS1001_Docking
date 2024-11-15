POINT_SPACING=0.375
NB_POINT=126
GRID_SIZE = NB_POINT * POINT_SPACING
AUTODOCK_EXEC = '../autodock4'
AUTOGRID_EXEC = '../autogrid4'

# Importing the necessary libraries
import subprocess, multiprocessing
import plotly.graph_objects as go
import numpy as np
import random
import os

from datetime import datetime

# Lancer le docking : autodock4 -p galactose.dpf -l galactose.dlg & 
# Lancer la generation des cartes  : autogrid4 -p 3thc.gpf -l 3thc.glg & 

def get_center_boxes(molecule_dims):
    x_min, y_min, z_min = molecule_dims['xmin'], molecule_dims['ymin'], molecule_dims['zmin']
    x_size, y_size, z_size = molecule_dims['size']
    
    # Adding 2 arstorm units to min and max values to make sure the molecule is inside the boxes
    x_min -= 3
    y_min -= 3
    z_min -= 3
        
    x_step = x_size / 3
    y_step = y_size / 3
    z_step = z_size / 3
    
    centers = []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                center_x = x_min + (i * x_step) + (GRID_SIZE / 2)
                center_y = y_min + (j * y_step) + (GRID_SIZE / 2)
                center_z = z_min + (k * z_step) + (GRID_SIZE / 2)
    
                centers.append((center_x, center_y, center_z))
    return centers

def get_molecule_dimensions(pdb_file):
    with open(pdb_file, 'r') as file:
        x_coords = []
        y_coords = []
        z_coords = []

        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)

        x_min, x_max = min(x_coords), max(x_coords)
        y_min, y_max = min(y_coords), max(y_coords)
        z_min, z_max = min(z_coords), max(z_coords)

        dimensions = {
            'xmin': x_min,
            'ymin': y_min,
            'zmin': z_min,
            'xmax': x_max,
            'ymax': y_max,
            'zmax': z_max,
            'size': ((abs(x_max - x_min)), (abs(y_max - y_min)), (abs(z_max - z_min)))
        }

        return dimensions
    
def show_molecule_in_3d(pdb_file):
    # Show the molecule in 3D with plotly
    with open(pdb_file, 'r') as file:
        x_coords = []
        y_coords = []
        z_coords = []
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)
        
        fig = go.Figure(data=[go.Scatter3d(
            x=x_coords,
            y=y_coords,
            z=z_coords,
            mode='markers',
            marker=dict(
                size=2,
                color=z_coords,                # set color to an array/list of desired values
                colorscale='Viridis',   # choose a colorscale
                opacity=0.8
            )
        )])

        # tight layout
        fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
        fig.show()

def show_molecule_and_boxes_in_3d(pdb_file, centers):
    fig = go.Figure()

    # Add molecule points
    with open(pdb_file, 'r') as file:
        x_coords = []
        y_coords = []
        z_coords = []
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)
        
        fig.add_trace(go.Scatter3d(
            x=x_coords,
            y=y_coords,
            z=z_coords,
            mode='markers',
            marker=dict(
                size=2,
                color=z_coords,
                colorscale='Viridis',
                opacity=0.8
            ),
            name='Molecule'
        ))

    # Add boxes
    for center in centers:
        x_center, y_center, z_center = center
        # choose a different color for each box
        color = f'rgb({random.randint(0, 255)}, {random.randint(0, 255)}, {random.randint(0, 255)})'
        fig.add_trace(go.Mesh3d(
            x=[x_center - GRID_SIZE / 2, x_center + GRID_SIZE / 2, x_center + GRID_SIZE / 2, x_center - GRID_SIZE / 2, x_center - GRID_SIZE / 2, x_center + GRID_SIZE / 2, x_center + GRID_SIZE / 2, x_center - GRID_SIZE / 2],
            y=[y_center - GRID_SIZE / 2, y_center - GRID_SIZE / 2, y_center + GRID_SIZE / 2, y_center + GRID_SIZE / 2, y_center - GRID_SIZE / 2, y_center - GRID_SIZE / 2, y_center + GRID_SIZE / 2, y_center + GRID_SIZE / 2],
            z=[z_center - GRID_SIZE / 2, z_center - GRID_SIZE / 2, z_center - GRID_SIZE / 2, z_center - GRID_SIZE / 2, z_center + GRID_SIZE / 2, z_center + GRID_SIZE / 2, z_center + GRID_SIZE / 2, z_center + GRID_SIZE / 2],
            alphahull=0,
            opacity=0.05,
            color=color,
            name='Box'
        ))

    fig.update_layout(scene=dict(aspectmode='data'))
    fig.show()

def write_gpf_file(pdb_file, centers_boxes):
    box = 1
    gpf_files = []
    for center_coord in centers_boxes:
        pdb_file = pdb_file.split('/')[-1]
        x_center, y_center, z_center = center_coord
        # Keep only 4 decimal places
        x_center = round(x_center, 4)
        y_center = round(y_center, 4)
        z_center = round(z_center, 4)
        gpf_file = f"output/{pdb_file[:-4]}_box{box}.gpf"
        gpf_files.append(gpf_file.split('/')[-1])
        with open(gpf_file, 'w') as file:
            
            file.write(f"npts {NB_POINT} {NB_POINT} {NB_POINT}       # num.grid points in xyz\n")
            file.write(f"gridfld {pdb_file[:-4]}_box{box}.maps.fld            # grid_data_file\n")
            file.write(f"spacing {POINT_SPACING}                     # spacing(A)\n")
            file.write("receptor_types A C H HD N OA SA              # receptor atom types\n")
            file.write("ligand_types C HD OA                         # ligand atom types\n")
            file.write(f"receptor {pdb_file[:-4]}.pdbqt              # receptor pdbqt file\n")
            file.write(f"gridcenter {x_center} {y_center} {z_center} # xyz-coordinates or auto\n")
            file.write("smooth 0.5                                   # store minimum energy w/in rad(A)\n")
            file.write(f"map {pdb_file[:-4]}_box{box}.C.map                   # atom-specific affinity map\n")
            file.write(f"map {pdb_file[:-4]}_box{box}.HD.map                  # atom-specific affinity map\n")
            file.write(f"map {pdb_file[:-4]}_box{box}.OA.map                  # atom-specific affinity map\n")
            file.write(f"elecmap {pdb_file[:-4]}_box{box}.e.map               # electrostatic potential map\n")
            file.write(f"dsolvmap {pdb_file[:-4]}_box{box}.d.map              # desolvation potential map\n")
            file.write("dielectric -0.1465                           # <0, AD4 distance-dep.diel;>0, constant\n")
        box += 1
    return gpf_files

# Run autogrid4  
def run_autogrid(gpf_file):
    subprocess.run(cwd='output', args=[AUTOGRID_EXEC, '-p', gpf_file, '-l', f'{gpf_file[:-4]}.glg'])
def calculate_maps_of_molecule(gpf_files, nb_workers=5):    
    print(f"Started molecule maps calculation with {nb_workers} worker processes...")
    # Start asynchonous process with a pool of x workers
    pool = multiprocessing.Pool(nb_workers)
    pool.map(run_autogrid, gpf_files)
        
def write_dpf_receptor(nb_boxes, ligand, receptor, nb_ga_runs=50):
    dpf_files = []
    for i in range(nb_boxes):
        dpf_file = f"output/{receptor}_box{i+1}.dpf"
        dpf_files.append(dpf_file.split('/')[-1])
        ligand_name = f"{ligand}_box{i+1}"
        with open(dpf_file, 'w') as file:
            file.write("autodock_parameter_version 4.2       # used by autodock to validate parameter set\n")
            file.write("outlev 1                             # diagnostic output level\n")
            file.write("intelec                              # calculate internal electrostatics\n")
            file.write("seed pid time                        # seeds for random generator\n")
            file.write("ligand_types C HD OA                 # atoms types in ligand\n")
            file.write(f"fld {ligand_name}.maps.fld                    # grid_data_file\n")
            file.write(f"map {ligand_name}.C.map                       # atom-specific affinity map\n")
            file.write(f"map {ligand_name}.HD.map                      # atom-specific affinity map\n")
            file.write(f"map {ligand_name}.OA.map                      # atom-specific affinity map\n")
            file.write(f"elecmap {ligand_name}.e.map                   # electrostatics map\n")
            file.write(f"desolvmap {ligand_name}.d.map                 # desolvation map\n")
            file.write(f"move {receptor}.pdbqt                 # small molecule\n")
            file.write("about -23.8676 -7.2514 17.0284       # small molecule center\n")
            file.write("tran0 random                         # initial coordinates/A or random\n")
            file.write("quaternion0 random                   # initial orientation\n")
            file.write("dihe0 random                         # initial dihedrals (relative) or random\n")
            file.write("torsdof 6                            # torsional degrees of freedom\n")
            file.write("rmstol 2.0                           # cluster_tolerance/A\n")
            file.write("extnrg 1000.0                        # external grid energy\n")
            file.write("e0max 0.0 10000                      # max initial energy; max number of retries\n")
            file.write("ga_pop_size 150                      # number of individuals in population\n")
            file.write("ga_num_evals 250000                  # maximum number of energy evaluations\n")
            file.write("ga_num_generations 27000             # maximum number of generations\n")
            file.write("ga_elitism 1                         # number of top individuals to survive to next generation\n")
            file.write("ga_mutation_rate 0.02                # rate of gene mutation\n")
            file.write("ga_crossover_rate 0.8                # rate of crossover\n")
            file.write("ga_window_size 10                    # \n")
            file.write("ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution\n")
            file.write("ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution\n")
            file.write("set_ga                               # set the above parameters for GA or LGA\n")
            file.write("sw_max_its 300                       # iterations of Solis & Wets local search\n")
            file.write("sw_max_succ 4                        # consecutive successes before changing rho\n")
            file.write("sw_max_fail 4                        # consecutive failures before changing rho\n")
            file.write("sw_rho 1.0                           # size of local search space to sample\n")
            file.write("sw_lb_rho 0.01                       # lower bound on rho\n")
            file.write("ls_search_freq 0.06                  # probability of performing local search on individual\n")
            file.write("set_psw1                             # set the above pseudo-Solis & Wets parameters\n")
            file.write("unbound_model bound                  # state of unbound ligand\n")
            file.write(f"ga_run {nb_ga_runs}                            # do this many hybrid GA-LS runs\n")
            file.write("analysis                             # perform a ranked cluster analysis\n")
    return dpf_files

# Run autodock4
def run_autodock(dpf_file):
    subprocess.run(cwd='output', args=[AUTODOCK_EXEC, '-p', dpf_file, '-l', f'{dpf_file[:-4]}.dlg'])
def calculate_docking(dpf_files, nb_workers=5):
    print(f"Started docking calculation with {nb_workers} worker processes...")
    dlg_files = []
    # Start asynchonous process with a pool of x workers
    pool = multiprocessing.Pool(nb_workers)
    pool.map(run_autodock, dpf_files)
    for dpf_file in dpf_files:
        dlg_files.append(f'{dpf_file[:-4]}.dlg')
    return dlg_files

def extract_best_pose(dlg_file, output_pdb):
    with open("output/"+dlg_file, 'r') as dlg:
        lines = dlg.readlines()

    best_score = float('inf')
    best_pose_start = None
    current_pose_start = None
    current_score = None

    for i, line in enumerate(lines):
        if 'DOCKED: MODEL' in line:
            current_pose_start = i
        elif 'DOCKED: USER    Estimated Free Energy of Binding' in line:
            current_score = float(line.split()[-3])
            if current_score < best_score:
                best_score = current_score
                best_pose_start = current_pose_start

    if best_pose_start is not None:
        with open(output_pdb, 'w') as out_pdb:
            for line in lines[best_pose_start:]:
                if line.startswith('DOCKED: ENDMDL'):
                    break
                if line.startswith('DOCKED: ATOM') or line.startswith('DOCKED: HETATM'):
                    out_pdb.write(line.replace('DOCKED: ', ''))
    return output_pdb, best_score 
def extract_best_poses(dlg_files):
    # Create the output folder if it doesn't exist
    if not os.path.exists("output/best_poses"):
        os.makedirs("output/best_poses")
    best_poses = []
    for dlg_file in dlg_files:
        output_pdb = "output/best_poses/"+dlg_file.replace('.dlg', '_best.pdb')
        result = extract_best_pose(dlg_file, output_pdb)
        best_poses.append(result)
    return best_poses
def extract_best_pose_of_best_poses(best_poses):
    best_score = float('inf')
    best_pose = None
    for result in best_poses:
        score = result[1]
        pose = result[0]
        if score < best_score:
            best_score = score
            best_pose = pose
    return best_pose, best_score

def clear_output():
    os.system("rm -rf output/*box*")
    
def store_output_files():
    date = datetime.now()
    timestamp = date.strftime("%m_%d_%Y_%H_%M_%S")
    if not os.path.exists(f"output/{timestamp}"):
        os.makedirs(f"output/{timestamp}")
    os.system(f"mv output/*box* -t output/{timestamp}")
    os.system(f"cp -r output/best_poses -t output/{timestamp}")


# Necessary files :
#   ligand.pdb
#   receptor.pdb
#   ligand.pdbqt
#   receptor.pdbqt
# Fix the CONSTANTS AUTODOCK_EXEC and AUTOGRID_EXEC with the path to the executables

# Usage example
receptor = 'galactose'
ligand = '3thc'
pdb_file = f'output/{ligand}.pdb'
clear_output()

dimensions = get_molecule_dimensions(pdb_file)
grid_centers = get_center_boxes(dimensions)
print(f"Number of boxes: {len(grid_centers)}")

show_molecule_in_3d(pdb_file) # Active to check the 3D representation
show_molecule_and_boxes_in_3d(pdb_file, grid_centers) # Active to check the 3D representation
"""
print("Generating the gpf files... ",end="")
gpf_files = write_gpf_file(pdb_file, grid_centers)

print("done.\nCalculating the maps...",end="")
calculate_maps_of_molecule(gpf_files, 27)

print("done.\nGenerating the dpf files...",end="")
dpf_files = write_dpf_receptor(len(grid_centers), ligand, receptor)

print("done.\nCalculating the docking...",end="")
dlg_files = calculate_docking(dpf_files, 27)

print("done.\nExtracting the best poses...",end="")
best_poses = extract_best_poses(dlg_files)

print("done.\nExtracting the best pose of best poses...",end="")
best_pose = extract_best_pose_of_best_poses(best_poses)
print("done.\nBest pose:", best_pose[0], "score:", best_pose[1])

store_output_files()
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 08:55:44 2022

@author: afisher
"""

from struct import unpack
import pandas as pd

def aggregate_gdf(files, skip_particle=False):
    ''' Compute aggregations on the particle data. Combine with the params
    dataframe into a single dataframe where the index is timestep. If file 
    is a list of files, concatenate the dataframes ensuring no overlap in avgz.'''
    
    if not isinstance(files, list):
        files = [files]
    
    gdfs = []
    particles = []
    for file in files:
        _, params, particle = load_gdf(file)

        # Compute aggregations on particle data
        grouped = particle.groupby(level=0)
        counts = grouped.z.count().rename('counts')
        means = grouped.mean().add_prefix('avg')
        stds = grouped.std().add_prefix('std')
        gdf = pd.concat([params, counts, means, stds], axis=1)
        
        # Remove screen outputs (where t is defined) and append
        gdf = gdf.loc[gdf.avgt.isnull()].dropna(axis=1, how='all')
        gdfs.append(gdf)
        
        # Do same for particles, optionally skip
        if not skip_particle:
            particle = particle[particle.t.isnull()].dropna(axis=1, how='all')
            particles.append(particle)
    

    # Concatenate with our without particles
    if len(particles)==0:
        return concatenate_gdfs(gdfs)
    else:
        return concatenate_gdfs(gdfs, particles)
    

def concatenate_gdfs(gdfs, particles=[]):
    ''' Concatenates gdfs. Truncate so no overlaps in avgz. '''
    # Remove overlaps in avgz
    for j in range(len(gdfs)-1):
        overlap = gdfs[j].avgz >= gdfs[j+1].avgz[0]
        gdfs[j] = gdfs[j][~overlap]
        if len(particles)>0:
            p_idx = particles[j].index.get_level_values(level=0).unique()[~overlap]
            particles[j] = particles[j].loc[(p_idx, slice(None))]
            
    # Adjust indices for concatenation
    idx_start = 0
    for j,gdf in enumerate(gdfs):
        new_index = range(idx_start, idx_start+len(gdf))
        idx_start += len(gdf)

        gdf.index = new_index
        if len(particles)>0:
            particles[j].index = particles[j].index.set_levels(new_index, level=0)
            
            
    # Concatenate and return
    if len(particles)==0:
        gdf_full = pd.concat(gdfs).reset_index(drop=True)
        return gdf_full
    else:
        gdf_full = pd.concat(gdfs).reset_index(drop=True)
        particle_full = pd.concat(particles)
        return gdf_full, particle_full


def load_gdf(filename, arrays_to_load=[['x','y','z','Bx','By','Bz','G','t']]):
    ''' GDF files contain parameter and array data for each timestep. If 
    arrays_to_load=[], load_gdf returns non-varying parameters as a 
    Series object and varying parameters are a DataFrame object. arrays_to_load 
    specifies a list of lists containing equal length arrays to group into 
    DataFrames (all arrays in a DataFrame must have same length). For example, 
    arrays_to_load=[['z','G'], ['freq00','A00']] will return two dataframes in 
    addition to the constant series object and varying parameter dataframe. 
    The first contains the particle data 'z' and 'G' while the second contains 
    the fel mode data 'freq00' and 'A00'.
    
    Parameter dataframe structure: 
        Index: Timestep
        Columns: Parameter Name
        
    Array dataframe structure:
        MultiIndex: 
            Level0: Timestep
            Level1: Array elements (particles or modes from example)
        Columns: Array Name
    '''
    
    encoding = 'windows-1252'
    arrays_dict = {ele:iarr 
                   for iarr, arr in enumerate(arrays_to_load) 
                   for ele in arr }
    
    # Create containers for data
    param_cnt = []
    arrays_cnt = [ [] for _ in range(len(arrays_to_load))]
    arrays = [ {} for _ in range(len(arrays_to_load))]
    
    info, params = {}, {}
    params['level1'] = {}
    elements_stored = 0
    level = 1


    GDFNAMELEN = 16
    GDFID = 94325877 
    t_ascii = int('0001', 16)
    t_s32 = int('0002', 16)
    t_dbl = int('0003', 16)
    t_nul = int('0010', 16)

    with open(filename, 'rb') as f:
        
        # Read file header
        info['ID'] = int.from_bytes(f.read(4), 'little')
        info['cretime'] = int.from_bytes(f.read(4), 'little')
        info['creator'] = f.read(GDFNAMELEN).decode(encoding).rstrip('\x00')
        info['destin'] = f.read(GDFNAMELEN).decode(encoding).rstrip('\x00')
        info['gdfmaj'] = int.from_bytes(f.read(1), 'little')
        info['gdfmin'] = int.from_bytes(f.read(1), 'little')
        info['cremaj'] = int.from_bytes(f.read(1), 'little')
        info['cremin'] = int.from_bytes(f.read(1), 'little')
        info['desmaj'] = int.from_bytes(f.read(1), 'little')
        info['desmin'] = int.from_bytes(f.read(1), 'little')
        info['dummy'] = int.from_bytes(f.read(2), 'little')
    
        # Read data blocks
        while True:
            name_in_bytes = f.read(GDFNAMELEN)
            name = name_in_bytes.decode(encoding).rstrip('\x00')

            # End loop if no bytes left
            if name_in_bytes == b'':
                break
                
            # Parse info in block header
            block_type = bin(int.from_bytes(f.read(4), 'little'))[2:].zfill(32) #binary
            block_size = int.from_bytes(f.read(4), 'little')
            data_type = int(block_type[-8:], 2)
            start_dir = (block_type[2*8 + 7]=='1')
            end_dir = (block_type[2*8 + 6]=='1')
            data_is_param = (block_type[2*8 + 5]=='1')
            data_is_array = (block_type[2*8 + 4]=='1')
            
            # If start of directory, increase level
            if start_dir:
                level +=1
                params[f'level{level}'] = params[f'level{level-1}']
            
            # If end of directory, append data and reduce level
            if end_dir:
                for array_cnt, array in zip(arrays_cnt, arrays):
                    array_cnt.append(pd.DataFrame(array))
                param_cnt.append( params[f'level{level}'].copy() )
                elements_stored +=1
                arrays = [ {} for _ in range(len(arrays_to_load))]
                level -= 1
            
            # When data is a parameter
            if data_is_param:
                if data_type==t_dbl:
                    value = unpack('d', f.read(8))[0]
                    params[f'level{level}'][name] = value
                elif data_type==t_nul:
                    pass
                elif data_type==t_ascii:
                    value = f.read(block_size).decode(encoding)
                    params[f'level{level}'][name] = value 
                elif data_type==t_s32:
                    value = int.from_bytes(f.read(4), 'little')
                    params[f'level{level}'][name] = value 
                else:
                    f.close()
                    print('Error: Unknown datatype.')
            
            # When data is an array
            if data_is_array:
                if data_type == t_dbl:
                    if block_size%8 != 0:
                        f.close()
                        print('Error: wrong size for double array')
                    N = int(block_size/8)
                    values = unpack('d'*N, f.read(8*N))
                    
                    if name in arrays_dict.keys():
                        arrays[arrays_dict[name]][name] = values
                    elif len(arrays_dict)==0:
                        arrays[name] = values
                else:
                    f.close()
                    print('Error: Datatype not double')
                    
    # Create dataframes from containers
    array_dfs = [pd.concat(array_cnt, keys=range(len(array_cnt))) 
                 for array_cnt in arrays_cnt]
    param_df = pd.DataFrame(param_cnt)
    
    # Separate constants from param_df
    constant_tf = (param_df.iloc[0]==param_df).all()
    constants = param_df.columns[constant_tf]
    const_series = param_df.loc[0, constants]
    param_df.drop(columns=constants, inplace=True)
    
    # Add cputime and numderivs to const_series
    const_series['cputime'] = params['level1']['cputime']
    const_series['numderivs'] = params['level1']['numderivs']
    
    return const_series, param_df, *array_dfs






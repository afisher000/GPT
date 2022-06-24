# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 08:55:44 2022

@author: afisher
"""

from struct import unpack
import pandas as pd


def load_gdf(filename, arrays_to_load, elements_to_load=None):
    encoding = 'windows-1252'
    arrays_dict = {ele:iarr 
                   for iarr,arr in enumerate(arrays_to_load) 
                   for ele in arr }
    last_element = -1 if elements_to_load is None else elements_to_load[-1]
    param_df = []
    arrays = [ {} for _ in range(len(arrays_to_load))]
    array_dfs = [ [] for _ in range(len(arrays_to_load))]
    
    info, params = {}, {}
    params['level1'] = {}
    current_element = 1
    elements_stored = 0
    level = 1

    # Using little endian
    GDFNAMELEN = 16
    GDFID = 94325877 
    
    t_ascii = int('0001', 16)
    t_s32 = int('0002', 16)
    t_dbl = int('0003', 16)
    t_nul = int('0010', 16)

    with open(filename, 'rb') as f:
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
    
        # Parse data
        while True:
            name_in_bytes = f.read(GDFNAMELEN)
            name = name_in_bytes.decode(encoding).rstrip('\x00')
            if name_in_bytes == b'':
                break
                
            block_type = bin(int.from_bytes(f.read(4), 'little'))[2:].zfill(32) #binary
            block_size = int.from_bytes(f.read(4), 'little')
            
            data_type = int(block_type[-8:], 2)
            start_dir = (block_type[2*8 + 7]=='1')
            end_dir = (block_type[2*8 + 6]=='1')
            data_is_param = (block_type[2*8 + 5]=='1')
            data_is_array = (block_type[2*8 + 4]=='1')
            
            if start_dir:
                level +=1
                params[f'level{level}'] = params[f'level{level-1}']
            
            if end_dir:
                if len(arrays)==0:
                    if elements_to_load is None:
                        elements_stored +=1
                        param_df.append( params[f'level{level}'].copy() )
                else:
                    if (elements_to_load is None) or (current_element in elements_to_load):
                        for df, array in zip(array_dfs, arrays):
                            df.append(pd.DataFrame(array))
                        param_df.append( params[f'level{level}'].copy() )
                        elements_stored +=1
                        if current_element == last_element:
                            break
                    arrays = [ {} for _ in range(len(arrays_to_load))]
                    current_element = current_element + 1
                level -= 1
            
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
            
            if data_is_array:
                if data_type == t_dbl:
                    if block_size%8 != 0:
                        f.close()
                        print('Error: wrong size for double array')
                    N = int(block_size/8)
                    values = unpack('d'*N, f.read(8*N))
                    
                    if name in arrays_dict.keys():
                        arrays[arrays_dict[name]][name] = values
                else:
                    f.close()
                    print('Error: Datatype not double')
                    
    # Create dataframes
    array_dfs = [pd.concat(df, keys=range(len(df))) for df in array_dfs]
    param_df = pd.DataFrame(param_df)
    return param_df, *array_dfs






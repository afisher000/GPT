# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 08:55:44 2022

@author: afisher

% Load binary gdf file from GPT in matlab structure array
%
% Input 
%       filename    File to load
%       arrays_to_load      (Optional) arrays to load from the gdf file,
%                   specify as cell array {'name1';'name2'} etc
%                   For example {'x';'Bx'}
%                   All arrays are loaded when this option is not given,
%                   or when '' is passed. 
%       elements_to_load    (Optional) elements to load from the gdf file, 
%                   specified as an array. For example when you have a gdf
%                   file with particle output at t=0 and t=1 as a function
%                   of a scanned parameters, and you are only interested in
%                   t=1 for now, specify elements = [2 4 6 ...]
%                   When you are using this option and want to load all
%                   array elements (see above), pass '' as the second
%                   argument.
%                   When not given, all elements are loaded
%
% Output
%       gdf         1 x N structure array with data with N elements and 
%                   two fields for each element:
%                   .d  structure with data in single group. Fields:
%                     ID  Particle identification numbers
%                     x   x coordinate [m]
%                     y   y coordinate [m]
%                     z   z coordinate [m]
%                     Bx  Normalized velocity beta_x = v_x / c
%                     By  Normalized velocity beta_y = v_y / c
%                     Bz  Normalized velocity beta_z = v_y / c
%                     G   Lorentz factor
%                     rxy Distance to z-axis [m], sqrt(x^2+y^2)
%                     When using screen output command in GPT:
%                     t   Time when particles cross the screen [s] 
%                     When using tout (time output) command in GPT:
%                     fEx Ex field at the particle coordinates [V/m]
%                     fEy Ey field at the particle coordinates [V/m]
%                     fEz Ez field at the particle coordinates [V/m]
%                     fBx Bx field at the particle coordinates [T]
%                     fBy By field at the particle coordinates [T]
%                     fBz Bz field at the particle coordinates [T]
%                   .p  structure with params of that group
%       info        structured array with main header information
%
% Example
%       load_gdf('data.gdf')           %Load whole file
%       load_gdf('data.gdf',{'x';'y'}) %Load only x and y data
%       load_gdf('data.gdf','', [2 4]) %Load only elements 2 and 4
%
% Structure gdf file: (each line contains a block of data)
% Header (level = 1 after this step)
% param @logo
% params scanned by mr (level = 2 after this step)
% param creator
% param @logo
% param position / time (level = 3 after this step)
% dat particle data (level = 2 after this step)
% param position / time (possibly) (level = 3 after this step)
% dat particle data     (possibly) (level = 2 after this step)
% param numderivs (only loaded into structure when elements_to_load is not specified as input)
% param cputime   (only loaded into structure when elements_to_load is not specified as input)
% (level = 1 after this step). All data belonging to these mr settings are read
% params scanned by mr
% param @logo
% param position / time
% dat particle data
% param position / time (possibly)
% dat particle data     (possibly)
% param numderivs (only loaded into structure when elements_to_load is not specified as input)
% param cputime   (only loaded into structure when elements_to_load is not specified as input)
% etc
"""

from struct import unpack
import pandas as pd

# Using little endian (first byte is least significant)
GDFNAMELEN = 16
GDFID = 94325877 

t_ascii = int('0001', 16)
t_s32 = int('0002', 16)
t_dbl = int('0003', 16)
t_nul = int('0010', 16)

t_dir = int('0100', 16)
t_edir = int('0200', 16)
t_param = int('0400', 16)
t_data = int('0800', 16)

filename = 'test.gdf'
encoding = 'windows-1252'

gdfd = []
gdfp = []
info, params, arrays = {}, {}, {}
params['level1'] = {}
# May need to use n-1 for these converting matlab to python indexing
current_element = 1
elements_stored = 0
level = 1

# elements_to_load and arrays_to_load are inputs
elements_to_load = []
arrays_to_load = []
last_element_to_load = -1


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
        byte_name = f.read(GDFNAMELEN)
        if byte_name == b'':
            break

        name = byte_name.decode(encoding).rstrip('\x00')
            
        block_type = bin(int.from_bytes(f.read(4), 'little'))[2:].zfill(32) #binary
        block_size = int.from_bytes(f.read(4), 'little')
        print(block_size)
        
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
                if len(elements_to_load)==0:
                    if isinstance(arrays_to_load, str):
                        elements_stored +=1
                    gdfp.append(params[f'level{level}'])
            else:
                if len(elements_to_load)==0 or (current_element in elements_to_load):
                    gdfd.append(arrays)
                    gdfp.append(params[f'level{level}'])
                    elements_stored +=1
                    if current_element == last_element_to_load:
                        break
                arrays = {}
                current_element = current_element+1
                print(current_element)
            level -= 1
        
        if data_is_param:
            if data_type==t_dbl:
                value = int.from_bytes(f.read(8), 'little')
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
                
                if len(arrays_to_load)==0 or (name in arrays_to_load):
                    arrays[name] = values #might break if non-allowed char exist
            else:
                f.close()
                print('Error: Datatype not double')
                






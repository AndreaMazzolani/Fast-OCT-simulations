function out_file = alloc_file_matfile(filename, dsetname, max_size_pix, voxel_size_metres, dset_type, varargin)
%% out_file = alloc_file_matfile(filename, dsetname, max_size_pix, voxel_size_metres, dset_type)
% out_file = alloc_file_matfile(filename, dsetname, max_size_pix, voxel_size_metres, dset_type, 'complex', true|false)
%
% Creates a new *.mat file to contain data to be written/read in chunks.
%
% Uses the 'matfile' class present in Matlab >R2011b to allocate data
% storage on disk, and read/write access the data once it's created.
% Note that this means it uses HDF5 as the underlying storage mechanism.
%
% INPUT:
% ------
% filename: STRING - Path to the *.mat file that will be used to store
% the OCT data. NOTE! if the file already exists, *it will be COMPLETELY
% overwritten!*
%
% dsetname: STRING or CELL ARRAY of STRINGS - name (or names) of the
% dataset(s) to be stored in this file.
%
% max_size_pix: STRUCTURE or ARRAY of STRUCTURES with .x, .y, .z
% members - positive integers. Stores the maximum number of pixels that
% will be present in (each) dataset.
%
% voxel_size_metres: STRUCTURE with .x, .y, .z members - positive real
% numbers. Stores the optical size in metres of the individual voxels in
% (each) dataset. NOTE!! ASSUMES all datasets stored in this file have
% the same voxel sizes!
%
% dset_type: STRING or CELL ARRAY of STRINGS - the matlab classtype for
% (each) of the datasets, e.g., 'single', 'double' etc.
%
% 'complex': (OPTIONAL) LOGICAL or ARRAY of LOGICALs - if 'true', then
% the corresponding dataset will be allocated to hold complex numbers.
% Defaults to 'false'.
%
% OUTPUT:
% ------
% out_file: a writable 'matfile' object allowing access to the OCT data
% and associated metadata stored in the *.mat file.
%
% out_file.voxel_size_metres: STRUCTURE with .x, .y, .z members - copied
% from voxel_size_metres.
%
% out_file."dsetname": 3-D array(s) holding the data.
%
% EXAMPLE:
% --------
% Allocate a .mat file to hold two datasets, 'Dset1', single precision
% complex floating point, and 'Dset2', double precision reals.
%
% dsetname = {'Dset1', 'Dset2'};
% max_size_pix(1) = struct('x', 100, 'y', 100, 'z', 50);
% max_size_pix(2) = struct('x', 100, 'y', 100, 'z', 40);
% voxel_size_metres = struct('x', 10e-6, 'y', 10e-6, 'z', 5e-6);
% dset_type = {'single', 'double'};
% cplx = [true, false];
% out_file = alloc_file_matfile('testfile.mat', dsetname, ...
%                               max_size_pix, voxel_size_metres, ...
%                               dset_type, 'complex', cplx);
% ...
% size_set1 = size(out_file, 'Dset1'); % [100, 100, 50]
% size_set2 = size(out_file, 'Dset2'); % [100, 100, 40]
% dataset1_enface = out_file.Dset1(:, :, 10);
% dataset2_enface = out_file.Dset2(:, :, 10);

%% Validate the inputs
input_isvalid = false;
dsetname_valid = @(x) ischar(x) || (iscell(x) && ischar([x{:}]));

% Operate on a possible array of inputs (must be finite +ve integers)
pix_size_valid = @(x) all(isfinite(x) & ( x == floor(x) ) & x > 0);
% NB: concatenate the subfields of max_size_pix into arrays before
% validating them through pix_size_valid, otherwise p_s_v(s.x) -> e.g.,
% p_s_v(1, 2, 3), which is invalid and causes an error
max_size_valid = @(s) all(isfield(s, {'x', 'y', 'z'})) && ...
    pix_size_valid([s.x]) && pix_size_valid([s.y]) && pix_size_valid([s.z]);

% No need to do this with the voxel sizes, since there should be only
% one voxel structure for the either dataset
vox_size_valid = @(x) all(isfinite(x) & x > 0);
vox_size_m_valid = @(s) all(isfield(s, {'x', 'y', 'z'})) && ...
    vox_size_valid(s.x) && vox_size_valid(s.y) && vox_size_valid(s.z);

dset_type_valid = @(x) ischar(x) || (iscell(x) && ischar([x{:}]));

p = inputParser();
addRequired(p, 'filename', @ischar);
addRequired(p, 'dsetname', dsetname_valid);
addRequired(p, 'max_size_pix', max_size_valid);
addRequired(p, 'voxel_size_metres', vox_size_m_valid);
addRequired(p, 'dset_type', dset_type_valid);
p.addParamValue('complex', false(size(max_size_pix)), @islogical);
p.parse(filename, dsetname, max_size_pix, voxel_size_metres, dset_type, ...
        varargin{:});

if ischar(dsetname) && ischar(dset_type) && numel(max_size_pix) == 1 ...
        && numel(p.Results.complex) == 1
    input_isvalid = true;
    num_dsets = 1;
elseif iscell(dsetname) && iscell(dset_type) && ...
        (numel(dsetname) == numel(dset_type)) && ...
        (numel(dsetname) == numel(max_size_pix)) && ...
        (numel(dsetname) == numel(p.Results.complex))
    input_isvalid = true;
    num_dsets = numel(dsetname);
end

if ~input_isvalid
    error('%s: Error in inputs: number of dataset parameters (dsetname, max_size_pix, dset_type, complex) must match!', mfilename);
end

%% Delete any existing file
if exist(filename, 'file')
    delete(filename);
end

%% Open a new matfile for writing
% Adds the attributes etc. that matlab requires to recognise this hdf5
% file as a .mat file
out_file = matfile(filename, 'Writable', true);

%% Save the metadata
% Save a structure holding the optical size in metres of each voxel in
% the 3-D OCT volume
vsm = struct('x', voxel_size_metres.x, 'y', voxel_size_metres.y, 'z', ...
             voxel_size_metres.z);
out_file.voxel_size_metres = vsm;

%% Pre-allocate the output volumes
% The matfile interface doesn't allow for setting of the chunk_size, a
% very important parameter that determines the performance of read/write
% access to the array, so we need to bypass the matfile layer and
% allocate the dataset(s) using low level HDF5 calls
%
% close the matfile layer first
clear out_file;
% allocate the datasets using the low level HDF5 interface
import OCT_OCE.FileIO.*;

if num_dsets == 1
    % Allocate room for a single dataset
%     fprintf('%s: Alloc: %s, [%d, %d, %d], type %s, cplx: %d\n', mfilename,...
%             dsetname, max_size_pix.x, max_size_pix.y, max_size_pix.z, ...
%             dset_type, p.Results.complex);
    alloc_matlab_hdf5_dataset(...
        filename, dsetname, ...
        [max_size_pix.x, max_size_pix.y, max_size_pix.z], ...
        dset_type, 'complex', p.Results.complex);
else
    % Allocate room for multiple datasets
    for idx = 1:num_dsets
%         fprintf('%s: %02d: Alloc: %s, [%d, %d, %d], type %s, cplx: %d\n', mfilename,...
%                 idx, ...
%                 dsetname{idx}, max_size_pix(1).x, max_size_pix(1).y, max_size_pix(1).z, ...
%                 dset_type{idx}, p.Results.complex(idx));
        alloc_matlab_hdf5_dataset(...
            filename, dsetname{idx}, ...
            [max_size_pix(idx).x, max_size_pix(idx).y, max_size_pix(idx).z], ...
            dset_type{idx}, 'complex', p.Results.complex(idx));
    end
end
% Reopen the file using the matlab interface for writing
out_file = matfile(filename, 'Writable', true);

function alloc_matlab_hdf5_dataset(filename, dset_name, size_pix, type, varargin)
%% alloc_matlab_hdf5_dataset(filename, dset_name, size_pix, type)
%% alloc_matlab_hdf5_dataset(filename, dset_name, size_pix, type, 'complex', true)
% Allocates an HDF5 dataset of name given by 'dset_name' in the HDF5
% file given by 'filename', of a size given by 'size_pix' (an n-d array
% of integers), with matlab class type given by 'type' (a string
% describing the matlab datatype, e.g., 'single', 'double'). Optionally,
% this can be a complex datatype, enabled by passing in the optional
% parameter 'complex' to true.
%
% The resulting dataset is designed to be readable by Matlab as a
% variable in a .mat file.
%
% NOTE! The hdf5 file *must be closed* before calling this function!
%
% NOTE! Only tested for simple arrays of numbers!

p = inputParser();
p.addParamValue('complex', false, @islogical);
p.parse(varargin{:});

%% Create a dataset manually in order to set the chunksize properly
fid = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');

%% type_id: derived from the matlab datatype
[type_id, type_size_bytes, type_fill] = matlab_to_hdf5_type(type, p.Results.complex);

%% space_id: the current and maxsize of the dataset.
% Note that HDF5 arrays are in reverse order as compared to matlab
% arrays.
dims = size_pix;
h5_dims = fliplr(dims);
% Match the Matlab default of setting the maxsize to UNLIMITED
h5_maxdims = repmat(H5ML.get_constant_value('H5S_UNLIMITED'), size(h5_dims));
space_id = H5S.create_simple(3, h5_dims, h5_maxdims);

%% dcpl_id: set dataset parameters, especially the chunk size
dcpl_id = H5P.create('H5P_DATASET_CREATE');
% Set the default value of uninitialised data
H5P.set_fill_value(dcpl_id, type_id, type_fill);
% set chunk size to roughly 64k bytes (arbitrary threshold), i.e.,
% product of chunk_dims * element size ~ 65536 bytes.
%
% prod(chunk_dims) * type_size_bytes = chunk_max_bytes;
% => prod(chunk_dims) = (chunk_max_bytes; / type_size_bytes) = chunk_size
% => prod( dims ./ dimen_factor ) = chunk_size
% => prod( dims ) ./ (dimen_factor ^ num_dims) = chunk_size
% => prod( dims ) ./ chunk_size = dimen_factor ^ num_dims
% => dimen_factor = (prod(dims) ./ chunk_size) ^ (1/num_dims)
dimen_factor = round( nthroot( prod(dims) / (65536/type_size_bytes), numel(dims))); 


% NOTE: The chunksize is set to an isotropic fraction of the total array
% size, this is to allow for (relatively) quick reslicing of the data
% along any axes
chunk_dims = round(dims/dimen_factor); % wallen 20170309 changed from round
chunk_dims = max(chunk_dims, ones(size(chunk_dims)));
h5_chunk_dims = fliplr(chunk_dims);
H5P.set_chunk(dcpl_id, h5_chunk_dims);
% NOTE: No compression! When saving individual B-scans to a 3-D dataset,
% compression slows down saving by a factor of >100, whilst improving
% the final file size only by ~10-20%
%
% set compression
% Matlab defaults level 3 deflate (out of 9)
% zlib (standard compression library) defaults to level 6 deflate
%H5P.set_shuffle(dcpl_id); % enhances the effectiveness of compression
%H5P.set_deflate(dcpl_id, 3);

% For boolean data, enable HDF5's nbit compression, since there's only
% 1-bit of actual data in each 8-bit data-type.
if strcmpi(type, 'logical')
    H5P.set_nbit(dcpl_id);
end

%% Create the dataset
dset_id = H5D.create(fid, dset_name, type_id, space_id, dcpl_id);

%% attr_id: Set attributes for Matlab interpretation
% Matlab requires an attribute called 'MATLAB_class', set to a string
% corresponding to the matlab class of the dataset
attrtype_id = H5T.copy('H5T_C_S1');
% BUG! The matlab hdf5 low level interface doesn't allow for variable
% length strings in attributes (although it implies it does, Matlab
% 2012b)
%H5T.set_size(attrtype_id, 'H5T_VARIABLE');
H5T.set_size(attrtype_id, numel(type));
H5T.set_strpad(attrtype_id, 'H5T_STR_NULLTERM');

attrspace_id = H5S.create('H5S_SCALAR');
attr_id = H5A.create(dset_id, 'MATLAB_class', ...
                     attrtype_id, attrspace_id, 'H5P_DEFAULT');
H5A.write(attr_id, 'H5ML_DEFAULT', type);

H5A.close(attr_id);
H5S.close(attrspace_id);
H5T.close(attrtype_id);

% For boolean data, add an additional attribute, 'MATLAB_int_decode' to
% tell matlab to decode this data as a boolean.
if strcmpi(type, 'logical')
    attrtype_id = H5T.copy('H5T_STD_U8LE');
    attrspace_id = H5S.create('H5S_SCALAR');
    attr_id = H5A.create(dset_id, 'MATLAB_int_decode', ...
                         attrtype_id, attrspace_id, 'H5P_DEFAULT');
    H5A.write(attr_id, 'H5ML_DEFAULT', uint8(1));

    H5A.close(attr_id);
    H5S.close(attrspace_id);
    H5T.close(attrtype_id);
end

%% Close all the file handles
H5D.close(dset_id);
H5P.close(dcpl_id);
H5S.close(space_id);
H5T.close(type_id);
H5F.close(fid);
end

function [type_id, type_size_bytes, type_fill] = matlab_to_hdf5_type(matlab_type, iscomplex)
% Generates an HDF5 type based on a matlab type
    if strcmpi(matlab_type, 'logical')
        % Special datatype to handle booleans
        h5typename = 'H5T_NATIVE_UCHAR';
        type_id = H5T.copy(h5typename);
        % Only 1-bit precision
        H5T.set_precision(type_id, 1);
        type_size_bytes = H5T.get_size(type_id);
        type_fill = cast(0, 'uint8');
    else
        switch lower(matlab_type)
          case 'double'
            h5typename = 'H5T_NATIVE_DOUBLE';
          case 'single'
            h5typename = 'H5T_NATIVE_FLOAT';
          case 'int64'
            h5typename = 'H5T_NATIVE_LLONG';
          case 'uint64'
            h5typename = 'H5T_NATIVE_ULLONG';
          case 'int32'
            h5typename = 'H5T_NATIVE_INT';
          case 'uint32'
            h5typename = 'H5T_NATIVE_UINT';
          case 'int16'
            h5typename = 'H5T_NATIVE_SHORT';
          case 'uint16'
            h5typename = 'H5T_NATIVE_USHORT';
          case 'int8'
            h5typename = 'H5T_NATIVE_SCHAR';
          case 'uint8'
            h5typename = 'H5T_NATIVE_UCHAR';
          case 'char'
            h5typename = 'H5T_C_S1';
        end

        if iscomplex
            subtype_id = H5T.copy(h5typename);
            subtype_size = H5T.get_size(subtype_id);
            H5T.close(subtype_id);

            type_size_bytes = 2*subtype_size;

            type_id = H5T.create('H5T_COMPOUND', 2*subtype_size);
            H5T.insert(type_id, 'real', 0,            h5typename);
            H5T.insert(type_id, 'imag', subtype_size, h5typename);

            type_fill.real = cast(0, matlab_type);
            type_fill.imag = cast(0, matlab_type);
        else
            type_id = H5T.copy(h5typename);
            type_size_bytes = H5T.get_size(type_id);
            type_fill = cast(0, matlab_type);
        end
    end
end

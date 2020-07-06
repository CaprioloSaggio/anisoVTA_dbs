function [stiff, diinsy, cols, sysmat] = simbio_stiff_matrix(cond, vol)

% based on SB_CALC_STIFF found in the FieldTrip package
%
% 
% INPUT
% cond: conductance tensor represented as a Nx6 array of conductivity
%       values
% vol: volume mesh structure containing nodes in the field called 'pos' and
%      the elements connectivity in a field called 'tet'. It must have N
%      elements
% 
% 
% OUTPUT
% stiff: sparse matrix containing in the element (i,j) the conductivity
%        between nodes i and j
% diinsy:
% cols:
% sysmat:

% TODO: check if I need all those outputs

if(~(size(vol.pos,2)==3))
    if(size(vol.pos,1)==3)
        node = vol.pos';
        warning('Dimension of the nodes of the volume in input should be #nodes x 3!')
    else
        error('The volume structure in input has node field with wrong dimensions!')
    end
else
    node = vol.pos;
end
npnt = size(node,1);
npnt = int32(npnt);

if isfield(vol,'tet')
    if size(vol.tet,1) == 4
        mele = size(vol.tet,1);
        elem = vol.tet;
    elseif size(vol.tet,2) == 4
        mele = size(vol.tet,2);
        elem = vol.tet';
    else
        error('The volume structure in input has tet (tetrahedrals) field with wrong dimensions!')
    end
    elem = [elem; zeros(4,size(elem,2))];
elseif isfield(vol,'hex')
    if size(vol.hex,1) == 8
        mele = size(vol.hex,1);
        elem = vol.hex;
    elseif size(vol.hex,2) == 8
        mele = size(vol.hex,2);
        elem = vol.hex';
    else
        error('vol.hex has wrong dimensions!')
    end
else
    error('Could not find connectivity information!')
end

if min(min(elem(1:mele,:))) == 0
    elem = elem + 1;
    warning('Numbering of nodes in vol.tet/vol.hex must start at 1 (Fortran numbering)!')
elseif min(min(elem(1:mele,:))) < 0
    error('No negative indices for conectivity information allowed!')
end

if size(vol.tet, 1) == size(cond, 1)
    cond = double(cond);
else
    error("The number of conductivity values doesn't match the number of elements in the mesh")
end

mele = int32(mele);
elem = int32(elem);

% check whether the nodes have right orientation

if isfield(vol,'tet')
    if ~sb_test_ori(node,elem(1:4,:)')
        error('Elements have wrong orientation, consider exchanging node 3 and 4');
    end
elseif isfield(vol,'hex')
    if ~sb_test_ori(node,elem')
        error('Elements have wrong orientation or are degenerated');
    end
end

try
    [diinsy,cols,sysmat] = calc_stiff_matrix_val(node,elem,cond,mele);
catch err
    if ispc && strcmp(err.identifier,'MATLAB:invalidMEXFile')
        error('Error executing mex-file. Microsoft Visual C++ 2008 Redistributables and Intel Visual Fortran Redistributables are required.')
    else
        rethrow(err)
    end
end
npnt = double(npnt);
diinsy = double(diinsy);
cols = double(cols);
rows = sb_sparse_to_mat(diinsy);
stiff = sparse(rows,cols,sysmat,npnt,npnt,length(sysmat));
end

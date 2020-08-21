function [b_val, b_vect] = extract_b(folder, save_option)

% TODO:
% - check if the case b_value==0 is correctly handled (at the moment a [0;0;0] 
%   array is added to b_vect)

% INPUT: 
% folder: full path of the folder containing all the dicom files of the 
% acquired volume
% save_option: option that, if true, saves b_val and b_vect vairables in 2 
% homonimous .txt files

% !!! --> the script reads all the files in the directory and takes into
% account only the dicom ones

% OUTPUT
% b_val: 1xN array containing the b-value used during the N acquisitions
% b_vect: 3xN array containing the gradient orientations used during the N 
% acquisitions
% OPTIONAL: you can save b_val and b_vect in 2 .txt files with 'space' 
%           separator


if nargin == 1
    save_option = 1;
elseif nargin == 2 && not(save_option ~= 0 || save_option ~= 1)
    error("Invalid value given as an input to save: legal assignments are 0 or 1")
elseif nargin > 2
    error("Too many input parameters: legal number of variables in input is 2 --> directory name and save option")
else
    if not(exist(folder) == 7)  % exists returns 7 if the argument is an existing folder
        error("First input is not found or it is not a valid directory, please check it")
    end
end

% initialize variables
b_val = [];
b_vect = [];
files_struct = dir(folder);

% scan the folder looking for the dicom files, from which extract b_val and
% b_vect
for f = 3:length(files_struct)
    try
        file = files_struct(f).name;
        path = fullfile(folder, file);
        header = dicominfo(path);
        b_ = header.Private_0019_100c;
        b_val = [b_val b_];
        if b_ == 0
            b_vect = [b_vect [0; 0; 0]];
        else
            b_vect = [b_vect header.Private_0019_100e];
            
        end
    catch
        warning("No DICOM header could be found for file: %s", file)
    end
end

% save into txt file
if save_option
    writematrix(b_val, fullfile(folder, 'bvals'), 'Delimiter', ' ');  % if the file already exists, it is overwritten 
    writematrix(b_vect, fullfile(folder, 'bvecs'), 'Delimiter', ' ');
end

end
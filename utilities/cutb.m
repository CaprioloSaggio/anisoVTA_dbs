function [bvals_short, bvecs_short, mrtrix_scheme] = cutb(varargin)
%
% INPUT 
% bvals: name of the file containing bvals according to fsl representation
%        of the gradient scheme
% bvecs: name of the file containing bvecs according to fsl representation
%        of the gradient scheme
% bvalsName: (optional) name for the text file to print in output, default 
%            is bvals_short.txt
% bvecsName: (optional) name for the text file to print in output, default 
%            is bvecs_short.txt
% grad_mrtrix: (optional) name for the text file to print in output, default 
%              is grad_mrtrix.txt
% 
% OUTPUT
% bvals_short: array containing the truncated bvals
% bvecs_short: array containing the truncated bvecs
% the function also saves two text files named bvals_short.txt and
% bvecs_short.txt
% 

% handle inputs to the function
if length(varargin)<3 || length(varargin)>6
    error("Wrong number of input. Type 'help cutb' to get information on how to use it")
end

bvals = varargin{1};
bvecs = varargin{2};
n = varargin{3};

if length(varargin)>5
    grad_mrtrix = varargin{6};
    bvecsName = varargin{5};
    bvalsName = varargin{4};
elseif length(varargin)>4    
    grad_mrtrix = 'grad_mrtrix.txt';
    bvecsName = varargin{5};
    bvalsName = varargin{4};
elseif length(varargin)>3
    grad_mrtrix = 'grad_mrtrix.txt';
    bvecsName = 'bvecs_short.txt';
    bvalsName = varargin{4};
else
    grad_mrtrix = 'grad_mrtrix.txt';
    bvecsName = 'bvecs_short.txt';
    bvalsName = 'bvals_short.txt';
end

% cut bvals at size n
bvals = importdata(bvals);
bvals_short = bvals(:,1:n); 

% save bvals into text file
fid = fopen( bvalsName, 'wt' );
for i = 1:size(bvals_short, 1)
    input = num2str(bvals_short(i,:));
    fprintf( fid, '%s\n', input);
end
fclose(fid);

% cut bvecs at size n
bvecs = importdata(bvecs);
bvecs_short = bvecs(:,1:n);

% save bvecs into text file
fid = fopen( bvecsName, 'wt' );
for i = 1:size(bvecs_short, 1)
    input = num2str(bvecs_short(i,:));
    fprintf( fid, '%s\n', input);
end
fclose(fid);


% along with the composition of the bvals and bvecs files the script also
% provides an mrtrix format gradient scheme

% build the gradient scheme from the fsl representation
mrtrix_scheme = [bvecs_short', bvals_short'];

% save mrtrix gradient scheme into text file
fid_mrtrix = fopen(grad_mrtrix, 'wt');
for i = 1:size(mrtrix_scheme, 1)
    input = num2str(mrtrix_scheme(i,:));
    fprintf( fid, '%s\n', input);
end
fclose(fid_mrtrix);

end
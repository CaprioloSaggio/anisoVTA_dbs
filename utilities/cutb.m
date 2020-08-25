function [bvals_short, bvecs_short] = cutb(bvals, bvecs, n)
%
% INPUT 
% bvals: name of the file containing bvals
% bvecs: name of the file containing bvecs
% 
% OUTPUT
% bvals_short: array containing the truncated bvals
% bvecs_short: array containing the truncated bvecs
% the function also saves two text files named bvals_short.txt and
% bvecs_short.txt
% 


% cut bvals at size n
bvals = importdata(bvals);
bvals_short = bvals(:,1:n);

% save bvals into text file
fid = fopen( 'bvals_short.txt', 'wt' );
for i = 1:size(bvals_short, 1)
    input = num2str(bvals_short(i,:));
    fprintf( fid, '%s\n', input);
end
fclose(fid);

% cut bvecs at size n
bvecs = importdata(bvecs);
bvecs_short = bvecs(:,1:n);

% save bvecs into text file
fid = fopen( 'bvecs_short.txt', 'wt' );
for i = 1:size(bvecs_short, 1)
    input = num2str(bvecs_short(i,:));
    fprintf( fid, '%s\n', input);
end
fclose(fid);

end
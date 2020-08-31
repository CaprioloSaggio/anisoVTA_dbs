function isotropic_conductivity_stn = teststn(vol, cond, side, options)

% Here's how to check the conductivity values of the STN, expected to be
% around 0.33, as it is grey matter

try
    load stn
catch
    % find STN data
    load([options.root, options.patientname, filesep, 'atlases', filesep, options.atlasset, filesep, 'atlas_index.mat']);
    stnl = atlases.fv{1, 1}.vertices;
    stnr = atlases.fv{1, 2}.vertices;
    save('stn.mat', 'stnl', 'stnr');
    % load stn points
    load stn
end

% select side
if side==1
    stn = stnl;
elseif side==2
    stn=stnr;
end

% find conductivity in the STN
test = unique(knnsearch(vol.ctr, stn));
test_c = cond(test, :);
test_cm = mean(test_c);

isotropic_conductivity_stn = isoCond(test_cm);

end  % function


function I = isoCond(condArray)
% 
% Returns the equivalent isotropic conductivity I of the array condArray 
% describing the ellipsoid for the anisotropic conductivity
% 

S = makeSymmetric(condArray);
[~,D] = eig(S);
D = diag(D);
I = (prod(D))^(1/3);

end  % subfunction


function S33 = makeSymmetric(array6)

S33 = zeros(3);
S33(1,1) = array6(1);
S33(1,2) = array6(2);
S33(2,1) = array6(2);
S33(1,3) = array6(3);
S33(3,1) = array6(3);
S33(2,2) = array6(4);
S33(2,3) = array6(5);
S33(3,2) = array6(5);
S33(3,3) = array6(6);

end  % subfunction



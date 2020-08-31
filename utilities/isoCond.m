function I = isoCond(condArray)
% 
% Returns the equivalent isotropic conductivity I of the array condArray 
% describing the ellipsoid for the anisotropic conductivity
% 

S = makeSymmetric(condArray);
[~,D] = eig(S);
D = diag(D);
I = (prod(D))^(1/3);

end

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

end
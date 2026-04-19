function grad_of = gradOf(u, nodes, elements)
% Compute gradient of compliance w.r.t. cross-sectional areas x
% grad_of(j) = - u_el' * (dK/dx_j) * u_el
% We build element stiffness matrices for unit area so that
% dK/dx_j = kelG_unit(:,:,j)
n_el = size(elements,1);

[~, ~, ~, kelG_unit] = buildStiffnessMatrix(nodes, elements, ones(n_el,1));

grad_of = zeros(n_el,1);
for ii = 1:n_el
    n1 = elements(ii,1);
    n2 = elements(ii,2);
    dofel = [ 2*n1-1, 2*n1, 2*n2-1, 2*n2];
    uelG = u(dofel);
    grad_of(ii) = - uelG' * kelG_unit(:,:,ii) * uelG;
end
end

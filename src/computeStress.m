function sigma = computeStress(u, x, nodes, elements)
% INPUTS:
% u, vector of nodes displacements (dim: 2*n_nodes x 1)
% x, vector of cross sectional areas
% nodes, matrix of node coordinates
% elements, matrix of nodes constituting elements

% OUTPUTS:
% sigma, vector of element stresses

	[~,nelem,R,kelG] = buildStiffnessMatrix(nodes,elements,x);
	sigma = zeros(nelem,1);
	for ii = 1:nelem
		dofel = [	2*elements(ii,1) - 1, 2*elements(ii,1), ...
							2*elements(ii,2) - 1, 2*elements(ii,2)]';
		uelG = u(dofel);
		FLocal = R(:,:,ii)*kelG(:,:,ii)*uelG;
		sigma(ii,1) = FLocal(3) / x(ii);
	end
end
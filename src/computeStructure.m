function [u, F, dof] = computeStructure(constr, ext_load, K)
% INPUTS:
% constr, matrix representing the constraints
% ext_load, matrix representing the external loads
% K, global stiffness matrix of the truss

% OUTPUTS:
% u, vector of nodes displacements (dim: 2*n_nodes x 1)
% F, vector of nodes forces (dim: 2*n_nodes x 1)
% dof, logical vector of nodes dof (dim: 2*n_nodes x 1). True if free

	n_constr = size(constr,1);
	n_load = size(ext_load,1);
	sz = size(K,1);
	dof = true(sz,1);

	for ii = 1:n_constr
		dof(2*constr(ii,1)+constr(ii,2)-2) = false;
	end
	Kfree = K(dof, dof);
	F = zeros(sz,1);
	for ii = 1:n_load
		F(2*ext_load(ii,1)+ext_load(ii,2)-2) = ext_load(ii,3);
	end
	u = zeros(size(K,1),1);
	u(dof) = Kfree\F(dof);
end
function [K, n_el, R, kelG] = buildStiffnessMatrix(nodes, elements, varargin)
% INPUTS:
% nodes, matrix of node coordinates
% elements, matrix of nodes constituting elements
% x (optional), vector of cross sectional areas

% OUTPUTS:
% K, global stiffness matrix of the truss
% n_el, number of elements
% R, matrix 4x4xn_el, collection of rotation matrices
% KelG, matrix 4x4xn_el, collection of element stiffness matrices in global CSYS

	n_el = size(elements,1);
	if nargin > 3
		error('Number of inputs cannot exceed 3.')
	end
	if nargin < 2
		error('Number of inputs must be at least 2.')
	end
	if nargin == 2
		x = 100*ones(n_el,1);
	else
		x = varargin{1};
		if setdiff(length(x),n_el)
			error('Vector x must contains %d elements.',n_el)
		end
	end

	lengths = zeros(n_el,1);
	kelL = zeros(4,4,n_el);
	for ii = 1:n_el
		x1 = nodes(elements(ii,1),1); 
		x2 = nodes(elements(ii,2),1);
		y1 = nodes(elements(ii,1),2); 
		y2 = nodes(elements(ii,2),2);
		lengths(ii) = sqrt((y2-y1)^2 + (x2-x1)^2);
		stiff = x(ii)/lengths(ii)*210E03;
		kelL(1,1,ii) = stiff;
		kelL(1,3,ii) = -stiff;
		kelL(3,1,ii) = -stiff;
		kelL(3,3,ii) = stiff;
	end

	alphas = zeros(n_el,1);
	R = zeros(4,4,n_el);
	kelG = kelL;
	for ii = 1:n_el
		x1 = nodes(elements(ii,1),1); 
		x2 = nodes(elements(ii,2),1);
		y1 = nodes(elements(ii,1),2); 
		y2 = nodes(elements(ii,2),2);
		alphas(ii) = atan2((y2-y1),(x2-x1));
		c = cos(alphas(ii)); 
		s = sin(alphas(ii));
		R(1,1,ii) = c;	R(1,2,ii) = s;
		R(2,1,ii) = -s;	R(2,2,ii) = c;
		R(3,3,ii) = c;	R(3,4,ii) = s;
		R(4,3,ii) = -s;	R(4,4,ii) = c;
		kelG(:,:,ii) = R(:,:,ii)'*kelL(:,:,ii)*R(:,:,ii);
	end

	K = zeros(2*size(nodes,1));
	for ii = 1:n_el
		dofel = [	2*elements(ii,1) - 1, 2*elements(ii,1), ...
							2*elements(ii,2) - 1, 2*elements(ii,2)];
		idx = ismember(1:2*size(nodes,1), dofel);
		Ktemp = zeros(2*size(nodes, 1));
		Ktemp(idx,idx) = kelG(:,:,ii);
		K = K + Ktemp;
	end
end
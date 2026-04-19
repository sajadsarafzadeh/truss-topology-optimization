function [V, lengths] = computeVolume(x, nodes, elements)
% INPUTS:
% nodes, matrix of node coordinates
% elements, matrix of nodes constituting elements
% x, vector of cross sectional areas

% OUTPUTS:
% V, truss volume
% lengths, lengths of truss elements

	n_el = size(elements,1);
	lengths = zeros(n_el,1);
	for ii = 1:n_el
		x1 = nodes(elements(ii,1),1); 
		x2 = nodes(elements(ii,2),1);
		y1 = nodes(elements(ii,1),2); 
		y2 = nodes(elements(ii,2),2);
		lengths(ii) = sqrt((y2-y1)^2 + (x2-x1)^2);
	end
	
	V = lengths'*x;
end
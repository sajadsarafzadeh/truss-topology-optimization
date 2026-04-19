function grad_constr = gradConstr(nodes, elements)
% INPUTS:
% nodes, matrix of node coordinates
% elements, matrix of nodes constituting elements

% OUTPUTS:
% grad_constr, vector of the gradient of the constraint (volume) computed
% with respect to x (cross sectional areas)

% g1(x) = V(x)-V_0
% V(x) = l^T*x       x: section area of each element
% dg1(x)/dx_j = l_j
% it means that the derivative of each function with respect to the design
% variable (x) is equal to the length


    
    n_el = size(elements,1);
    grad_constr = zeros(n_el,1); % This vector will store the lengths

    for ii = 1:n_el
        % Get coordinates of the two nodes for this element
        x1 = nodes(elements(ii,1),1); 
        x2 = nodes(elements(ii,2),1);
        y1 = nodes(elements(ii,1),2); 
        y2 = nodes(elements(ii,2),2);
        
        % Calculate length (which is the gradient) 
        grad_constr(ii) = sqrt((y2-y1)^2 + (x2-x1)^2);
    end

end
   

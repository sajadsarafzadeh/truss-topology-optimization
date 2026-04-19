function [g, x_star] = gradPhi(lambda, x_k, problem)
% INPUTS:
% lambda, value of Lagrange multiplier
% x_k, point with respect to which the problem is linearised
% problem, struct with problem data

% OUTPUTS:
% x_star, value of the optimal x, computed as a function the value of
% lambda supplied
% g, value of the constraint evaluated at x_star (must be = 0 when 
% lambda = lambda_star!!)

	% -------------------------
	% Impose consistency on lambda i.e. ensure that it is >= 0
	%
	lambda = max(1e-10,lambda);

	% -------------------------
	% Collect gradients from problem data
	%
	grad_of = problem.grad_of;
	grad_constr = problem.grad_constr;

	% -------------------------
	% Compute x_star(lambda) depending on the sign of the gradients
	%
	x_star = zeros(problem.n_el,1);
	for ii = 1:problem.n_el
		if grad_of(ii) > 0 && grad_constr(ii) > 0
			x_star(ii) = problem.x_lb(ii);
				;
		elseif grad_of(ii) > 0 && grad_constr(ii) < 0
			x_star(ii) = x_k(ii)*sqrt(-lambda*grad_constr(ii)/grad_of(ii));
				;
		elseif grad_of(ii) < 0 && grad_constr(ii) > 0
			x_star(ii) = x_k(ii)*sqrt(-grad_of(ii)/(lambda*grad_constr(ii)));
				;
		elseif grad_of(ii) < 0 && grad_constr(ii) < 0
			x_star(ii) = problem.x_ub(ii);
				;
		end
	end
	% -------------------------
	% Impose box constraints on x_star
	%
	x_star = max(problem.x_lb, min(problem.x_ub, x_star));

	% -------------------------
	% Evaluate the constraint violation
	%
	g = problem.computeVolume(x_star, problem.nodes, problem.elements) - problem.V;
end
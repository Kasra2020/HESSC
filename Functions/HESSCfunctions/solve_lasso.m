function C = solve_lasso(X, Y, lambda)
%  This function solves the following lasso problem
%  min_{c_j}||c_j||_1 + lambda/2 ||y_j - X c_j||_2^2
%  where y_j is the j-th column of Y. The output is C = [c_1, c_2, ...]
%  We use the SPAMS package http://spams-devel.gforge.inria.fr/ for solving
%  the lasso problem. You can also replace it with your favorate solver.
param.lambda = 1 / lambda;
C = mexLasso(Y, X, param);
%if you want rowmedian the same height as S, replace [] by size(S, 1)
% lambda=1/lambda;
% C=lasso(Y,X,'Lambda',lambda);
% C=C';
end

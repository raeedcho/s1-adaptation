% 
% Compute principal angles between plabes A and B
%
%   [angles, varargout] = principal_angles(A,B)
%
% Inputs:
%   A           : n-by-p matrix that defines a p-dimensional plane in an
%               n-dimensional Euclidean space
%   B           : n-by-q matrix that defines a q-dimensional plane in an
%               n-dimensional Euclidean space
%
% Outputs (opt):
%   angles      : principal angles
%   (C)         : singular values of Qa'*Qb
%   (U)         : matrix S when doing [U C V] = svd(Qa'*Qb)
%   (V)         : matrix S when doing [U C V] = svd(Qa'*Qb)
%
%
% Syntax:
%   angles      = principal_angles(A,B)
%   [angles, C]  = principal_angles(A,B)
%   [angles, U, C, V]  = principal_angles(A,B)
%
%
% Note: inspired by this post: http://goo.gl/eMclEw
%


function [angles, varargout] = principal_angles(A,B)

[Qa,~]          = qr(A,0);
[Qb,~]          = qr(B,0);
[U,C,V]         = svd(Qa'*Qb,0);
% make sure none of the elements in C is > 1
angles          = acos( min( ones(length(diag(C)),1), diag(C) ) );

if nargout == 2
    varargout{1} = C;
end
if nargout == 4
    varargout{1} = U;
    varargout{2} = C;
    varargout{3} = V;
end
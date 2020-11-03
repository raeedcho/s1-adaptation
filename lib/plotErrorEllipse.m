function plotErrorEllipse(mu, Sigma, p, varargin)
% PLOTERRORELLIPSE - function to plot confidence ellipse for bivariate gaussian
%   Based on: https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
%   Inputs:
%       mu - mean of distribution
%       Sigma - covariance matrix of distribution
%       p - confidence level
%       varargin - list of additional plotting parameters
%
%   Explanation:
%       For X ~ N(mu,Sigma), goal is to find ellipse that contains fraction p of data.
%       For uncorrelated X, ellipse is axis-aligned and at boundary,
%       sum((x_i/sigma_i)^2) = s
%       for some value s, which follows a chi-square distribution with degrees of
%       freedom equal to number of dimensions in X. Need to find s_p such that
%       P(s<s_p) = p, which will give ellipse boundary. This requires inversion of
%       the chi-square CDF.
%       For general case of correlated X, we use eigenvectors of Sigma as the
%       axes of the ellipse, and the eigenvalues as the variances along the axes.

    s_p = chi2inv(p,size(mu,2));

    [V, D] = eig(Sigma * s_p);

    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];

    plot(a(1, :) + mu(1), a(2, :) + mu(2),varargin{:});
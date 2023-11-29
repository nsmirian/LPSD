function y = hlc_calc_gaussian(x, para)
% Calculate a Gaussian distribution for a given set of parameters.
%
% y = hlc_calc_gaussian(x, parameters)
%   Returns a Gaussian distribution for the given set of x values and parameters.
%   The Gaussian is calculated as follows:
%   y = offset + amplitude * exp(-(x-mu)^2 / 2 / sigma^2)
%
% PARAMETERS
%   parameters = [ offset, amplitude, mu, sigma ]
%
% RETURN VALUE
%   An array containing y(x) for every entry in x.
%
% See also hlc_fit_gaussian.
%
% -- Lars Froehlich, 2007

    offset = para(1);
    amp = para(2);
    mu = para(3);
    sigma = para(4);

    y = offset + amp .* exp(-(x-mu).^2/2/sigma.^2);
end

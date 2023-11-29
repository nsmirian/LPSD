function fit_data = hlc_fit_gaussian(x, y)
% Fit a Gaussian distribution to the given data set.
%
% fit_data = hlc_fit_gaussian(x, y)
%   Fits a Gaussian distribution to the given data set. The Gaussian
%   is calculated as follows:
%   y = baseline + height * exp( -(x - mu)^2 / 2 / sigma^2 )
%
% PARAMETERS
%   x - x values
%   y - y values
%
% RETURN VALUE
%   fit_data - Array containing the fit results
%              [baseline, height, mu, sigma]
%
% See also hlc_calc_gaussian.
%
% -- Lars Froehlich, 2007

    fit_data = [];
    if (isempty(x) || isempty(y))
        return;
    end

    % Ignore array orientation
	x = x(:);
	y = y(:);

    [fit_data_pos, sumsq_error_pos] = fit_positive_gaussian(x, y);
    [fit_data_neg, sumsq_error_neg] = fit_positive_gaussian(x, -y);
    
    if (sumsq_error_neg < sumsq_error_pos)
        fit_data = fit_data_neg .* [-1, -1, 1, 1];
    else
        fit_data = fit_data_pos;
    end

    fit_data = fminsearch(@err_gaussian, fit_data, ...
        optimset('display', 'off', 'maxiter', 1500, 'maxfunevals', 1500), [x, y]);
end


function [fit_data, sumsq_error] = fit_positive_gaussian(x, y)
	% determine starting parameters for the fit
	baseline = min(y);
	height = max(y) - baseline;
	pb = y - baseline;
	pbx = pb .* x;
	pbxsq = pbx .* x;
    sum_pb = sum(pb);
    if (sum_pb == 0)
        sum_pb = 1e-3;
    end
	mu = sum(pbx) / sum_pb;
	sigma = sqrt( sum(pbxsq) / sum_pb - mu^2 );

    [fit_data, sumsq_error] = fminsearch(@err_gaussian, [baseline, height, mu, sigma], ...
        optimset('display', 'off', 'maxiter', 1500, 'maxfunevals', 1500), [x, y]);
end



%
% err = err_gaussian(parameters, data)
%
% Returns the total quadratic deviation of a Gaussian function
% (specified by the parameters) from the given data.
%
% parameters = [offset, amplitude, mu, sigma]
% data(:, 1) - x values
% data(:, 2) - y values
%
function err = err_gaussian(parameters, data)
	err = (data(:,2) - hlc_calc_gaussian(data(:,1), parameters)).^2;
	err = sum(err);
end

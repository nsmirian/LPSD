function [com, rms, fwhm, a, p_out] = get_profile_stats(p, c)

    % p: profile
    % c: calibration factor
    
    n       = size(p,2); % number of bins
    a       = c * ( (1:n) - n/2 ); % axis
        
    com     = zeros([size(p,1), 1]);
    rms     = zeros([size(p,1), 1]);
    fwhm    = zeros([size(p,1), 1]);
    p_out   = p;
    for ii = 1:size(p,1)
        
        p_sum           = sum(p(ii,:));
        
        if (p_sum ~= 0)
    
            p_n             = p(ii,:)/p_sum ; % normalized profile
            p_n(p_n<0)      = 0;

            com(ii)         = a*p_n(:);
            rms(ii)         = sqrt( (a-com(ii)).^2*p_n(:) );

            p_fhwm          = p_n;
            p_fhwm(p_fhwm<max(p_fhwm)/2) = 0;
            tmp             = find(p_fhwm ~= 0);
            fwhm(ii)        = abs(diff(a(tmp([1,end]))));

            p_out(ii,:)     = 1/abs(c) * p_n(:);
        
        else
            
            com(ii)         = NaN;
            rms(ii)         = NaN;
            fwhm(ii)        = NaN;
            
        end

    end

end
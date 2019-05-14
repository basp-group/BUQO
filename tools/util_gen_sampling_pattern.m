function [u, v] = util_gen_sampling_pattern(pattern, param)

if ~isfield(param, 'sigma')
    % variance of the gaussian for generating measurements
    param.sigma = pi/3;
end

if strcmp(pattern, 'gaussian')
    sigma_m = param.sigma;
    Nm = round(param.p * param.N);
    
    u = sigma_m * randn(Nm, 1);
    v = sigma_m * randn(Nm, 1);

    % discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    while length(sf) < Nm
        Nmextra = 2 * (Nm - length(sf));
        u = [u; sigma_m * randn(Nmextra, 1)];
        v = [v; sigma_m * randn(Nmextra, 1)];
        
        
        % discard points outside (-pi,pi)x(-pi,pi)
        sfu = find((u<pi) & (u>-pi));
        sfv = find((v<pi) & (v>-pi));
        sf = intersect(sfu, sfv);
    end
    
    vw = v(sf(1:Nm));
    uw = u(sf(1:Nm));
    
    clear u;
    clear v;
    
    u = uw;
    v = vw;
    
    Nm = length(uw);
end
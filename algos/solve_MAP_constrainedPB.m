function result = solve_MAP_constrainedPB(param_data,param_algo)
% Condat-Vu primal-dual algorithm to solve 
% minimize  || Psit x ||_1 
%  s.t. x positive and  || Phi x - y || <= l2bound

%% Initialization


lips = param_data.normPhi + param_algo.normPsi ;

if ~isfield(param_algo,'tau'), tau = 0.99/sqrt(lips) ;
else tau = param_algo.tau ; end
if ~isfield(param_algo,'gamma'), gamma = 0.99 / (lips*tau) ; 
else gamma = param_algo.gamma ; end


if gamma * tau * lips >= 1
    disp('Error convergence parameters: default parameters chosen')
    tau = 0.99/sqrt(lips) ;
    gamma = tau ;
end

x = param_algo.x0 ;
Phix = param_data.Phi(x) ;
v = 0*param_algo.Psit(x) ;
u = 0*Phix ;

lamb_conv = 1 ;


%% Iterations

Crit(1) = sum(abs(param_algo.Psit(x))) ;
condl2(1) = sqrt(sum(abs(Phix(:)-param_data.y(:)).^2)) ;

disp(' ')
disp('*********************************')
disp('Initialization')
disp(['l1 norm: ',num2str(Crit(end))])
disp(['l2 norm: ',num2str(condl2(end)),' vs. bound: ',num2str(param_data.l2bound)])
disp(['tau: ',num2str(tau)])
disp(['gamma: ',num2str(gamma)])



for it = 1:param_algo.NbIt
Phixold = Phix ;
xold = x ;
uold = u ;
vold = v ;

start_it = tic;

xtmp = x - tau * ( real(param_data.Phit( u )) + param_algo.Psi(v) ) ;
x = proj_min_max(xtmp, 0, Inf) ;
x = lamb_conv * x + (1-lamb_conv) * xold ;

Phix = param_data.Phi(x) ;
utmp = u+ gamma* (2*Phix - Phixold) ; 
p1 = proj_l2ball(utmp/gamma, param_data.l2bound, param_data.y) ;
u = utmp - gamma * p1 ;
u = lamb_conv * u + (1-lamb_conv) * uold ;

vtmp = v+ gamma* param_algo.Psit(2*x-xold) ;
p2 = prox_l1(vtmp/gamma, param_algo.lambda/gamma) ;
v = vtmp - gamma * p2 ;
v = lamb_conv * v + (1-lamb_conv) * vold ;
  
Crit(it+1) = sum(abs(param_algo.Psit(x))) ;
condl2(it+1) = sqrt(sum(abs(Phix(:)-param_data.y(:)).^2)) ;
norm_it(it) = norm(x(:) - xold(:))/norm(x(:)) ;
Time_it(it) = toc(start_it) ;

if mod(it, param_algo.display)==0
disp('*********************************')
disp(['it ',num2str(it)])
disp(['l1 norm         = ',num2str(Crit(end))])
disp(['l2 norm         = ',num2str(condl2(end))])
disp(['   vs. bound    = ',num2str(param_data.l2bound)])
disp(['norm it         = ',num2str(norm_it(end))])
disp(['   vs stop norm = ',num2str(param_algo.stop_norm)])
disp(['crit it         = ',num2str(abs(Crit(end)-Crit(end-1))/Crit(end))])
disp(['   vs stop crit = ',num2str(param_algo.stop_crit)])
disp(['cond it         = ',num2str(abs(condl2(it+1) - condl2(it))/condl2(it))])
disp(['   vs stop cond = ',num2str(param_algo.stop_cond)])

figure(101)
imagesc(log10(x)), axis image; colorbar; colormap jet ; caxis([-3.5 0])
pause(0.1)
end

if norm_it(end)<param_algo.stop_norm ...    % variation iterates
        &&  abs(Crit(end)-Crit(end-1))/Crit(end) < param_algo.stop_crit ...   %variation criterion
        &&  condl2(end) <= param_data.l2bound ... % l2 ball 
        &&  abs(condl2(it+1) - condl2(it)) < param_algo.stop_cond * condl2(it)
    break;
end

end

result.x = x ;
result.crit = Crit ;
result.condl2 = condl2 ;
result.norm_it = norm_it ;
result.Time_it = Time_it ;
disp('*********************************')
disp(['it ',num2str(it)])
disp(['Rec. time       = ',num2str(sum(Time_it))])
disp(['l1 norm         = ',num2str(Crit(end))])
disp(['l2 norm         = ',num2str(condl2(end))])
disp(['   vs. bound    = ',num2str(param_data.l2bound)])
disp(['norm it         = ',num2str(norm_it(end))])
disp(['   vs stop norm = ',num2str(param_algo.stop_norm)])
disp(['crit it         = ',num2str(abs(Crit(end)-Crit(end-1))/Crit(end))])
disp(['   vs stop crit = ',num2str(param_algo.stop_crit)])
disp(['cond it         = ',num2str(abs(condl2(it+1) - condl2(it))/condl2(it))])
disp(['   vs stop cond = ',num2str(param_algo.stop_cond)])
disp('*********************************')

end

%%
function p = proj_min_max(x,xmin,xmax) 
p = max(min(x,xmax),xmin) ;
end

function p = prox_l1(z, T) 
p = sign(z).*max(abs(z)-T, 0);
end

function p = proj_l2ball(x, eps, y)
% projection of x onto the l2 ball centered in y with radius eps
p = x-y ;
p = p* min(eps/norm(p(:)),1) ;
p = p+y ;

end

function crit = compute_crit(x, sigma2, Phi, y, Psit, lambda)
res = Phi(x) - y ;
fid = (1/(2*sigma2)) * sum( abs(res(:)).^2 ) ;
Px = Psit(x);
reg = lambda * sum( abs( Px(:) ) ) ;
crit = fid + reg ;
end

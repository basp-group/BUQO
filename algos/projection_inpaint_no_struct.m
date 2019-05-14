function [x,v] = projection_inpaint_no_struct(xhpd_in, xC_out, v, param)

if ~isfield(param,'display'), param.display = 100 ; end
if ~isfield(param,'final_display'), param.final_display = 1 ; end
%% Parameters

beta = param.normL + 1 ;

gamma_ = 2 ;
gamma = gamma_ / (2*param.normL) ;
tau = 1.9 / (beta + gamma_ ) ;

if tau * (beta/2 + gamma * param.normL) >=1
    disp('Error convergence parameters for HPD constraint')
end

%% Initialization 


u=param.x0;

x = 0*param.Mask(:) ;
x(param.Mask==0) = u ;
x(param.Mask>0) = param.L*u ;
x = reshape(x, size(param.Mask,1), size(param.Mask,2)) ;


res = sum(abs(param.L*u - xhpd_in).^2) + sum(abs(u - xC_out).^2) ;
crit(1) = sum(abs(res).^2) ;
Lx = param.L*u ;
l2_cond(1) = sqrt(sum(abs(Lx - param.l2_mean).^2)) ;
pos_cond(1) = min(u(:)) ;


grad =@(u) param.L'*(param.L*u - xhpd_in) + (u - xC_out) ;

%% Iterations


for it = 1:param.NbIt
xold = x ;
vold = v ;
uold = u ;
    

v_ = v + gamma * Lx ;
v = v_ - gamma * proj_l2ball( gamma^(-1) * v_, param.l2_bound, param.l2_mean ) ;

u_ = u - tau* ( grad(u) +  param.L'*(2*v-vold) ) ;
u = max(u_, 0) ;

x = 0*param.Mask(:) ;
x(param.Mask==0) = u ;
x(param.Mask>0) = param.L*u ;
x = reshape(x, size(param.Mask,1), size(param.Mask,2)) ;

res = sum(abs(param.L*u - xhpd_in).^2) + sum(abs(u - xC_out).^2) ;
crit(it+1) =sum(abs(res).^2) ;
Lx = param.L*u ;
l2_cond(it+1) = sqrt(sum(abs(Lx - param.l2_mean).^2)) ;
pos_cond(it+1) = min(u(:)) ;

% if mod(it, param.display)==0
%     figure(101)
%     subplot 311
%     hold off, plot(crit,'k'), 
%     xlabel('it'), ylabel('crit')
%     subplot 312
%     hold off, semilogy(abs(l2_cond-param.l2_bound),'r'), 
%     xlabel('it'), ylabel('dist to MAP inpaint')
%     subplot 313
%     hold off, plot(pos_cond,'r'), hold on, plot(0*ones(size(pos_cond)),'b--')
%     xlabel('it'), ylabel('min x')
%     pause(1)
% end

if norm(u(:)-uold(:))<1e-5*norm(u(:)) ...         % variation iterates
        && abs(crit(it+1)-crit(it))<1e-5*crit(it+1) ...  % variation critere
        && pos_cond(it+1) >= 0 ...                % positivity constraint
        && l2_cond(it+1) - param.l2_bound < 1e-5         % smoothness constraint
    if param.final_display == 1
    disp(['Projection onto S, it max: ',num2str(it)])
    end
    break;
end

end


%% Final display 

%     figure(101)
%     subplot 311
%     hold off, plot(crit,'k'), 
%     xlabel('it'), ylabel('crit')
%     subplot 312
%     hold off, semilogy(l2_cond,'r'), hold on, semilogy(param.l2_bound*ones(size(l2_cond)),'b--')
%     xlabel('it'), ylabel('dist to inpaint')
%     subplot 313
%     hold off, plot(pos_cond,'r'), hold on, plot(0*ones(size(pos_cond)),'b--')
%     xlabel('it'), ylabel('min x')
% pause(1)

x = 0*param.Mask(:) ;
x(param.Mask==0) = u ;
x(param.Mask>0) = param.L*u ;
x = reshape(x, size(param.Mask,1), size(param.Mask,2)) ;


end




%%
function p = proj_l2ball(x, eps, y)
% projection of x onto the l2 ball centered in y with radius eps
p = x-y ;
p = p* min(eps/norm(p(:)),1) ;
p = p+y ;

end


function [x,z1,z2] = projection_inpaint_constrainedPB(xbar, z1,z2, param)

if ~isfield(param,'display'), param.display = 100 ; end
%% Parameters

beta = param.normLbar ;

gamma_ = 2 ;
gamma = gamma_ / (2*param.normL) ;
tau = 1.9 / (beta + gamma_ ) ;

disp(['gamma = ',num2str(gamma)])
disp(['tau = ',num2str(tau)])


if tau * (beta/2 + gamma * param.normL) >=1
    disp('Error convergence parameters for HPD constraint')
end

%% Initialization 

% % % Ls =@(x) param.L * x(param.Mask==0) - x(param.Mask>0) ;
% % % Lst =@(z)  Lt(z,param.L,param.Mask) ;

x=param.x0;
crit(1) = sum(abs(x(:)-xbar(:)).^2) ;
Psix = param.Psit(x) ;
l1_cond(1) = sum(abs(Psix)) ;
Lsx = Ls(x) ;
smooth_cond(1) = min(abs(Lsx)) ;
pos_cond(1) = min(x(:)) ;

[Ny,Nx] = size(x) ;

%% Iterations


for it = 1:param.NbIt
z1old = z1 ;
z2old = z2 ;
xold = x ;
    
z1_ = z1 + gamma * Psix ;
z1 = z1_ - gamma * real(oneProjector( gamma^(-1) * z1_, param.l1bound )) ;

z2_ = z2 + gamma * Lsx ;
z2 = z2_ ; % Projection onto singleton {0} 

x_ = x - tau* ( (x-xbar) + param.Psi(2*z1-z1old) + Lst(2*z2-z2old) ) ;
x = max(x_, 0) ;

crit(it+1) =sum(abs(x(:)-xbar(:)).^2) ;
Psix = param.Psit(x) ;
l1_cond(it+1) = sum(abs(Psix)) ;
Lsx = Ls(x) ;
smooth_cond(it+1) = min(abs(Lsx)) ;
pos_cond(it+1) = min(x(:)) ;

if mod(it, param.display)==0
    figure(101)
    subplot 221
    hold off, plot(crit,'k'), 
    xlabel('it'), ylabel('crit')
    subplot 222
    hold off, plot(l1_cond,'r'), hold on, plot(param.l1bound*ones(size(l1_cond)),'b--')
    xlabel('it'), ylabel('l1 ball constraint')
    subplot 223
    hold off, semilogy(smooth_cond,'r'), hold on, semilogy(0*ones(size(smooth_cond)),'b--')
    xlabel('it'), ylabel('smooth cond')
    subplot 224
    hold off, plot(pos_cond,'r'), hold on, plot(0*ones(size(pos_cond)),'b--')
    xlabel('it'), ylabel('min x')
    pause(1)
end

if norm(x(:)-xold(:))<1e-4*norm(x(:)) ...   % variation iterates
        && abs(crit(it+1)-crit(it))<1e-4*crit(it+1) ...  % variation critere
        && l1_cond(it+1) - param.l1bound<1e-4 ...   % l1 ball constraint
        && pos_cond(it+1) >= 0 ...                % positivity constraint
        && abs(smooth_cond(it+1)) < 1e-4   % smoothness constraint
    disp(['inpainting projection, it max: ',num2str(it)])
    break;
end

end


%% Final display 

figure(101)
subplot 221
hold off, plot(crit,'k'), 
xlabel('it'), ylabel('crit')
subplot 222
hold off, plot(l1_cond,'r'), hold on, plot(param.l1bound*ones(size(l1_cond)),'b--')
xlabel('it'), ylabel('l1 ball constraint')
subplot 223
hold off, semilogy(smooth_cond,'r'), hold on, semilogy(0*ones(size(smooth_cond)),'b--')
xlabel('it'), ylabel('smooth cond')
subplot 224
hold off, plot(pos_cond,'r'), hold on, plot(0*ones(size(pos_cond)),'b--')
xlabel('it'), ylabel('min x')
pause(1)

x = reshape(x, size(param.Mask,1), size(param.Mask,2)) ;

end




%%
function x = Lt(z,L,Mask)
x = 0*Mask ;
x(Mask>0) = -z ;
x(Mask==0) = L'*z ;
end


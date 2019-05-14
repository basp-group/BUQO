function [x,v1,v2] = projection_HPD_constrainedPB(xbar, v1,v2, param, param_data)

if ~isfield(param,'display'), param.display = 100 ; end
%%

gamma_ = 1 ;


param.beta = 1e-4 * sqrt(2*(param.normPsi+param_data.normPhi)) ;
gamma = gamma_/sqrt(2*(param.normPsi+param_data.normPhi)) ;
tau = 1.99 / (param.beta + gamma_*sqrt(2*(param.normPsi+param_data.normPhi))) ;

if tau * (param.beta + 2*gamma * (param.normPsi+param_data.normPhi)) >=2
    disp('Error convergence parameters for HPD constraint')
end

%%
x=param.x0;
crit(1) = sum(abs(x(:)-xbar(:)).^2) ;
Phix = param_data.Phi(x) ;
Psix = param.Psit(x) ;
condl2(1) = sqrt(sum(abs( Phix(:)-param_data.y(:) ))) ;
condl1(1) = sum(abs(Psix)) ;
condpos(1) = min(x(:)) ;


%%


for it = 1:param.NbIt
Psixold = Psix ;
Phixold = Phix ;
xold = x ;
x = max( x - tau * (param.beta * (x-xbar) + param.Psi(v1) + real(param_data.Phit(v2))) , 0 ) ;

Psix = param.Psit(x) ;
v1_ = v1 + gamma * ( 2*Psix - Psixold ) ;
v1 = v1_ - gamma * real(oneProjector( gamma^(-1) * v1_, param_data.HPDconstraint )) ;

Phix = param_data.Phi(x) ;
v2_ = v2 + gamma * ( 2*Phix - Phixold ) ;
v2 = v2_ - gamma * proj_l2ball( gamma^(-1) * v2_, param_data.l2bound, param_data.y ) ;


crit(it+1) =sum(abs(x(:)-xbar(:)).^2) ;
condl2(it+1) = sqrt(sum(abs( Phix(:)-param_data.y(:) ).^2)) ;
condl1(it+1) = sum(abs(Psix)) ;
condpos(it+1) = min(x(:)) ;

% if mod(it, param.display)==0
%     figure(101)
%     subplot 221
%     hold off, plot(crit,'k'), 
%     xlabel('it'), ylabel('crit')
%     subplot 222
%     hold off, plot(condl2,'r'), hold on, plot(param_data.l2bound*ones(size(condl2)),'b--')
%     xlabel('it'), ylabel('l2 ball')
%     subplot 223
%     hold off, plot(condl1,'r'), hold on, plot(param_data.HPDconstraint*ones(size(condl1)),'b--')
%     xlabel('it'), ylabel('l1 ball - HPD')
%     subplot 224
%     hold off, plot(condpos,'r'), hold on, plot(0*ones(size(condpos)),'b--'), 
%     xlabel('it'), ylabel('min x')
%     pause(1)
% end

if norm(x(:)-xold(:))<param.cond_stop*norm(x(:)) ...   % variation iterates
        && abs(crit(it+1)-crit(it))<param.cond_stop*crit(it+1) ...  % variation critere
        && condl2(it+1) - param_data.l2bound<param.cond_stop ...   % l2 ball constraint
        && condpos(it+1) >= 0 ...                % positivity constraint
        && condl1(it+1) - param_data.HPDconstraint < param.cond_stop   % HPD constraint
    disp(['Projection onto Calpha, it max: ',num2str(it)])
    break;
end

end

%%

% figure(101)
% subplot 221
% hold off, plot(crit,'k'), 
% xlabel('it'), ylabel('crit')
% subplot 222
% hold off, plot(condl2,'r'), hold on, plot(param_data.l2bound*ones(size(condl2)),'b--')
% xlabel('it'), ylabel('l2 ball')
% subplot 223
% hold off, plot(condl1,'r'), hold on, plot(param_data.HPDconstraint*ones(size(condl1)),'b--')
% xlabel('it'), ylabel('l1 ball - HPD')
% subplot 224
% hold off, plot(condpos,'r'), hold on, plot(0*ones(size(condpos)),'b--') 
% xlabel('it'), ylabel('min x')
% pause(1)


end




function p = proj_l2ball(x, eps, y)
% projection of x onto the l2 ball centered in y with radius eps
p = x-y ;
p = p* min(eps/norm(p(:)),1) ;
p = p+y ;

end

function [result] = BUQO_struct(xmap,xmap_S, param_algo, param_data, param_C, param_struct)



%% Initialisation

x = xmap ;


% projections in S and C
xS = xmap_S ;
xC = xmap_S ;

% initialisation sub-it for set S
v = 0*param_struct.L*xC(param_struct.Mask==0) ;

% Information energy in the structure
En_struct_mask = sqrt(sum( abs( xmap(param_struct.Mask>0) - xmap_S(param_struct.Mask>0) ).^2 )) ;
En_xC_mask = sqrt(sum( abs( xmap(param_struct.Mask>0) ).^2 )) ;
disp('**************************************************************')
disp('Initialization')
disp('--------------------------------------')
disp(['Energy M(xmap)    = ',num2str(En_xC_mask)])
disp(['Energy M(xmap_S)  = ',num2str(sqrt(sum( abs( xmap_S(param_struct.Mask>0) ).^2 )))])
disp('--------------------------------------')
disp(['HPD constraint    = ',num2str(param_data.HPDconstraint)])
disp(['l1 norm (xmap)    = ', num2str( sum(abs(param_C.Psit(xmap))) )])
disp(['l1 norm (xmap_S)  = ', num2str( sum(abs(param_C.Psit(xmap_S))) )])
disp('--------------------------------------')
disp(['l2 bound          = ',num2str(param_data.l2bound)])
disp(['l2 norm (xmap)    = ',num2str( norm( param_data.Phi(xmap) - param_data.y) )])
disp(['l2 norm (xmap_S)  = ',num2str( norm( param_data.Phi(xmap_S) - param_data.y) )])
disp('**************************************************************')
disp('**************************************************************')
disp(' ')

display_images(1000, 0, xmap, xS, xC) 

%% Iterations 

for it_glob = 1:param_algo.NbIt
    
xS_old = xS ; 
xC_old = xC ;


%% Projection onto set without strong structure 
%  (smooth image into the Mask)



xC_out = xC(param_struct.Mask==0) ;
xC_in  = xC(param_struct.Mask>0) ;

if it_glob > 1
    disp('***************************************')
    disp('Projection onto S ...')
    param_struct.x0 = xC_out ; 
    [xS,v] = projection_inpaint_no_struct(xC_in, xC_out, v, param_struct) ;
    disp('... done')
    disp('***************************************')
end


% Saving data
En_xS_mask = sqrt(sum( abs( xS(param_struct.Mask>0) ).^2 )) ;
Phi_xinp = param_data.Phi(xS) ;
condl2(it_glob) = sqrt(sum(abs(Phi_xinp(:)-param_data.y(:)).^2)) ;
crit_map(it_glob) = sum(abs(param_C.Psit(xS))) ;

% Display
disp('***************************************')
disp(['it (xS): ', num2str(it_glob)])
disp('--------------------------------------')
disp(['Energy M(xS)      = ',num2str(En_xS_mask)])
disp(['Energy M(xC)      = ',num2str(En_xC_mask)])
disp('--------------------------------------')
disp(['|| Phi(x) - y ||  = ',num2str(condl2(it_glob)),' vs. l2 bound: ',num2str(param_data.l2bound)])
disp(['|| Psit(x) ||_1   = ',num2str(crit_map(it_glob)),' vs. HPD bound: ',num2str(param_data.HPDconstraint)])
disp('--------------------------------------')
disp(['d(xinp, xhpd)     = ', num2str(norm(xS(:)-xC(:)))])
disp('***************************************')

% stopping criteria
if crit_map(it_glob) - param_data.HPDconstraint < 1e-5 ...
        && condl2(it_glob) - param_data.l2bound < 1e-5
    disp('Intersection nonempty: xS in Calpha')
    dist(it_glob) = norm(xS(:)-xC(:)) ;
    x = xS ;
    break;
else
    display_images(1000, it_glob, xmap, xS, xC)
end


%% Projection into the HPD set

param_C.x0 =  xS ;
if it_glob == 1
    v1 = 0*param_C.Psit(xS) ;
    v2 = 0*param_data.y ;
end
disp('***************************************')
disp('Projection onto Calpha ...')
[xC,v1,v2] = projection_HPD_constrainedPB(xS, v1,v2, param_C, param_data) ;
disp('... done')
disp('***************************************')

% Saving data
hpd_mask = xC ; hpd_mask(param_struct.Mask==0) = 0 ;
En_xC_mask = sqrt(sum( abs( hpd_mask(:) ).^2 )) ;


%% Updates

dist(it_glob) = norm(xS(:)-xC(:)) ;
rho_alpha(it_glob) = dist(it_glob) / En_struct_mask ;
norminp(it_glob) = norm(xS(:)-xS_old(:)) / norm(xS(:)) ;
normhpd(it_glob) = norm(xC(:)-xC_old(:)) /  norm(xC(:)) ;
if it_glob>1
cond_dist = abs(dist(it_glob)-dist(it_glob-1)) / dist(it_glob) ;
end

disp(' ')
disp('***************************************')
disp(['it (xC): ', num2str(it_glob)])
disp('--------------------------------------')
disp(['Energy M(xS)      = ',num2str(En_xS_mask)])
disp(['Energy M(xC)      = ',num2str(En_xC_mask)])
disp('--------------------------------------')
disp(['|| Phi(x) - y ||  = ',num2str(condl2(it_glob)),' vs. l2 bound: ',num2str(param_data.l2bound)])
disp(['|| Psit(x) ||_1   = ',num2str(crit_map(it_glob)),' vs. HPD bound: ',num2str(param_data.HPDconstraint)])
disp('--------------------------------------')
disp(['d(xinp, xhpd) = ', num2str(dist(it_glob))])
disp(['rho_alpha = ', num2str(rho_alpha(it_glob))])
if it_glob>1
disp('--------------------------------------')
disp('Conditions for stopping algorithm:')
disp(['variations xhpd = ', num2str(normhpd(it_glob))])
disp(['variations xinp = ', num2str(norminp(it_glob))])
disp(['     vs. stop cond = ',num2str(param_algo.stop_norm)])
disp(['variation distance = ', num2str(cond_dist)])
disp(['     vs. stop cond = ',num2str(param_algo.stop_dist)])
end
display_images(1000, it_glob, xmap, xS, xC)
disp('**************************************************************')
disp('**************************************************************')
disp(' ')

%Stopping criteria
if it_glob>1 && cond_dist < param_algo.stop_dist 
    disp(['d(S, C) converged, with cond_dist=',num2str(cond_dist)])
    x = xC ;
    break;
end
if it_glob>1 && norminp(it_glob) < param_algo.stop_norm ...
        && normhpd(it_glob) < param_algo.stop_norm
    x = xC ;
    disp('xS and xC converged')
    disp(['var. xS = ',num2str(norminp(it_glob))])
    disp(['var. xC = ',num2str(normhpd(it_glob))])
    break;
end    


end


%% Results

display_images(1000, it_glob, xmap, xS, xC)

result.xS = xS ;
result.xC = xC ;
result.x = x ;
result.dist = dist ;
result.rho_alpha = rho_alpha ;
result.crit_map = crit_map ;

result.norminp = norminp ;
result.normhpd = normhpd ;
end



%% 

function display_images(no, it, xmap, xinp, xhpd)


figure(no)
subplot 131
imagesc(log10(xmap)), axis image; colorbar, colormap jet, caxis([0,max(xmap(:))]), xlabel('x map')
subplot 132
imagesc(log10(xinp)), axis image; colorbar, colormap jet, caxis([0,max(xinp(:))]), xlabel(['xS - it ',num2str(it)])
subplot 133
imagesc(log10(xhpd)), axis image; colorbar, colormap jet, caxis([0,max(xhpd(:))]), xlabel(['xC - it ',num2str(it)])
pause(1)

end

%%









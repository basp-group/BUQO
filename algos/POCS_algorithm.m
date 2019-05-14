function [result] = POCS_algorithm(xmap,xmap_inp, param_algo, param_data, param_hpd, param_struct)



%% Initialisation
flag_dual_init = 0 ;

x = xmap ;

if strcmp(param_struct.choice_set,'back0')==0
gam = 1.9/(param_struct.normL+1) ;
end


xinp = xmap_inp ;
xhpd = xmap_inp ;
% save iterates
Xinp(1:size(xinp,1), 1:size(xinp,2), 1) = xinp ;
Xhpd(1:size(xhpd,1), 1:size(xhpd,2), 1) = xhpd ;


map_mask = xmap ; map_mask(param_struct.Mask==0) = 0 ;
map_inp_mask = xmap_inp ; map_inp_mask(param_struct.Mask==0) = 0 ;
En_map_mask = sqrt(sum( abs( map_mask(:) ).^2 )) ;
En_map_inp_mask = sqrt(sum( abs( map_inp_mask(:) ).^2 )) ;
En_hpd_mask = En_map_mask ;
En_inp_mask = En_map_inp_mask ;
disp('**************************************************************')
disp('Initialization')
disp(['Energy M(xmap)     = ',num2str(En_map_mask)])
disp(['Energy M(sm. xmap) = ',num2str(En_map_inp_mask)])
disp(['HPD constraint    = ',num2str(param_data.HPDconstraint)])
disp(['l1 norm (map)     = ', num2str( sum(abs(param_hpd.Psit(xmap))) )])
disp(['l1 norm (sm. map) = ', num2str( sum(abs(param_hpd.Psit(xmap_inp))) )])
disp(['l2 bound          = ',num2str(param_data.l2bound)])
disp(['l2 norm (map)     = ',num2str( norm( param_data.Phi(xmap) - param_data.y) )])
disp(['l2 norm (sm. map) = ',num2str( norm( param_data.Phi(xmap_inp) - param_data.y) )])
disp('**************************************************************')
disp('**************************************************************')
disp(' ')

display_images(1000, 0, xmap, xinp, xhpd) 

%% Iterations 

for it_glob = 1:param_algo.NbIt
    
xinp_old = xinp ; 
xhpd_old = xhpd ;


%% Projection onto set without strong structure 
%  (smooth image into the Mask)



xhpd_out = xhpd(param_struct.Mask==0) ;
xhpd_in  = xhpd(param_struct.Mask>0) ;

if it_glob > 1
switch param_struct.choice_set
    case 'smooth'
        xinp = proj_inpaint_smooth(xhpd_out, xhpd_in, param_struct.L, gam, param_struct.Mask) ;
    case 'pos'
        xinp = proj_inpaint_pos(xhpd_out, xhpd_in, param_struct.L, gam, param_struct.Mask) ;
    case 'l1'
        param_struct.x0 = proj_inpaint_pos(xhpd_out, xhpd_in, param_struct.L, gam, param_struct.Mask) ;
        if flag_dual_init == 0
            z1 = 0*param_hpd.Psit(xinp) ;
            z2 = 0*param_struct.Lbar*xinp(:) ;
            flag_dual_init = 1 ;
        end
        [xinp,z1,z2] = projection_inpaint_constrainedPB(xhpd, z1,z2, param_struct) ;
    case 'l2_const'
        param_struct.x0 = xhpd_out ; 
        if flag_dual_init == 0
            v = 0*param_struct.L*xhpd_out ;
            flag_dual_init = 1 ;
        end
        [xinp,v] = projection_inpaint_no_struct(xhpd_in, xhpd_out, v, param_struct) ;
    case 'back0'
        xinp = proj_back_zero(xhpd, param_struct.Mask) ;
end
end


% Saving data
Xinp(1:size(xhpd,1), 1:size(xhpd,2), it_glob) = xinp ;
inp_mask = xinp ; inp_mask(param_struct.Mask==0) = 0 ;
En_inp_mask = sqrt(sum( abs( inp_mask(:) ).^2 )) ;
Phi_xinp = param_data.Phi(xinp) ;
condl2(it_glob) = sqrt(sum(abs(Phi_xinp(:)-param_data.y(:)).^2)) ;
crit_map(it_glob) = sum(abs(param_hpd.Psit(xinp))) ;

% Display
disp('***************************************')
disp(['it (inp): ', num2str(it_glob)])
disp(['Energy M(xmap)     = ',num2str(En_map_mask)])
disp(['Energy M(sm. xmap) = ',num2str(En_map_inp_mask)])
disp(['Energy M(xinp) = ',num2str(En_inp_mask)])
disp(['Energy M(xhpd) = ',num2str(En_hpd_mask)])
disp(['|| Phi(x) - y || = ',num2str(condl2(it_glob)),' vs. l2 bound: ',num2str(param_data.l2bound)])
disp(['|| Psit(x) ||_1 = ',num2str(crit_map(it_glob)),' vs. HPD bound: ',num2str(param_data.HPDconstraint)])
disp(['d(xinp, xhpd) = ', num2str(norm(xinp(:)-xhpd(:)))])
disp('***************************************')

% stopping criteria
if crit_map(it_glob) - param_data.HPDconstraint < 1e-5 ...
        && condl2(it_glob) - param_data.l2bound < 1e-5
    disp('intersection nonempty: inpainted image inside the HPD region')
    dist(it_glob) = norm(xinp(:)-xhpd(:)) ;
    x = xinp ;
    break;
else
    disp('inpainted image outside the HPD region')
    display_images(1000, it_glob, xmap, xinp, xhpd)
end


%% Projection into the HPD set

param_hpd.x0 =  xinp ;
if it_glob == 1
    v1 = 0*param_hpd.Psit(xinp) ;
    v2 = 0*param_data.y ;
end
[xhpd,v1,v2] = projection_HPD_constrainedPB(xinp, v1,v2, param_hpd, param_data) ;

% Saving data
Xhpd(1:size(xhpd,1), 1:size(xhpd,2), it_glob) = xhpd ;
hpd_mask = xhpd ; hpd_mask(param_struct.Mask==0) = 0 ;
En_hpd_mask = sqrt(sum( abs( hpd_mask(:) ).^2 )) ;


%% Updates

dist(it_glob) = norm(xinp(:)-xhpd(:)) ;
norminp(it_glob) = norm(xinp(:)-xinp_old(:)) / norm(xinp(:)) ;
normhpd(it_glob) = norm(xhpd(:)-xhpd_old(:)) /  norm(xhpd(:)) ;

disp(' ')
disp('***************************************')
disp(['it (hpd): ', num2str(it_glob)])
disp('- - - - - - - - - - - - - - - - - - - -')
disp(['Energy M(xmap)     = ',num2str(En_map_mask)])
disp(['Energy M(sm. xmap) = ',num2str(En_map_inp_mask)])
disp(['Energy M(xinp) = ',num2str(En_inp_mask)])
disp(['Energy M(xhpd) = ',num2str(En_hpd_mask)])
disp('- - - - - - - - - - - - - - - - - - - -')
disp(['|| Phi(x) - y || = ',num2str(condl2(it_glob)),' vs. l2 bound: ',num2str(param_data.l2bound)])
disp(['|| Psit(x) ||_1 = ',num2str(crit_map(it_glob)),' vs. HPD bound: ',num2str(param_data.HPDconstraint)])
disp('- - - - - - - - - - - - - - - - - - - -')
disp(['d(xinp, xhpd) = ', num2str(dist(it_glob))])
disp(['variations xinp = ', num2str(norminp(it_glob))])
disp(['variations xhpd = ', num2str(normhpd(it_glob))])
display_images(1000, it_glob, xmap, xinp, xhpd)
disp('**************************************************************')
disp('**************************************************************')
disp(' ')

%Stopping criteria
if it_glob>1
cond_dist = abs(dist(it_glob)-dist(it_glob-1)) / dist(it_glob) ;
end
if it_glob>1 && cond_dist < param_algo.stop_dist 
    disp(['d(xinp, xhpd) converged, with cond_dist=',num2str(cond_dist)])
    x = xhpd ;
    break;
end
if it_glob>1 && norminp(it_glob) < param_algo.stop_norm ...
        && normhpd(it_glob) < param_algo.stop_norm
    x = xhpd ;
    disp('xinp and xhpd converged')
    disp(['var. inp = ',num2str(norminp(it_glob))])
    disp(['var. hpd = ',num2str(normhpd(it_glob))])
    break;
end    


end


%% Results

display_images(1000, it_glob, xmap, xinp, xhpd)

result.xinp = xinp ;
result.xhpd = xhpd ;
result.x = x ;
result.dist = dist ;
result.crit_map = crit_map ;
result.Xinp = Xinp ;
result.Xhpd = Xhpd ;
result.norminp = norminp ;
result.normhpd = normhpd ;
end



%% 

function display_images(no, it, xmap, xinp, xhpd)


figure(no)
subplot 231
imagesc(xmap), axis image; colorbar, colormap jet, caxis([0,max(xmap(:))]), xlabel('map')
subplot 232
imagesc(xinp), axis image; colorbar, colormap jet, caxis([0,max(xinp(:))]), xlabel(['S it ',num2str(it)])
subplot 233
imagesc(xhpd), axis image; colorbar, colormap jet, caxis([0,max(xhpd(:))]), xlabel(['Calpha it ',num2str(it)])
subplot 234
imagesc(log10(xmap)), axis image; colorbar, colormap jet, caxis([-3.5, log10(max(xmap(:)))]), xlabel('map log')
subplot 235
imagesc(log10(xinp)), axis image; colorbar, colormap jet, caxis([-3.5, log10(max(xinp(:)))]), xlabel('S log')
subplot 236
imagesc(log10(xhpd)), axis image; colorbar, colormap jet, caxis([-3.5, log10(max(xhpd(:)))]), xlabel('Calpha log')
pause(1)

end

%%

function xinp = proj_inpaint_smooth(xhpd_out, xhpd_in, L, gam, Mask)

xout = xhpd_out ;
for ii = 1:500
xout_ = xout ;
xout = xout - gam * L'*(L*xout - xhpd_in) ...
            - gam * (xout - xhpd_out) ;
if norm(xout(:)-xout_(:))<1e-6 * norm(xout(:)) 
disp(['init sub-it: ', num2str(ii)])
    break;
end
end
xinp_out = xout ;
xinp(Mask==0) = xinp_out ; 
xinp(Mask>0) = L * xinp_out ; 

xinp = reshape(xinp, size(Mask,1), size(Mask,2)) ;

end

function xinp = proj_inpaint_pos(xhpd_out, xhpd_in, L, gam, Mask)

xout = xhpd_out ;
for ii = 1:500
xout_ = xout ;
xout = xout - gam * L'*(L*xout - xhpd_in) ...
            - gam * (xout - xhpd_out) ;
xout = max(xout,0) ;
if norm(xout(:)-xout_(:))<1e-6 * norm(xout(:)) 
disp(['init sub-it: ', num2str(ii)])
    break;
end
end
xinp_out = xout ;
xinp(Mask==0) = xinp_out ; 
xinp(Mask>0) = L * xinp_out ; 

xinp = reshape(xinp, size(Mask,1), size(Mask,2)) ;

end


function xinp = proj_back_zero(xhpd, Mask)
xinp = xhpd ;
xinp(Mask>0) = 0 ;
end









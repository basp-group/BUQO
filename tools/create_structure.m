function [param_struct, xmap_S] = create_structure(xmap, cropx, cropy, mincrop, maxcrop)



if maxcrop<Inf  % background
Mask = ones(size(xmap)) ;
Mask(xmap<maxcrop) = 0 ;
Mask = imdilate(Mask,strel('disk',7));
Mask = abs(Mask-1);
else            % structures
Mask = zeros(size(xmap));
Mask(cropy,cropx) = 1;
Mask = Mask.*(xmap>mincrop);
Mask = Mask.*(xmap<maxcrop);
Mask = imdilate(Mask,strel('disk',3));
end
param_struct.Mask = Mask ;

tmp = xmap ;
tmp(Mask>0) = 0 ;
figure,
subplot 221, imagesc(xmap), axis image; colorbar, caxis([0,1]), colormap jet
xlabel('xmap (lin scale)')
subplot 222, imagesc(log10(xmap)), axis image; colorbar, caxis([-3.5 0]), colormap jet
xlabel('xmap (log scale)')
subplot 223, imagesc(log10(tmp)), axis image; colorbar, caxis([-3.5 0]), colormap jet
xlabel('xmap without struct (log scale)')
subplot 224, imagesc(Mask), axis image; colorbar
xlabel('mask')

  


if maxcrop<Inf      % structures
param_struct.choice_set = 'back0' ;

xmap_S = xmap ;
xmap_S(param_struct.Mask>0) = 0 ;
figure, imagesc(log10(xmap_S)), axis image ; colorbar, colormap jet, caxis([-4.5 0])
else                % background
param_struct.choice_set = 'l2_const' ;
param_struct.Size_Gauss_kern = [3,7,11] ;
disp('create inpainting operator...')
disp('     -- recursive - can take few minutes... --')
param_struct.L = create_inpainting_operator_test(param_struct.Mask, param_struct.Size_Gauss_kern, xmap) ;
disp('...done')


param_struct.Lbar = sparse([-speye(sum(param_struct.Mask(:))), param_struct.L]) ;
param_struct.Nout = sum(param_struct.Mask(:)==0) ;
param_struct.normL = op_norm(@(x) param_struct.L*x, @(x) param_struct.L'*x, [param_struct.Nout,1], 1e-4, 200, 0);  
param_struct.normLbar = op_norm(@(x) param_struct.Lbar*x, @(x) param_struct.Lbar'*x, [numel(Mask),1], 1e-4, 200, 0);  

param_struct.Li =@(u) [param_struct.L*u ; u] ;
param_struct.Lit =@(x) param_struct.L'*x(param_struct.Mask>0) + x(param_struct.Mask==0) ;
param_struct.normLi = op_norm(param_struct.Lit, param_struct.Li, [numel(param_struct.Mask),1], 1e-4, 200, 0);  

xmap_S = xmap ;
xmap_S(param_struct.Mask>0) = param_struct.L*xmap(param_struct.Mask==0) ;
figure, imagesc(log10(xmap_S)), axis image ; colorbar, colormap jet, caxis([-3.5 0])

param_struct.l2_mean = 0 ;
param_struct.l2_bound = sqrt(sum(abs(xmap_S(param_struct.Mask>0)).^2)) ; % 1e-3 * sqrt(sum(param_inpaint.Mask(:))) ; 
end


param_struct.NbIt = 3000 ;
param_struct.display = 10 ;



end
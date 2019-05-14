function [L, xinp_tot,xinp_multscale]  = create_inpainting_operator_test(Mask, Size_Gauss_kern, xmap)



%% Initialization

Nx = size(Mask,2) ; Ny = size(Mask,1) ;

% Define the variance for the Gaussian weights
Sigma_Gauss_kern = 0.5* sqrt(Size_Gauss_kern) ; 

% save the inpainted images for each scale
xinp_multscale = cell(1,1,length(Size_Gauss_kern)) ;

%% create matrix

Nm = sum(Mask(:)) ; % dim variable inside mask
Nc = Nx*Ny - Nm ;   % dim variable outside mask

L = sparse(Nm, Nc) ;
s = 0 ;

%%
for Sgauss = Size_Gauss_kern
s = s+1 ;
% build the Gaussian weights
sigmgauss = Sigma_Gauss_kern(s) ; 
wG = fspecial('gaussian',Sgauss*[1,1],sigmgauss) ;
wG( floor(Sgauss/2)+1, floor(Sgauss/2)+1) = 0 ;
wG = wG/( wG( floor(Sgauss/2), floor(Sgauss/2)+1)) ;

Maskt = Mask ;
xinp = xmap ; xinp(Mask>0)=0 ;
ro = 0 ;
indx = [] ; indy = [] ;
Sbord = 0 ;
% edges of Mask at last iteration
indx_old = indx ;
indy_old = indy ;
Sbord_old = Sbord ;

Ltmp = sparse(Nm, Nc) ;


%%
while sum(Maskt(:))>0
ro = ro+1 ;


%% mask update

% search indices in the Mask
indx_ = [] ;
indy_ = [] ;
for i = 1:Nx
for j = 1:Ny
if Maskt(j,i) ==1
indx_= [indx_; i] ;
indy_= [indy_; j] ;
end
end
end
Ind = [indy_,indx_];

% search the edges of the mask
Maskx = Maskt ; Maskx(2:end, 2:end) = Maskt(2:end,2:end) - Maskt(1:end-1,2:end) ;
for i=1:Nx
for j=1:Ny
if Maskx(j,i) == -1
Maskx(j,i) = 0 ;
Maskx(j-1,i) = 1 ;
end
end
end
Masky = Maskt ; Masky(2:end, 2:end) = Maskt(2:end,2:end) - Maskt(2:end,1:end-1) ;
for i=1:Nx
for j=1:Ny
if Masky(j,i) == -1
Masky(j,i) = 0 ;
Masky(j,i-1) = 1 ;
end
end
end
Maskg = min(abs(Maskx)+abs(Masky),1);
Sbord = sum(Maskg(:)) ; % number of elements in the edges

% search indices of the edges of the Mask
indx = [] ;
indy = [] ;
for i = 1:Nx
for j = 1:Ny
if Maskg(j,i) ==1
indx= [indx, i] ;
indy= [indy, j] ;
end
end
end

xtmp = xmap ; xtmp(Maskt>0) = 0 ;


%%
for t=1:Sbord % each pixel in the edges of the mask
    
i=indx(t);
j=indy(t);

%%
% truncated neighborhood if the pixel of interest is too close to the edge 
% of the global image
Voisx = max(i - floor(Sgauss/2), 1) : min(i + floor(Sgauss/2), Nx) ;
Voisy = max(j - floor(Sgauss/2), 1) : min(j + floor(Sgauss/2), Ny) ;
Vois = [] ;
for u = 1:length(Voisx)
Vois = [Vois ; [Voisy', Voisx(u)*ones(length(Voisy),1)] ] ;
end

% keep neighboors only outside the Mask
ind_rem = ismember(Vois, Ind, 'rows') ;
Vois(ind_rem,:) = [] ;

%%
% truncated weights if the pixel of interest is too close to the edge of
% the global image
wGtmp = wG ;
if i - floor(Sgauss/2)<1
wGtmp(:,1:abs(i - floor(Sgauss/2)+1)) = [] ;
end
if j - floor(Sgauss/2)<1
wGtmp(1:abs(j - floor(Sgauss/2)+1),:) = [] ;
end
if i + floor(Sgauss/2)>Nx
wGtmp(:,end-abs(i + floor(Sgauss/2)-Nx)+1:end) = [] ;
end
if j + floor(Sgauss/2)>Ny
wGtmp(end-abs(j + floor(Sgauss/2)-Ny)+1:end,:) = [] ;
end
WGt = zeros(size(Mask)) ;
WGt(Voisy,Voisx) = wGtmp/sum(wGtmp(:)) ;

WG = zeros(size(Mask)) ;
for u = 1:numel(Vois(:,1)) 
WG(Vois(u,1), Vois(u,2)) = WGt(Vois(u,1), Vois(u,2)) ;
end
WG = WG / sum(WG(:)) ;

%%
% weighted mean value
% % % values = [] ;
% % % sum_weigths = 0 ;
% % % for u = 1:numel(Vois(:,1)) 
% % % values = [values, WGt(Vois(u,1), Vois(u,2)) * xtmp(Vois(u,1), Vois(u,2))] ;
% % % sum_weigths = sum_weigths + WGt(Vois(u,1), Vois(u,2)) ;
% % % end
% % % xinp(j,i) = sum( values )/sum_weigths ;

% % % values = WG.*xtmp ;
% % % xinp(j,i) = sum( values(:) )/sum(WG(:)) ;


%% remove coefficients from WG within the mask
%  (write them w.r.t. the other coefficients)

[S_up, iv, isup] = intersect(Vois, [indy_old', indx_old'], 'rows', 'legacy') ;

for tt = 1:size(S_up,1)    % go through all the neighbors within the mask already treated
    ii = S_up(tt,2) ; %indx_old(isup(tt)) ;
    jj = S_up(tt,1) ; %indy_old(isup(tt)) ;
    % save the weight
    wt = WG(jj,ii) ; 
    % find the corresponding pixel n in L(n,:)
    tmp = zeros(Ny,Nx) ;
    tmp(jj,ii) = 1 ;
    tmp = tmp(Mask>0) ;
    tmp2 = 1:size(tmp,1) ;
    n = tmp2(tmp>0) ;    
    WGt = zeros(Ny*Nx,1) ;
    WGt(Mask==0) = Ltmp(n,:) ;
    WGt = reshape(WGt,Ny, Nx) ;
    % update WG
    WG(jj,ii) = 0 ;  
    WG = WG + wt * WGt ;
end


values = WG.*xtmp ;
xinp(j,i) = sum( values(:) ) ;


%% build matrix for coeff (i,j)

WGt = WG(Mask==0) ;

% index correspondance within the mask
tmp = zeros(Ny,Nx) ;
tmp(j,i) = 1 ;
tmp = tmp(Mask>0) ;
tmp2 = 1:size(tmp,1) ;
n = tmp2(tmp>0) ;

Ltmp(n,:) = sparse(WGt)' ;




end



%%
% update mask
Maskt = Maskt - Maskg ;
% edges of Mask at last iteration
indx_old = [indx_old, indx] ;
indy_old = [indy_old, indy] ;
Sbord_old = Sbord_old+Sbord ;

end

%% update multiscaled inpainted images 
xinp_multscale{s} = xinp ;
L = L+Ltmp ;

end
%% global inpainted image
xinp_tot_ = (1/length(Size_Gauss_kern)) * squeeze(sum(cell2mat(xinp_multscale),3)) ;
xinp_tot = xmap ;
xinp_tot(Mask>0) = xinp_tot_(Mask>0) ;


L = L/length(Size_Gauss_kern) ;

end
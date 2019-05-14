%  POCS for uncertainty quantification
%  Solving for the constrained problem
%  minimize || Psit(x) ||_1 
%      s.t. x positive
%       and || Phi(x) - y ||_2 <= epsilon
%  ------------------------------------------------------------------------
%% ------------------------------------------------------------------------

clear all
clc
close all


addpath(genpath('data'))
addpath(genpath('algos'))
addpath(genpath('tools'))

% Necessitate the use of irt package
% Download at http://web.eecs.umich.edu/~fessler/irt/fessler.tgz
irtdir = '../irt' ; % SET THE SOPT PATH
setup(irtdir);


%%
SNR =@(x,xtrue) 20 * log10(norm(xtrue(:))/norm(xtrue(:)-x(:)));


%% Load image

load('im_W28_256')

%% Radio-astro measurements


param_data.p = 0.5 ; 
param_data.sigma2 = 0.03 ; 
param_data.sigma_noise = sqrt(param_data.sigma2) ;
param_data.sigma = pi/4 ;
param_data.N = N ;

% Generate Gaussian random u-v sampling
[u, v] = util_gen_sampling_pattern('gaussian', param_data);

% Initialize nuFFT operator
% Generate measurment operator with nufft
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft
[A, AT, Gw, ~] = op_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);

param_data.Phi =@(x) Gw * A(x) ;
param_data.Phit =@(y) AT(Gw' * y) ;

% norm of the measurement operator
param_data.normPhi = op_norm(param_data.Phi, param_data.Phit, [Ny, Nx], 1e-4, 200, 0);    


% Generate noisy measurements
y0 = param_data.Phi(im);
noise = param_data.sigma_noise*(randn(size(y0)) + 1i*randn(size(y0))) ;
param_data.y = y0 + noise;

param_data.M = length(y0) ;

% param_data.l2bound = 1.1 * norm(noise) ;
param_data.l2bound = sqrt(2*param_data.M + 2* sqrt(4*param_data.M)) *  param_data.sigma_noise ;

% Sparsity basis for l1 regularization
[Psi, Psit] = op_sp_wlt_basis({'db8'}, 4, Ny, Nx);
normPsi = 1 ; %if normalized

%% Run MAP estimate



param_algo_map.NbIt = 3000 ;
param_algo_map.stop_crit = 1e-6 ;
param_algo_map.stop_norm = 1e-5 ;
param_algo_map.stop_cond = 1e-6 ;
param_algo_map.Psit = Psit ;
param_algo_map.Psi = Psi ;
param_algo_map.normPsi = normPsi ;
param_algo_map.x0 = 0 * real(param_data.Phit(param_data.y)) ;
param_algo_map.x0 = max(param_algo_map.x0 / max(param_algo_map.x0(:)), 0) ;
param_algo_map.lambda = 1e-3 ; 
param_algo_map.gammat = param_data.normPhi/1e4 ; 
param_algo_map.display = 200 ;

result_MAP = solve_MAP_constrainedPB(param_data,param_algo_map) ;


figure, 
subplot 221
imagesc(log10(result_MAP.x)), axis image; colorbar, colormap jet, caxis([-3.5,0])
xlabel('x MAP (log scale)')
subplot 222
imagesc(log10(im)), axis image; colorbar, colormap jet, caxis([-3.5,0])
xlabel('x true (log scale)')
subplot 223
imagesc(result_MAP.x), axis image; colorbar, colormap jet, caxis([0,1])
xlabel('x MAP (lin scale)')
subplot 224
imagesc(im), axis image; colorbar, colormap jet, caxis([0,1])
xlabel('x true (lin scale)')





%% HPD definition and parameters

xmap = result_MAP.x ;

param_C.lambda_t = param_data.N / sum(abs(param_algo_map.Psit(xmap))) ;


alpha = 1e-2 ; 
talpha = sqrt( (16*log(3/alpha)) / param_data.N );
HPDconstraint = param_C.lambda_t* sum(abs(param_algo_map.Psit(xmap))) ...
                + param_data.N*(1+talpha);
param_data.HPDconstraint = HPDconstraint/param_C.lambda_t ;




param_C.NbIt = 3000 ;
param_C.Psit = param_algo_map.Psit ;
param_C.Psi = param_algo_map.Psi ;
param_C.normPsi = param_algo_map.normPsi ;
param_C.lambda = param_algo_map.lambda ;
param_C.display = 200 ;
param_C.perc = 1 ;
param_C.cond_stop = 1e-4 ;



%%
Test_choice = 'Struct2' ;  
param_algo.test_choice = Test_choice ;

switch Test_choice
    case 'Struct1'      % large faint structure bottom
        cropx = 138 : 178 ;
        cropy = 200 : 230 ;
        mincrop = 10^(-3.2) ;
        maxcrop = Inf ;
    case 'Struct2'      % faint point source left
        cropx = 10 : 20 ;
        cropy = 160 : 175 ;
        mincrop = 10^(-3.0) ;
        maxcrop = Inf ;
    case 'Struct3'      % important point source right
        cropx = 247 : 256 ;
        cropy = 98 : 107 ;
        mincrop = 10^(-3.5) ;
        maxcrop = Inf ;
end

[param_struct, xmap_S] = create_structure(xmap, cropx, cropy, mincrop, maxcrop) ;

%% POCS algorithm

param_algo.NbIt = 5000 ;
param_algo.stop_dist = 5e-6 ;
param_algo.stop_norm = 5e-6 ;



l1_inpaint = sum(abs(param_algo_map.Psit(xmap_S))) ;
l2_inpaint = sqrt(sum( abs( param_data.Phi(xmap_S) - param_data.y ).^2 )) ;


disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp('*******************************************************')
disp('*******************************************************')
disp(['p = ',num2str(param_data.p)])
disp(['sigma2 = ',num2str(param_data.sigma2)])
disp(['l1 inpaint = ',num2str(l1_inpaint)])
disp(['HPD bound = ',num2str(param_data.HPDconstraint)])
disp(['l2 data inpaint = ',num2str(l2_inpaint)])
disp(['data bound = ',num2str(param_data.l2bound)])

if l1_inpaint <= param_data.HPDconstraint && l2_inpaint <= param_data.l2bound
    disp('Intersection between S and Calpha nonempty')
    disp('*******************************************************')
    result.xS = xmap_S ;
    result.xC = xmap_S ;
else
    disp('xmap_S OUTSIDE Calpha -> run alternating projections')
    disp('*******************************************************')
    disp(' ')
    result = BUQO_struct(xmap,xmap_S, param_algo, param_data, param_C, param_struct) ;
end

disp('  ')
disp('  ')
disp('  ')
disp('**********************************************************************')
disp('**********************************************************************')
disp('**********************************************************************')

disp(' ')
disp('*****************************************')
Test_choice
disp(['p = ',num2str(param_data.p)])
disp(['sigma noise = ',num2str(sqrt(param_data.sigma2))])
disp(['rho_alpha = ',num2str(result.rho_alpha(end))])
disp('*****************************************')




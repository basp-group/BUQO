function x = so_fft2_adj(X, N, No, scale)
% Computes the inverse scaled FTT2
%
% in:
% X[:]     - 2D fourier transform of the image
% N[2]     - size of image
% No[2]    - size of zero padded image
% scale[:] - scale factor precomputed by nufft_init
%
% out:
% x[:][:]  - inverse FFT2

% scale factor to cancel the scaling performed asscoiated with the 
% convolution kernel to the non uniform frequency domain
iscale = conj(scale);	

% compute the inverse fourier transform
X = reshape(X, No);
x = ifft2(X);

% scale the solution
x = No(1) * No(2) * x(:);

% reshape to oversampled image size
x = reshape(x, No);

% trim oversampled part to actual image size
x = x(1:N(1), 1:N(2));

% rescale
x = x .* iscale;

end

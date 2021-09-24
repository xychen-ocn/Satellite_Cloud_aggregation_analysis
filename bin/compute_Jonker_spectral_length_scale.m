function [L, LOgive] = compute_Jonker_spectral_length_scale(data_s, a, fthres)

% Reference: Janssens et al. (2021) Supplemental Materials:
% function: wil use the scripts that I had before to compute wavenumber
% spectrum from a 2D image. 

% 1. plane fit the input data to remove tilt. (tilt-compensated by
% subtracting the scene's best-fit plane.
cldfield = data_s.values;

% x = data_s.lon;
% y = data_s.lat;

x = 1:size(cldfield,2);
y = 1:size(cldfield,1);
[XX, YY ] = meshgrid(x,y);
inM = [XX(:), YY(:), cldfield(:)]; 
[n, ~, p]=affine_fit(inM);

% figure
% subplot(1,2,1)
% plot3(XX,YY,cldfield,'.r');
% hold on;
% plot3(p(1),p(2),p(3),'bo','markersize',15,'markerfacecolor','b');
% 
 best_fit_plane =  -(n(1)/n(3)*XX+n(2)/n(3)*YY-dot(n,p)/n(3));
% surf(XX, YY,best_fit_plane, 'facecolor','red','facealpha',0.5);
% axis('square');
% 
% 
% subplot(1,2,2)
 cldfield_ = cldfield - best_fit_plane;
% plot3(XX,YY,cldfield_+ mean(cldfield(:)),'.b');
% axis('square');


% 2. perform 2D FFT and find the omidirectional wavenumber spectrum S(k)
% could call a function or write out the steps. 

% construct cloudmask again
data_s.values = cldfield_ + mean(cldfield(:));

cldmask = construct_cloud_mask(data_s);
cldmask(cldmask<fthres) = 0;
%cldmask = cldmask>fthres;

[k, Sk, LOgive] = compute_wavenumber_spectrum_from_2Dimages(data_s, cldmask);



% 3. compute the spectral length scale as defined in the reference;
valid = ~isnan(Sk) & ~isnan(k);
Lambda = trapz(k(valid), (k(valid).^a).*Sk(valid))./trapz(k(valid), Sk(valid));

L = 1/Lambda;


return 
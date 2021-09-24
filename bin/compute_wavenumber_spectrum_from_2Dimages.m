function [wn_bins, Sk_filled, LOgive] = compute_wavenumber_spectrum_from_2Dimages(datasub, cldmask)

% Input: datasub --> a matlab structure
%        cldmask --> a logical matrix. (1= cloud, 0=not cloud)
% 
midlat = 0.5*(datasub.lat(1)+datasub.lat(end));
Lx = (datasub.lon - datasub.lon(1))*111E3*cosd(midlat);
Ly = (datasub.lat - datasub.lat(1))*111E3;
dx = diff(Lx);
dy = diff(Ly);


%% from here on, the code can be built into a matlab function. 

%% perform fft2:
Y = fft2(cldmask);
S = abs(fftshift(Y)).^2;

% 1. wavenumber spectrum;
[Ny, Nx] = size(cldmask);

% -- 1.a construct wavenumber space:
if mod(Nx,2)==0
  ix = linspace(-Nx/2+1,Nx/2,Nx);
else
  ix = linspace(-floor(Nx/2),floor(Nx/2),Nx);
end

if mod(Ny,2)==0
  iy = linspace(-Ny/2+1,Ny/2,Ny);
else
  iy = linspace(-floor(Ny/2),floor(Ny/2),Ny);
end

kx = 2*pi/(Nx*mean(dx,'omitnan')) .* ix;   % the last kx is the nyquist frequency
ky = 2*pi/(Ny*mean(dy,'omitnan')) .* iy; 

[KX, KY] = meshgrid(kx,ky);


% construct k-theta space, the matrix below is not on a regular grid.
K = sqrt(KX.^2 + KY.^2);
TH = atan2(KY, KX)*180/pi;

% --- bin the spectrum according to wavenumber and direction:
dx_nyquist = 2* sqrt(round(dx(1)/1E3).^2 + round(dy(1)/1E3).^2)*1E3;
k_nyq = 2*pi/dx_nyquist;
L = sqrt((Nx*dx(1)).^2 + (Ny*dy(1)).^2);
k0 =  2*pi/L;
exp_ind = (log10(k0)):0.05:(log10(k_nyq));
wn_edges = 10.^exp_ind;

[Indx_kbined, ~]=discretize(K, wn_edges);
[Indx_thbined, DirEdge]=discretize(TH,36);   % 36 bins;

wn_bins = 0.5*(wn_edges(1:end-1) + wn_edges(2:end));
theta_bins = 0.5*(DirEdge(1:end-1) + DirEdge(2:end));

S_kth = zeros(length(wn_edges)-1, length(DirEdge)-1);
for ik= 1:length(wn_edges)-1
    ids = find(Indx_kbined==ik);
    if ~isempty(ids)
        rec_k(ik) = length(ids);
        Sk_tmp = S(ids);
        theta_tmp = TH(ids);
        K_tmp = K(ids);
        
        for j = 1:length(DirEdge)-1
            Ysub =discretize(theta_tmp, DirEdge);
            for ith = 1:length(DirEdge)-1
                idth = find(Ysub==ith);
                if ~isempty(ith)
                    rec_th(ik,ith) = length(idth);
                    S_kth(ik, ith) = mean(Sk_tmp(idth),'omitnan');
                else
                    rec_th(ik,ith) = 0;
                    S_kth(ik, ith) = NaN;
                end
               
            end
        end
        
    else
        rec_k(ik) = 0;
    end
end


% ---- build a regular grid:
[THgrid, WNgrid] = meshgrid(theta_bins, wn_bins);
mask = ~isnan(S_kth);
theta_bins_valid = THgrid(mask);
wnbins_valid = WNgrid(mask);
S_kth_valid = S_kth(mask);

S_kth_filled = griddata(theta_bins_valid, wnbins_valid,S_kth_valid, THgrid, WNgrid);

% figure(1)
% pcolor(theta_bins, wn_bins, S_kth_filled);shading flat
% colorbar
% set(gca,'yscale','log');


% 3. integrate the 2D spectrum in the directional space:
%    using trapezoidal method:
Sk= trapz(theta_bins, S_kth_filled, 2); 

% do another interpolation:
mask0 = (Sk~=0);
Sk_valid = Sk(mask0);
wn_bins_valid = wn_bins(mask0);
Sk_filled = interp1(log10(wn_bins_valid), Sk_valid, log10(wn_bins),'linear');

Stot=[];
for i=2:length(wn_bins)
    Stot(i) = trapz(wn_bins(1:i), Sk(1:i));
end

%  --- 3.a find the critical wavenumber using Ogive. --- 
[uqid,ia, ic] = unique(log10(Stot));
kc_exponent = interp1(Stot(ia), log10(wn_bins(ia)), 1/3*Stot(end));
kc = 10^kc_exponent;   % 2/3 of the variance comes from wavenumber larger than kc.
LOgive = 2*pi/kc/1E3;


%% make plots to show the spectrum of the input data:
figure(9); clf;
% show the input cloudmask:
hsub(1)=subplot(2,2,1);
[XX, YY] = meshgrid(Lx,Ly);
pcolor(XX, YY, double(cldmask));shading flat;
colormap(hsub(1), gray)
colorbar;

% show the spectrum in kx-ky space:
hsub(2) = subplot(2,2,2);
pcolor(KX, KY, S);shading flat;
ylim([-0.25,0.25]*10^(-3))
xlim([-0.25,0.25]*10^(-3))
xlabel('kx (rad/m)');
ylabel('ky (rad/m)');
colormap(hsub(2), parula)
colorbar;
caxis([0 10^6])



% show the binned averaged and interpolated spectrum in wavenumber
% direction space:
hsub(3) = subplot(2,2,3);
pcolor(theta_bins, wn_bins, S_kth_filled);shading flat
colormap(hsub(3), parula)
colorbar
set(gca,'yscale','log');
caxis([0 10^6])
title('wavenumber-direction spectrum');
xlabel('direction (dgr)')
ylabel('wavenumber k (rad/m)');

% show the omnidirectional spectrum in the wavelength space:
subplot(2,2,4);
y_ideal =50* wn_bins.^(-5/3);
wvlen_bins = 2*pi./wn_bins /1E3;

plot(wvlen_bins, Sk_filled,'-k');
hold on;
hl(1)=plot(wvlen_bins, y_ideal, '--b');
yrange = get(gca,'ylim');
hl(2)=plot(2*pi./[kc,kc]/1E3,[10^5, 10^10],'--r');

set(gca,'xscale','log','yscale','log','xdir','reverse');
%set(gca,'xtick',fliplr([10^3, round(2*pi/kc/1E3), 10^2, 10^1, 10^0]));
text(2*pi/kc/1E3, 10^5*5, num2str(round(2*pi./kc/1E3)));
xlim([10^0, 10^3]);
%xlabel('wavenumber (rad/m)');
xlabel('wavelength (km)');
ylabel('S(k) (units?)');
grid on
lgd = legend(hl,{'k^{-5/3}','k_c (Ogive)'});



return 
function Iorg= compute_iorg(Scene_coord, Centroids_loc)

X = Centroids_loc;
nbins = 100000;
%[ny,nx]=size(Scene_coord);
% nx = length(Scene_coord.lon);
% ny = length(Scene_coord.lat);

Lx  = max(Scene_coord.lon) - min(Scene_coord.lon);
Ly = max(Scene_coord.lat) - min(Scene_coord.lat);

maxDist =  sqrt(Lx^2 + Ly^2);
bins = linspace(0, maxDist, nbins);

nnb_pair= find_nearest_neighbor_pairs(X);

figure(12);
h=histogram(nnb_pair.dist,bins,'Normalization','cdf');

cdf_dN = h.Values;
%binEdge = h.BinEdges;
binave = 0.5*(bins(1:end-1) + bins(2:end));

% cdf for Poisson process:
lam = size(X,1)/(Lx*Ly);
W = 1-exp(-lam *pi .* binave.^2);    % this is already CDF..

figure(12)
plot(binave, W);
hold on
plot(binave, cdf_dN);
hold off

Iorg = trapz(W,cdf_dN);

%figure()
%plot(W,cdf_dN);

return
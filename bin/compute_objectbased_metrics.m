function [M1, M2, center_cldinfo, hfig]= compute_objectbased_metrics(f, fthres,geo_coord, num_nan, circ_cen, highlight_flag)
% Purpose: this script will obtain the scene statistics of individual
% object properties from the regionprops matlab function;
%   mean object length;
%   maximum object length;
%   number of ojects (cloud number)
%   mean oject perimeter 
%
% Input: 
% . f: raw satellite field;
% . fthres: threshold value to generate a mask of the input cloud field.
% (threshold to identify the pixel as a cloud pixel or a non-cloud pixel)
% . geo_coord: geo_coord.lon, geo_coord.lat of the input satellite scene.
%
% Output:
%   

% use find_image_2d function;
% f input dimension: NLat x NLon (I think this should work)

%% to do: add optional input "highlight_flag" to highlight the cloud object within certain radius from the center of RHB.
switch nargin
    case 5
        highlight_flag = false;
    
end
[blob_coord,blob_image, fout, stats] = find_image_2d_updated(f', fthres);
% the 'blob_coord' here is not the real coordinate in degree. but it is the
% number of grid points;
% therefore, I need to convert blob_coord to the geographical coord. in
% degree. 
% ..
if length(stats.EqvDiam)>=1
    lon_indx = 1:length(geo_coord.lon);
    lat_indx = 1:length(geo_coord.lat);
    blob_GeoCoord.lon = interp1(lon_indx, geo_coord.lon, blob_coord(:,2));
    blob_GeoCoord.lat = interp1(lat_indx, geo_coord.lat, blob_coord(:,1));
    
    edgeblob_GeoCoord.lon = interp1(lon_indx, geo_coord.lon, stats.coord_EdgeClouds(:,2));
    edgeblob_GeoCoord.lat = interp1(lon_indx, geo_coord.lon, stats.coord_EdgeClouds(:,1));
    
    % plot results again:
    %hfig=0;
    hfig= figure(11);
    pcolor(geo_coord.lon, geo_coord.lat, f); shading flat;
    hb=colorbar;
    colormap(gray);
    hold on;
    plot(circ_cen(1), circ_cen(2), 'p','markersize',20, 'linewidth',2,'color','m');
    validx = ~isnan(geo_coord.lon) ;
    validy= ~isnan(geo_coord.lat);
    contour(geo_coord.lon(validx), geo_coord.lat(validy), fout(validx,validy)','c', 'linewidth',0.6);
    plot(blob_GeoCoord.lon, blob_GeoCoord.lat,'+r', 'linewidth',1.5,'markersize',8);
    %circle(circ_cen(1), circ_cen(2), 1, 'c',1.2);
    if highlight_flag
        
    end
    
    xlabel('longitude (^{\circ}E)')
    ylabel('latitude (^{\circ}N)');
    %set(get(hb,'xlabel'),'string','cloud mask');
    set(get(hb,'xlabel'),'string','radiance values');
    set(gca,'fontsize',12);
    hold off;
    %
    
    
    %% add a matric to identify whether the center is over a cloud pixel:
    center_cldinfo.flag = interp2(geo_coord.lon, geo_coord.lat, fout',circ_cen(1), circ_cen(2),'nearest');
    center_cldinfo.values = interp2(geo_coord.lon, geo_coord.lat, f,circ_cen(1), circ_cen(2),'nearest');
    
    
    
    %% start computing the basic object statistics: first group:
    
    % compute the area of the scene:
    Llon = max(geo_coord.lon) - min(geo_coord.lon);
    Llat = max(geo_coord.lat) - min(geo_coord.lat);
    dlon = mean(diff(geo_coord.lon),'omitnan');
    dlat = mean(diff(geo_coord.lat),'omitnan');
    
    Ascene = numel(f);    % units: pixel
    A = Ascene - num_nan;
    
    
    M1.MeanLen = mean(stats.EqvDiam,'omitnan');
    M1.MaxLen = max(stats.EqvDiam);
    M1.CloudNum = length(stats.EqvDiam);
    M1.MeanPerim = mean(stats.Perimeter,'omitnan');
    M1.CloudDensity = M1.CloudNum/(A-stats.Area_removed) * Ascene;         % units: #/pixel
    M1.CloudSize = mean(stats.CloudArea);
    M1.CloudCover = sum(stats.CloudArea)./(A-stats.Area_removed).*100;        % units: %   % removing the edge clouds or not doesn't affect cloud cover that much (assumed). 
    
    %% second group of indices:
    %% a. Simple Convective Aggregation Index (SCAI), Tobin et al. (2012) .
    %  (more like a "disaggregation index", the larger the index, the less aggregated the cloud field)
    % SCAI = No/Nmax * (D0/L) x 1000   (No = CloudNum)
    % Nmax = maximum number of pixels that can "exist" in the cloud scene;
    % D0: geometric mean of the distance of object pairs. confirmed with
    % the original paper.
    % L: characteristic length of the domain (it was taken as the length of the domain in Tobin et al. 2012)
    
    Nmax = (size(f,1) * size(f,2) - num_nan)/4;   % (use area of the entire scene / the area of the min. object size)
    No = M1.CloudNum;
    Np = No * (No - 1)/2;   %  number of non-repeated pairs of objects.
    % compute distance between centroid pairs:
    % this will need to use the blob_coordinate:
    dist_arry = zeros(Np,1);
    cnt =0;
    if No>1
        for i = 1:No-1
            for j = i+1:No
                cnt = cnt+1;
                dx = (blob_GeoCoord.lon(j) - blob_GeoCoord.lon(i))*111E3*cosd(blob_GeoCoord.lat(i));
                dy = (blob_GeoCoord.lat(j) - blob_GeoCoord.lat(i))*111E3;
                dist_arry(cnt) = sqrt( dx.^2 + dy.^2 );
            end
        end
        D0 = geomean(dist_arry,'omitnan');
        L = (max(geo_coord.lon) - min(geo_coord.lon))*111E3;   % characteristic length, doesn't need to be exact.
        M2.SCAI = No/Nmax * D0/L * 1000;   % per thousand.
    
   %% need some revision here:    
    else
        M2.SCAI = NaN;
%     else 
%         M2.SCAI = NaN;
    end
    
    %% b. Convective Organization Potential (COP) (White et al. 2018)
    cnt = 0;
    tmp = 0;
    for i =1:No-1
        for j = i+1:No
            cnt = cnt+1;
            tmp = tmp + (stats.EqvDiam(i) + stats.EqvDiam(j))/(sqrt(pi)*dist_arry(cnt));
        end
    end
    
    M2.COP = tmp/Np;
    
    %% c. maxRDF
    
    %% d. Degree Variance
    
    
    %% e. Organization Index (Iorg) & Iorg*
    % the normal Iorg (need to be tested).
    
    Centroids_loc=[blob_GeoCoord.lon, blob_GeoCoord.lat];
    
    if size(Centroids_loc,1)>2
        M2.Iorg = compute_iorg(geo_coord, Centroids_loc);
        
    %% need some revision here:    
%     elseif size(Centroids_loc,1)==1  % and the centroids is close to the center of the image.
%          dist2cen = sqrt((Centroids_loc(1) - circ_cen(1)).^2 -(Centroids_loc(2) - circ_cen(2)).^2);
%          if dist2cen<=0.5*(0.5*(max(geo_coord.lon)-min(geo_coord.lon))) && stats.EqvDiam>50
%              M2.Iorg = 1;
%          else
%              M2.Iorg =NaN;
%          end
    else
        M2.Iorg =NaN;
    end
    
    
    %% e2. inhibition Iorg:
    % this does not require removing the edge clouds. require info. of the
    % diameter of all the cloud object in the scene
    % I won't do periodic condition at this point;
    Centroids_loc=[blob_coord(:,2), blob_coord(:,1)];
    Centroids_loc2 = [stats.coord_EdgeClouds(:,2), stats.coord_EdgeClouds(:,1)];
    
    Centroids_all = [Centroids_loc; Centroids_loc2];
    
    R_all = [stats.EqvDiam, stats.EqvDiam_EdgeClouds]./2;
    
    if size(Centroids_loc,1)>2
        M2.Iorg_inhbited = compute_iorg_inhibition(size(f), Centroids_all, R_all);   %first arg: scene size [ NYx NX];
    else
        M2.Iorg_inhbited = NaN;
    end
    
    
    %% compute the clear sky metric in Janssens.
    f_cldmask = f>fthres;
    %f_cldmask(f_cldmask>fthres)=true;
    M2.Opensky_AreaFrac = compute_clear_sky_area(f_cldmask, A);
    
   
    
else
    M1_fn = fieldnames(M1);
    for iv = 1:length(M1_fn)
        M1.(M1_fn{iv}) = NaN;
    end
    
    M2.SCAI = NaN;
    M2.COP = NaN;
    M2.Iorg = NaN;
    M2.Iorg_inhibited = NaN;
    %M2.Opensky_AreaFrac = 
    center_cldinfo.flag=false;
    center_cldinfo.values = interp2(geo_coord.lon, geo_coord.lat, f,circ_cen(1), circ_cen(2),'nearest');
    hfig=NaN;
    
end

return
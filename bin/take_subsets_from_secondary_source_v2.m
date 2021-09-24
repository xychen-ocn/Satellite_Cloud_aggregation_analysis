function datasub = take_subsets_from_secondary_source_v2(NCFN, x0, y0, radius,channel )

% from Emultic2kmNC4_goes16 netcdf file;

% read in the netCDF file:
X2km = ncread(NCFN, 'X2km');
Y2km = ncread(NCFN, 'Y2km');

lat_proj_origin = ncreadatt(NCFN,'geos','latitude_of_projection_origin');
lon_proj_origin = ncreadatt(NCFN,'geos','longitude_of_projection_origin');

lat_pc = (Y2km-0)./111E3 + lat_proj_origin;
lon_pc = (X2km-0)./(111E3.*cosd(lat_pc)) + lon_proj_origin;


if strcmp(channel, 'VIS')
    %varn = 'refl_0_65um_nom';
    varn = 'VIS_006';
    
else
    %varn = 'temp_10_4um_nom';
    varn = 'IR_103';
    
end
dataval = ncread(NCFN, varn);    % nx, ny

%% 3. extract a subset of data points:
lon_subrange = [x0-radius, x0+radius];
lat_subrange = [y0-radius, y0+radius];

lon_mask = (lon_pc>=lon_subrange(1)) & (lon_pc <=lon_subrange(2));
lat_mask = (lat_pc>=lat_subrange(1)) & (lat_pc <=lat_subrange(2));


% use the index to extract the subset;
datasub.lon = lon_pc(lon_mask);
datasub.lat = lat_pc(lat_mask);
datasub.values = dataval(lon_mask, lat_mask)';

% 
% figure(1);clf;
% pcolor(datasub.lon, datasub.lat, datasub.values);
% hold on
% plot(x0, y0,'pm');
% shading flat;
% colormap(flipud(gray));
% colorbar;
% %caxis([250, 290]);
% 
% figure(2);
% histogram(datasub.values(datasub.values>290),'normalization','pdf');
% pd = fitdist(datasub.values(datasub.values>290), 'Normal')
% %pd = makedist('Normal','mu', 295 ,'sigma',0.4);
% x_values = [290:0.1:297];
% y = pdf(pd, x_values);
% hold on
% plot(x_values,y,'-r');
% 
% 
% figure(3)
% cldmask = datasub.values<292;
% pcolor(datasub.lon, datasub.lat, double(cldmask));
% shading flat
% colorbar;
% colormap(gray);
% 
% [blob_coord,blob_image, fout, stats] = find_image_2d_updated(double(cldmask), 0);
% 

%pause()


return
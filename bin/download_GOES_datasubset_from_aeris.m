function datasub = download_GOES_datasubset_from_aeris(SateDataDir, scene_center, scene_timevec, scene_radius)
%% purpose: download GOES data subset from the aeris remote server.
%  Input: center of the cloud scene, time_range_of_data, size of the cloud
%         scene.
%  Output: radiance product, TWP and CWP (if available.)
%

% urlbase ='/Volumes/Cumulonimbus/Research_backup/CU_postdoc_work/GOES16/';
% SateDataDir = [urlbase filesep 'observations.ipsl.fr/aeris/eurec4a-data/SATELLITES/GOES-E/2km_10min/2020'];
%default_source = '2km_10min';

%%%% ---- Preparation ---- %%%%
xcen = scene_center.lon;
ycen = scene_center.lat;
t0 = scene_timevec(1);


%%%% ---- Call another function to download data ---- %%%%
% try
%     %% -- combine the following script into one single function script, so that for the HALO circle, the data
%     %%% -- will only be read one time and the code will automatically
%     %%% switch the source netCDF depending on the local time of the
%     %%% snapshot of interest.
%     %%%
%     datasub=read_netCDF_subsets(scene_center.lon, scene_center.lat, scene_radius, t0, scene_timevec, SateDataDir);
%     datasub.source = '0.5km_01min';
% catch
    default_source = false;
    datasub = read_netCDF_subsets(xcen,ycen, scene_radius, t0, scene_timevec, SateDataDir, default_source);
    % datasub = read_netCDF_subsets(x0,y0, radius, t0, [t0, t0+1/24], SateDataDir, default_source);
    

    

    datasub.source = '2km_10min';   % this should be input.
    datasub.source ='2km_10min_fulldisk';
    
    
    
%end

%%%% ---- save the data: ----- %%%%



return
% This script is used to take out satellite data according to the location
% of the frontal segment identified on several ATOMIC days.

% try Jan 09:


clear all; close all;

addpath('/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_Cloud_aggregation_analysis/bin');

%% 1. prepare measurement location to get the subset of satellite cloud scene:
%   a. try the in-situ RHB measurement first.

RHB_locdir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/data';
%load([RHB_locdir filesep 'rhb_dailydata_grouped_by_SSTvariance.mat']);
load([RHB_locdir filesep 'rhb_daily_grouped_10min_data_0909latest.mat']);

% dates of interests:
groupn = 'strong_SSTvar';
DOI = rhbdates.(groupn);

% path to satellite data (on a hard drive currently)
urlbase ='/Volumes/Cumulonimbus/Research_backup/CU_postdoc_work/GOES16';  
%urlbase ='smb://10.0.0.150/GOES16';
%urlbase = 'https://';
SateDataDir = [urlbase filesep 'observations.ipsl.fr/aeris/eurec4a-data/SATELLITES/GOES-E'];

%% save data loc:
DataSvDirBase = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_Cloud_aggregation_analysis/test_RHB_GOES16_2km_10min';
% save this data locally:
DataSvDir = [DataSvDirBase filesep 'satellite_data_subsets/front_centered'];
if ~exist(DataSvDir,'dir')
    feval('mkdir', DataSvDir)
end

statsdir = [DataSvDir filesep 'stats_mat'];
if ~exist(statsdir,'dir')
    feval('mkdir', statsdir)
end



figsvdir = [DataSvDir filesep 'figs_out'];
if ~exist(figsvdir, 'dir')
     mkdir(figsvdir);
    % make subfolder in the parent folder:
     mkdir(figsvdir,'cloud_scenes');
     mkdir(figsvdir,'stats');
end

tres = 10/(24*60);
channel = 'VIS';              % or VIS.
sample_radius = 1;         %

for id = 1:length(DOI)       % loop through dates of interest.   % I know I don't have Jan09 and Jan18 for now, but I can obtain that through the server.
    
    DN = datestr(DOI(id),'mmmdd');
    load([RHB_locdir filesep DN filesep DN '_SST_front_segments.mat']);
    num_seg = length(seg);
    
    seg_stats = struct();

    for i = 5 %: num_seg
        disp(['-> Seg' num2str(i) ':']);
%         scene_center.lon = seg(i).front_loc(1);
%         scene_center.lat = seg(i).front_loc(2);
        scene_timevec_seg = seg(i).time;             % UTC time is used as input.
        time_window = 6;     % 3hr;
        scene_timevec = mean(scene_timevec_seg)-0.5*time_window/24:tres: mean(scene_timevec_seg)+0.5*time_window/24;

        scene_center.lon = repmat(seg(i).front_loc(1), 1, length(scene_timevec));
        scene_center.lat = repmat(seg(i).front_loc(2), 1, length(scene_timevec));

        %scene_timerange = [timevec(1), timevec(end)];
        
        
        %% 2.
        % select the size (radius) of the cloud scene of interest:
        scene_radius = 5;   % degree
        
        % call a function to "download" and subset satellite data:
        if strcmp(DN, 'Jan09') || strcmp(DN,'Jan18')
            data_source = '2km_10min_fulldisk';
        else
            data_source = '2km_10min';
        end
        disp(['data source:' data_source])
        
        try
        datasub{i} = download_GOES_datasubset_from_aeris_v2(SateDataDir,data_source, ...
        scene_center, scene_timevec, scene_radius, channel);
    
    % I should use either VIS or IR only. (do not blend them together.)
    % read IR and VIS separately, and then compute the statistics. 
    
    %%
     figsvdir_loc = [figsvdir filesep 'cloud_scenes' filesep DN];
        if ~exist(figsvdir_loc, 'dir')
            mkdir(figsvdir_loc);
        end
        
        % select a smaller region to compute the organization index,
        % Also, note whether there is a cloud over the observations.
        %for sample_radius = box_size
            %sample_radius = 0.25;  % degree (50km)  % try 0.25 and 0.1 degree (50km box and 20km box)
            sample_radius_km = round(sample_radius*111*cosd(mean(scene_center.lat)));
            
%             fn2 = ['object_metrics_r' num2str(sample_radius_km) 'km'];
%             fn1 = ['stats_moments_r' num2str(sample_radius_km) 'km'];
%             
            % then call subroutines to compute the three indices:
            % --a. call a function to compuate SCAI (I had this, need to build Iorg)
            
            
            fthres = 0.1;  % values larger than 10% of the value range. Let me express this more clearly.
            
            %for i = 1:num_seg
                Nt = size(datasub{i}.values,3);
                scene_center.lon = seg(i).front_loc(1);
                scene_center.lat = seg(i).front_loc(2);
                
                cldmask = construct_cloud_mask(datasub{i});
               % cldmask = construct_cloud_mask_following_BennerCurry98(datasub{i});
                fieldn = fieldnames(datasub{i});
                
                if strcmp(channel,'IR')
                    thres_val = max(datasub{i}.values(:)) -  fthres.*(max(datasub{i}.values(:)) - min(datasub{i}.values(:)));
                else
                    thres_val = min(datasub{i}.values(:)) +  fthres.*(max(datasub{i}.values(:)) - min(datasub{i}.values(:)));
                end
                
                
                thres_val = find_threshold_for_cloudy_pixels_following_BennerCurry98(datasub{i});
                % make movies out of both the raw data and the cloud mask
                % to demonstrate my point that it is possible that object
                % detection from cloud mask will need to be used together
                % with the data with reflectance or brightness temperature
                % values.
%                 for it =1:Nt
%                     % obtain identified clouds border to add to the plot
%                     % here.
%                     fin = datasub{i}.values(mask_lat,mask_lon,it);
%                     [blob_coord,blob_image, fout, stats] = find_image_2d_updated(fin, thres_val);
% 
%                     figure(1);
%                     subplot(2,2,[1,2])
%                     % raw data with colorbar range adjusted to reveal more
%                     % shallow clouds
%                     pcolor(datasub{i}.lon(mask_lon,it), datasub{i}.lat(mask_lat,it), fin);
%                     hold on
%                     plot(mean(datasub{i}.lon(:,it)), mean(datasub{i}.lat(:,it)),'pm','markersize',12);
%                     %plot(x0, y0,'pm');
%                     contour(geo_coord.lon(validx), geo_coord.lat(validy), fout(validx,validy)','c', 'linewidth',0.6);
%                     plot(blob_GeoCoord.lon, blob_GeoCoord.lat,'+r', 'linewidth',1.5,'markersize',8);
%                     %circle(circ_cen(1), circ_cen(2), 1, 'c',1.2);
%                     xlabel('longitude (^{\circ}E)')
%                     ylabel('latitude (^{\circ}N)');
%                     
%                     shading flat;
%                     colormap(gray);
%                     colorbar;
%                     caxis([0, 40]);
%                     axis('square');
% 
%                     
%                     subplot(2,2,[3,4])
%                     % cloud mask:
%                     pcolor(geo_coord.lon, geo_coord.lat, f); shading flat;
%                     hb=colorbar;
%                     colormap(gray);
%                     hold on;
%                     plot(circ_cen(1), circ_cen(2), 'p','markersize',20, 'linewidth',2,'color','m');
%                     validx = ~isnan(geo_coord.lon) ;
%                     validy= ~isnan(geo_coord.lat);
%                     contour(geo_coord.lon(validx), geo_coord.lat(validy), fout(validx,validy)','c', 'linewidth',0.6);
%                     plot(blob_GeoCoord.lon, blob_GeoCoord.lat,'+r', 'linewidth',1.5,'markersize',8);
%                     %circle(circ_cen(1), circ_cen(2), 1, 'c',1.2);
%                     xlabel('longitude (^{\circ}E)')
%                     ylabel('latitude (^{\circ}N)');
%                     set(get(hb,'xlabel'),'string','cloud mask');
%                     set(gca,'fontsize',12);
%                     hold off;
%                     
%                 end
                
                
                
                
                
                %% compute statistics
                clear M1 M2 cen_cldinfo M3
                %thres_val= thres_val*2;
                for it = 1:Nt
                    % construct mask to subsample the extracted satellite data:
                    mask_lon = abs(datasub{i}.lon(:,it) - scene_center.lon)<=sample_radius;
                    mask_lat = abs(datasub{i}.lat(:,it) - scene_center.lat)<=sample_radius;
                    
                    scene_cen=[scene_center.lon,scene_center.lat];
                    geo_coord.lon = datasub{i}.lon(mask_lon,it);
                    geo_coord.lat = datasub{i}.lat(mask_lat,it);
                    
                    num_nan = length(find(isnan(cldmask(mask_lat,mask_lon,it))==1));
                    
                    % vary this threshould for VIS data and for IR data:
                    %if strcmp(datasub.channel{it}, 'IR')
                    %    fthres=0.1;
                    %else
                    %end
                    
                   % [M1(it), M2(it), cen_cldinfo(it), hfig0]=compute_objectbased_metrics(cldmask(mask_lat,mask_lon,it), fthres, geo_coord, num_nan, scene_cen);
                    [M1(it), M2(it), cen_cldinfo(it), hfig0]=compute_objectbased_metrics(datasub{i}.values(mask_lat,mask_lon,it), thres_val, geo_coord, num_nan, scene_cen);

                    %% M3: compute spectral related indices:
                    % 1. the spectral length scale (as defined in Jonker et al. )
                    a=1;
                    for iv=1:length(fieldn)-1
                        vn = fieldn{iv};
                        if strcmp(vn, 'values')
                            data_s.(vn) = double(datasub{i}.(vn)(mask_lat,mask_lon,it));
                        elseif strcmp(vn,'channel')
                            data_s.(vn){1} = datasub{i}.(vn){it};
                        else
                            data_s.lon = double(datasub{i}.lon(mask_lon,it));
                            data_s.lat = double(datasub{i}.lat(mask_lat, it));
                        end
                        
                    end
                    
                    [M3(it).SpecLenScale, M3(it).OgiveLen] = compute_Jonker_spectral_length_scale(data_s, a, fthres);
                    
                    
                    figure(hfig0.Number); hold on
                    title({[datestr(datasub{i}.time(it)) ' (UTC)']; ['SCAI: ' num2str(M2(it).SCAI) '; Iorg: ' num2str(M2(it).Iorg) ...
                        '; Iorg_i: ' num2str(M2(it).Iorg_inhbited)]});
                    hold on
                    caxis([0 40])
                    plot(seg(i).lon, seg(i).lat,'--y','linewidth',1.5);
                    %pause
                    
                    % save figure:
                    % save the hfig0: 0 6 12 18
                    %pause
                    hr = hour(datetime(datestr(datasub{i}.time(it))));   % local time
                    %if mod(hr,6)==0 || hr==23
                    timehere = datestr(datasub{i}.time(it),'mmmddHHMM');
%                     figname = ['Seg' num2str(i,'%2.2i') '_Radius' num2str(sample_radius) 'dgrs_scene_over_RHB_' timehere '.jpg' ];
%                     xc_savefig(hfig0, figsvdir_loc, figname,[0 0 10 8]);
                    %end
                    frame(it) = getframe(hfig0);
                    
                    pause(0.5);
                    
                    close(hfig0)
                    
                    % save as movie:
                    
                    
                end
                

                % how we define the mask also makes a difference...(is
                % there a better definition of the mask??)
                
                
                % I can plot the statistics out for this segment, together with
                % time series? (Not sure what is the best approach.)
                seg_stats(i).M1 = M1;
                seg_stats(i).M2 = M2;
                seg_stats(i).M3 = M3;
                seg_stats(i).time_UTC = datasub{i}.time;  
                seg_stats(i).overhead_cldinfo = cen_cldinfo;
                
            %end
            
        end
            
        %end
        % note: highlight a cloud nearest to the center of front during
        % sampling.
        
        VideoName=[DN '_Seg' num2str(i,'%2.2i') '_Radius' num2str(sample_radius) 'dgrs_scene_over_RHB_channel' ...
             channel '_fthres'  num2str(thres_val, '%5.1f') '_6hour'  ];

        vdofn=[figsvdir_loc filesep VideoName '.mp4'];
        write_frames_into_video(frame, vdofn);
            
    
       % pause()
        
    end
    
    figure()
    subplot(2,1,1)
    plot(scene_timevec, [M3.OgiveLen],'linewidth',1.2);
    hold on
    plot(scene_timevec, [M3.SpecLenScale]./1E3,'linewidth',1.2);
    plot(scene_timevec, [M1.MeanPerim].*2,'linewidth',1.2);
    plot(scene_timevec, [M1.MeanLen].*2,'linewidth',1.2);

    legend('Ogive length','Spectral Length Scale','Mean Perimeter','Mean Diameter');
    datetick('x','hh:mm','keepticks')
    title([DN, '-seg', num2str(i,'%2.2i')]);
    ylabel('km');
    xlabel('time (UTC)');
    set(gca,'fontsize',14,'ytick',[0:10:140])
    
     subplot(2,1,2)
    yyaxis left
    plot(scene_timevec, [M2.SCAI],'linewidth',1.2);
    yyaxis right
    plot(scene_timevec, [M2.Iorg_inhbited],'linewidth',1.2);
    ylim([0 1])
%     plot(scene_timevec, [M1.MeanPerim].*2,'linewidth',1.2);
%     plot(scene_timevec, [M1.MeanLen].*2,'linewidth',1.2);

    legend('SCAI','Iorg_i');
    datetick('x','hh:mm','keepticks')
    title([DN, '-seg', num2str(i,'%2.2i')]);
    xlabel('time (UTC)');
    set(gca,'fontsize',14)

    
    % add code to advance data to the next day.
    
    % save data back into structure.
    save([statsdir filesep DN '_' channel '_cloud_statistics_over_segments_boxsize' num2str(sample_radius) ...
        '_' channel '_fthres'  num2str(thres_val, '%5.1f') '_seg05-06.mat'],'seg_stats');
    
    dataFN = [DN '_' channel '_local_time_CloudScene_Radius' num2str(scene_radius) 'dgrs_TimeWindow_' ...
             num2str(time_window) 'hrs_Centered_on_frontal_segments.mat'];
    save([DataSvDir filesep dataFN],'datasub');
    
end


figsvdir = [DataSvDir filesep 'figs_out'];
if ~exist(figsvdir, 'dir')
     mkdir(figsvdir);
    % make subfolder in the parent folder:
     mkdir(figsvdir,'cloud_scenes');
     mkdir(figsvdir,'stats');
end


% compute some cloud statistics out of it:
box_size=[2.5];  %[0.5, 0.25, 0.1, 0.05];
for id = 1:length(DOI)
    
        % load data:
        DN = datestr(DOI(id),'mmmdd');
        
        dataFN =[DN 'local_time_CloudScene_Radius' num2str(scene_radius) 'dgrs_TimeWindow_' ...
             num2str(time_window) 'hrs_Centered_on_frontal_segments.mat'];
        load([DataSvDir filesep dataFN],'datasub');
        
        figsvdir_loc = [figsvdir filesep 'cloud_scenes' filesep DN];
        if ~exist(figsvdir_loc, 'dir')
            mkdir(figsvdir_loc);
        end
        
        % select a smaller region to compute the organization index,
        % Also, note whether there is a cloud over the observations.
        for sample_radius = box_size
            %sample_radius = 0.25;  % degree (50km)  % try 0.25 and 0.1 degree (50km box and 20km box)
            sample_radius_km = round(sample_radius*111*cosd(mean(scene_center.lat)));
            
%             fn2 = ['object_metrics_r' num2str(sample_radius_km) 'km'];
%             fn1 = ['stats_moments_r' num2str(sample_radius_km) 'km'];
%             
            % then call subroutines to compute the three indices:
            % --a. call a function to compuate SCAI (I had this, need to build Iorg)
            
            
            fthres = 0.1;  % values larger than 10% of the value range. Let me express this more clearly.
            
            for i = 2:num_seg
                Nt = size(datasub{i}.values,3);
                scene_center.lon = seg(i).front_loc(1);
                scene_center.lat = seg(i).front_loc(2);
                
                cldmask = construct_cloud_mask(datasub{i});
                fieldn = fieldnames(datasub{i});
                
                clear M1 M2 cen_cldinfo M3
                for it = 1:15; %Nt
                    % construct mask to subsample the extracted satellite data:
                    mask_lon = abs(datasub{i}.lon(:,it) - scene_center.lon)<=sample_radius;
                    mask_lat = abs(datasub{i}.lat(:,it) - scene_center.lat)<=sample_radius;
                    
                    scene_cen=[scene_center.lon,scene_center.lat];
                    geo_coord.lon = datasub{i}.lon(mask_lon,it);
                    geo_coord.lat = datasub{i}.lat(mask_lat,it);
                    
                    num_nan = length(find(isnan(cldmask(mask_lat,mask_lon,it))==1));
                    
                    % vary this threshould for VIS data and for IR data:
                    %if strcmp(datasub.channel{it}, 'IR')
                    %    fthres=0.1;
                    %else
                    %end
                    
                    [M1(it), M2(it), cen_cldinfo(it), hfig0]=compute_objectbased_metrics(cldmask(mask_lat,mask_lon,it), fthres, geo_coord, num_nan, scene_cen);
                    
                    %% M3: compute spectral related indices:
                    % 1. the spectral length scale (as defined in Jonker et al. )
                    a=1;
                    for iv=1:length(fieldn)-1
                        vn = fieldn{iv};
                        if strcmp(vn, 'values')
                            data_s.(vn) = double(datasub{i}.(vn)(mask_lat,mask_lon,it));
                        elseif strcmp(vn,'channel')
                            data_s.(vn){1} = datasub{i}.(vn){it};
                        else
                            data_s.lon = double(datasub{i}.lon(mask_lon,it));
                            data_s.lat = double(datasub{i}.lat(mask_lat, it));
                        end
                        
                    end
                    
                    [M3(it).SpecLenScale, M3(it).OgiveLen] = compute_Jonker_spectral_length_scale(data_s, a, fthres);
                    
                    
                    figure(hfig0.Number); hold on
                    title({[datestr(datasub{i}.time(it)) ' (UTC)']; ['SCAI: ' num2str(M2(it).SCAI) '; Iorg: ' num2str(M2(it).Iorg) ...
                        '; Iorg_i: ' num2str(M2(it).Iorg_inhbited)]});
                    hold on
                    plot(seg(i).lon, seg(i).lat,'--y','linewidth',1.5);
                    
                    
                    % save figure:
                    % save the hfig0: 0 6 12 18
                    %pause
                    hr = hour(datetime(datestr(datasub{i}.time(it))));   % local time
                    %if mod(hr,6)==0 || hr==23
                    timehere = datestr(datasub{i}.time(it),'mmmddHHMM');
                    figname = ['Seg' num2str(i,'%2.2i') '_Radius' num2str(sample_radius) 'dgrs_scene_over_RHB_' timehere '.jpg' ];
                    xc_savefig(hfig0, figsvdir_loc, figname,[0 0 10 8]);
                    %end
                    close(hfig0)
                end
                
                
                % I can plot the statistics out for this segment, together with
                % time series? (Not sure what is the best approach.)
                seg_stats(i).M1 = M1;
                seg_stats(i).M2 = M2;
                seg_stats(i).M3 = M3;
                seg_stats(i).time_UTC = datasub{i}.time;  
                seg_stats(i).overhead_cldinfo = cen_cldinfo;
                
            end
            
            % save data back into structure.
            save([statsdir filesep 'stats' filesep DN '_cloud_statistics_over_segments.mat'],'seg_stats');
            
        end
        
    end
    

% note: I need to define the same threshold for all the fields so that the
% cloud object selected doesn't change much. 
% Is there a cloud mask variable that I can use directly from the satellite
% images. (check Judith Curry's Iorg inhibition paper again, they describe
% the selection of threshold in details, I believe.)
%

% tailor to these images of small fractions of distinct clouds over open
% ocean. 

% Benner and Curry: Small Tropical Cumulus Clouds. 
% visible radiance has to be higher than the visible threshold
% IR radiance has to be lower than the IR threshold. (brightness
% temperature); Due to the fact that a cloud should be more reflective and
% colder than the underlying ocean surface. 
%

% Visible threshold was determined from the histogram of all visible
% radiance values in the image, and similarly also for the IR threshold
% (Figure 1 in Benner and Curry 1998).
VIS_thres = ;
IR_thres = ;


% Also, it might be useful to select the cloud object cloesest to the front
% at the time when RHB was at the beginning of the front segment. 
% and see how it evolves downstream in time. (produce a video from the
% cloud scene.)
% 
% no cloud mask for the Emultic2kmNC4 data. 
% 

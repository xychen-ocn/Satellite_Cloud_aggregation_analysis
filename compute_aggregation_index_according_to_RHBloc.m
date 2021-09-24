% This script will compute the aggregation index in terms of SCAI, Iorg,
% and TWP mesoscale variance with the cloud scene centered at the RHB SST
% and flux measurements.
clear all; close all;

addpath('/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_Cloud_aggregation_analysis/bin');

%% 1. prepare measurement location to get the subset of satellite cloud scene:
%   a. try the in-situ RHB measurement first.
RHB_locdir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb';
%load([RHB_locdir filesep 'rhb_dailydata_grouped_by_SSTvariance.mat']);
load([RHB_locdir filesep 'rhb_dailydata_grouped_by_SSTvariance_added_aggregation_metrics_4boxsizes.mat']);

% dates of interests:
% groupn = 'strong_SSTvar';
% DOI = rhb_days.(groupn).datenumbers;

DOI = [rhb_days.strong_SSTvar.datenumbers rhb_days.moderate_SSTvar.datenumbers];


% path to satellite data (on a hard drive currently)
urlbase ='/Volumes/Cumulonimbus/Research_backup/CU_postdoc_work/GOES16';
SateDataDir = [urlbase filesep 'observations.ipsl.fr/aeris/eurec4a-data/SATELLITES/GOES-E'];

%% save data loc:
DataSvDirBase = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_Cloud_aggregation_analysis/test_RHB_GOES16_2km_10min';
% save this data locally:
DataSvDir = [DataSvDirBase filesep 'satellite_data_subsets'];
if ~exist(DataSvDir,'dir')
    feval('mkdir', DataSvDir)
end


figsvdir = [DataSvDirBase filesep 'figs_out'];
if ~exist(figsvdir, 'dir')
     mkdir(figsvdir);
    % make subfolder in the parent folder:
     mkdir(figsvdir,'cloud_scenes');
     mkdir(figsvdir,'stats');
end

%% start the core part of the program:
for id = 3:length(DOI)       % loop through dates of interest.   % I know I don't have Jan09 and Jan18 for now, but I can obtain that through the server.
    
    DN = datestr(DOI(id),'mmmdd');
    scene_center.lon = rhb_days.(groupn).(DN).lon;
    scene_center.lat = rhb_days.(groupn).(DN).lat;
    scene_timevec = rhb_days.(groupn).(DN).time;      % UTC time is used as input.
    %scene_timerange = [timevec(1), timevec(end)];
    

    %% 2.
    % select the size (radius) of the cloud scene of interest:
    scene_radius = 2;   % degree
    
    % call a function to "download" and subset satellite data:

    datasub = download_GOES_datasubset_from_aeris(SateDataDir,scene_center, scene_timevec, scene_radius);
    
    dataFN = [DN 'local_time_CloudScene_Radius' num2str(scene_radius) 'dgrs_Centered_RHB.mat'];
    save([DataSvDir filesep dataFN],'datasub');
    
end



box_size=[0.1,0.05];  %[0.5, 0.25, 0.1, 0.05];
for id = 1:length(DOI)
   % if DOI(id) >=datenum(2020,1,22)
    if DOI(id) ==datenum(2020,2,3)  || DOI(id) ==datenum(2020,1,28)  
        if id> 8
            groupn = 'moderate_SSTvar';
        else
            groupn = 'strong_SSTvar';
        end
        % load data:
        DN = datestr(DOI(id),'mmmdd');
        scene_center.lon = rhb_days.(groupn).(DN).lon;
        scene_center.lat = rhb_days.(groupn).(DN).lat;
        scene_timevec = rhb_days.(groupn).(DN).time;      % UTC time is used as input.
        scene_radius = 2;   % degree
        dataFN =[DN 'local_time_CloudScene_Radius' num2str(scene_radius) 'dgrs_Centered_RHB.mat'];
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
            
            fn2 = ['object_metrics_r' num2str(sample_radius_km) 'km'];
            fn1 = ['stats_moments_r' num2str(sample_radius_km) 'km'];
            
            % then call subroutines to compute the three indices:
            % --a. call a function to compuate SCAI (I had this, need to build Iorg)
            
            
            fthres = 0.1;  % values larger than 10% of the value range. Let me express this more clearly.
            
            Nt = size(datasub.values,3);
            M1{Nt} = [];    % first group: basic object statistics (statistical moments)
            M2{Nt} = [];    % second group: aggregation indics (object metrics )
            
            cldmask = construct_cloud_mask(datasub);
            for it = 1:Nt
                % construct mask to subsample the extracted satellite data:
                mask_lon = abs(datasub.lon(:,it) - scene_center.lon(it))<=sample_radius;
                mask_lat = abs(datasub.lat(:,it) - scene_center.lat(it))<=sample_radius;
                
                scene_cen=[scene_center.lon(it),scene_center.lat(it)];
                geo_coord.lon = datasub.lon(mask_lon,it);
                geo_coord.lat = datasub.lat(mask_lat,it);
                
                num_nan = length(find(isnan(cldmask(mask_lat,mask_lon,it))==1));
                
                % vary this threshould for VIS data and for IR data:
                %if strcmp(datasub.channel{it}, 'IR')
                %    fthres=0.1;
                %else
                %end
                
                [M1{it}, M2{it}, cen_cldinfo{it}, hfig0]=compute_objectbased_metrics(cldmask(mask_lat,mask_lon,it), fthres, geo_coord, num_nan, scene_cen);
                
                if isstruct(M2{it})
                    figure(hfig0.Number); hold on
                    title({[datestr(datasub.time(it)-4/24) ' (local time)']; ['SCAI: ' num2str(M2{it}.SCAI) '; Iorg: ' num2str(M2{it}.Iorg)]});
                    
                    % save figure:
                    % save the hfig0: 0 6 12 18
                    hr = hour(datetime(datestr(datasub.time(it)-4/24)));   % local time
                    if mod(hr,6)==0 || hr==23
                        timehere = datestr(datasub.time(it)-4/24,'mmmddHHMM');
                        figname = ['Radius' num2str(sample_radius) 'dgrs_scene_over_RHB_' timehere '.jpg' ];
                        xc_savefig(hfig0, figsvdir_loc, figname,[0 0 10 8]);
                    end
                    close(hfig0)
                    
                    M1_variable_names = fieldnames(M1{it});
                    M2_variable_names = fieldnames(M2{it});
                    
                    for iv = 1:length(M1_variable_names)
                        VarN = M1_variable_names{iv};
                        rhb_days.(groupn).(DN).cloud_aggregation.(fn1).(VarN)(it) = M1{it}.(VarN);
                    end
                    
                    for iv = 1:length(M2_variable_names)
                        VarN = M2_variable_names{iv};
                        rhb_days.(groupn).(DN).cloud_aggregation.(fn2).(VarN)(it) = M2{it}.(VarN);
                    end
                      
                else
                    for iv = 1:length(M1_variable_names)
                        VarN = M1_variable_names{iv};
                        rhb_days.(groupn).(DN).cloud_aggregation.(fn1).(VarN)(it) = NaN;
                    end
                    
                    for iv = 1:length(M2_variable_names)
                        VarN = M2_variable_names{iv};
                        rhb_days.(groupn).(DN).cloud_aggregation.(fn2).(VarN)(it) = NaN;
                    end
                    
                end
                    rhb_days.(groupn).(DN).cld_flag(it) = cen_cldinfo{it}.flag;
                    rhb_days.(groupn).(DN).cld_relative_refval(it) = cen_cldinfo{it}.values;
            end
            
            % save data back into structure.
            
        end
        
    end
    
end


%     for i =1:Nt
%         Iorg(i) = M2{i}.Iorg;
%         SCAI(i) = M2{i}.SCAI;
%     end
%     figure()
%     plot(scene_timevec, Iorg);
%     hold on
%     plot(scene_timevec, SCAI./max(SCAI));
%     
    
    % --c. call a function to compute TWP variance
    %TWPvar_ratio = compute_BB17_TWP_variance_ratio(datasub.TWP);   %this
    %type of calculation is not good enough
    
    
    % save data back into structure.
%     rhb_days.(groupn).(DN).cloud_aggregation.object_metrics = M2;
%     rhb_days.(groupn).(DN).cloud_aggregation.stats_moments = M1;
    %rhb_days.weak_SSTvar.(DN).cloud_aggregation.TWPvar_ratio = TWPvar_ratio;
    
%end

% save data:
dataname = '';
save([RHB_locdir filesep 'rhb_dailydata_grouped_by_SSTvariance_added_aggregation_metrics_4boxsizes.mat'], 'rhb_days');

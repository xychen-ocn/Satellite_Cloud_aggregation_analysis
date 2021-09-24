% ----------------------------------------------------------------------- %
% Purpose: this script is used to calculate different indices that
%          objectively quantifies the degree of self aggregation.
% Data: GOES 0.5km_1min ABI regridded product by Hauke in a selected area.
%       variable is outgoing radiance. (two different channels)
% 
% Functions: /bin/find_image_2d.m, this function uses the matlab "regionprops" 
%           function to identify objects.  
%
% Steps: 1) read in the satellite data from OpenDap;
%        2) call functions to compute object-based indices for aggregation.
%        3) save results by dates and/or time.
% 
% ----------------------------------------------------------------------- %

clear all; clc; close all;

% 1) construct data path:
AERIS_ROOT = 'https://observations.ipsl.fr/thredds/dodsC/EUREC4A/SATELLITES';
SateName = 'GOES-E';
ProdName = '0.5km_01min';
Year = '2020';

SateDataDir = [AERIS_ROOT filesep SateName];% filesep ProdName filesep Year];

% JOANNE circle mean product:
JOANNE_Dir ='/Users/cirrus/Documents/MATLAB/shallow_convection/Obsv/EUREC4A/Dropsondes';


%load([JOANNE_Dir filesep 'JNLv4_HALO.mat']);

% start day and end day of dropsonde circles (start with P3 and ends with
% HALO)
% day0 = datenum(2020,1,17,0,0,0);
% dayN = datenum(2020,2,15,0,0,0);
% dates_arry = day0:dayN;
% actually: 
% dates_arry = unique(floor(JNLv4.circle_time));
% num_of_days = length(dates_arry);

% format1 = 'yyyy_mm_dd';
% format2 = 'yyyymmdd';



% % information obtained from web search
% local_sunrise_hr = 7; 
% local_sunset_hr = 17; 
% % barbados in Atlantic Standard Time (GMT-4)
% UTCoffset = -4;

addpath('./bin');

%% data needed to take the subset in space and time;
radius = 5;  % units: degree
Twindow = 6;
tinc = 10;      % units: 1 instances = 1 min (or 1 unit of temporal resolution.)
tres = 1;      % units: minute;

testdirn = ['test_v0_dt' num2str(tinc*tres,'%2.2i') 'min'];
if ~exist(testdirn, 'dir')
    mkdir(testdirn);
end
figsvdir = [testdirn filesep 'figs_out'];
if ~exist(figsvdir, 'dir')
    mkdir(figsvdir);
    % make subfolder in the parent folder:
     mkdir(figsvdir,'cloud_scenes');
     mkdir(figsvdir,'stats');
end

datasvdir = [testdirn filesep 'satellite_data_subsets'];
if ~exist(datasvdir, 'dir')
    mkdir(datasvdir);
end




% if strcmp(platform, 'P3')
%     Twindow = 6;    % units: hour;
% else % HALO
%     Twindow = 12;
% end
%% input is P3:
platform = 'P3';
load([JOANNE_Dir filesep 'JNLv4_' platform '.mat']);
num_of_P3circs = length(JNLv4_P3.circle_time);

for ic = 8 %1:num_of_P3circs
    
    x0 = JNLv4_P3.circle_lon(ic);
    y0 = JNLv4_P3.circle_lat(ic);
    t0 = JNLv4_P3.circle_time(ic);
    
    % display circle segment id and mean time:
    circID = JNLv4_P3.segment_id{ic};
    disp([circID ':' datestr(t0)]);
    
   % time_range = (t0-Twindow/2/24) :tinc*tres/(24*60): (t0+Twindow/2/24);
    time_range = t0;
    
    % skip Feb 05 for now: (revise this in a better way later)
    if floor(t0) ~= datenum(2020,2,5,0,0,0)
        %% -- combine the following script into one single function script, so that for the HALO circle, the data
        %%% -- will only be read one time and the code will automatically
        %%% switch the source netCDF depending on the local time of the
        %%% snapshot of interest.
        %%%
        datasub=read_netCDF_subsets(x0,y0, radius, t0, time_range, SateDataDir);
        datasub.source = ProdName;
    else
        default_source = false;
        datasub = read_netCDF_subsets(x0,y0, radius, t0, time_range, SateDataDir, default_source);
       % datasub = read_netCDF_subsets(x0,y0, radius, t0, [t0, t0+1/24], SateDataDir, default_source);

        datasub.source = '2km_10min';
    end
        
        
        %% calculate aggregation indices with an object-based algorithm from Dias et al. (2012)
        %  refer to Janssens et al. (2021) paper, there are X object based
        %  metrics that we can use to quantify the degree of aggregation from
        %  the OLR.
        Nt = size(datasub.values,3);
        M1{Nt} = [];
        M2{Nt} = [];
        
        for it = 1:Nt
            channel = datasub.channel{it};
            range_data = max(max(datasub.values(:,:,it))) - min(min(datasub.values(:,:,it)));
            %% note: how to deal with NaN...
            % calculate area that is occupied by NaN and remove those area when
            % calculating the statistics.
            
            % scale the data values to be within range 0~1 (relative magnitude). the scaling is
            % different for VIS and IR input data;
            if strcmpi(channel, 'VIS')
                %datain = satedata(:,:,it)./max(max(satedata(:,:,it)));
                % cloudy pixels associates with larger reflectance. increasing
                % value --> clouds
                f = (datasub.values(:,:,it)-min(min(datasub.values(:,:,it))))./range_data;
                % scale data by its range:
            else
                % brightness temperature: lower --> cloudy
                % after normalization: 0=no cloud, 1, cloudy.
                f = (max(max(datasub.values(:,:,it)))-datasub.values(:,:,it))./range_data;
            end
            
            fthres = 0.1;  % values larger than 10% of the value range. Let me express this more clearly.
            geo_coord.lon = datasub.lon;
            geo_coord.lat = datasub.lat;
            
            circ_cen=[x0,y0];
            num_nan = length(find(isnan(f)==1));
            %f(isnan(f))=0;
            %ftemp=f(1:end-1, 1:end-1);
            
            [M1{it}, M2{it}, hfig0]=compute_objectbased_metrics(f, fthres, geo_coord, num_nan, circ_cen);
            
            % save figure:
            % save the hfig0:
            if it==1 || it==Nt || it==round(Nt/2)
                timehere = datestr(datasub.time(it),'mmddHHMM');
                figname = [JNLv4_P3.segment_id{ic} '_scene_' timehere '.jpg' ];
                xc_savefig(hfig0, [figsvdir filesep 'cloud_scenes'], figname,[0 0 10 8]);
            end
            
        end
        
        M1_names =fieldnames(M1{1});
        MG1{ic} = struct(M1_names{1}, 0);
        for iv = 1:length(M1_names)
            for it =1:Nt
                MG1{ic}.(M1_names{iv})(it)=M1{it}.(M1_names{iv});
                
            end
        end
        
        for it =1:Nt
            
            SCAI{ic}(it) = M2{it}.SCAI;
        end
        
        % make some plots:
        hfig=figure(1);clf;
        for iv = 1:length(M1_names)
            varn = M1_names{iv};
            subplot(3,2,iv);
            plot(datasub.time, MG1{ic}.(varn)(1:Nt),'linewidth',1.2);
            title(varn);
            datetick('x');
            hold on;
            ylimval = get(gca,'ylim');
            plot([t0, t0],ylimval,'--r','linewidth',1.2);
            set(gca,'fontsize',12);
        end
        
        subplot(3,2,5);
        plot(datasub.time, SCAI{ic}(1:Nt),'linewidth',1.2);
        title('SCAI');
        datetick('x')
        hold on;
        ylimval = get(gca,'ylim');
        plot([t0, t0],ylimval,'--r','linewidth',1.2);
        set(gca,'fontsize',12);
        
        %%
        %     figsvdir='fig_test_v0';
        figname = [circID '_stats.jpg'];
        %pause
        xc_savefig(hfig,[figsvdir filesep 'stats'], figname, [0 0 10 8]);
        
        
        % save datasub into a matlab structure to avoid reading in the cloud scene
        % again;
        % matfilename: platform_segmentID_radius_dt.mat'
        matfilename = ['P3_' circID '_SceneRadius_' num2str(radius) ...
            '_dt' num2str(tinc*tres,'%2.2i') 'min.mat'];
        save([datasvdir filesep matfilename],'datasub');
        
%     else
%         SCAI{ic} = NaN;
%         MG1{ic} = NaN;
%     end
end  %[4,8]
%% save the data
save([testdirn filesep 'P3_test_objective_aggregation_index_dt' num2str(tinc*tres,'%2.2i') 'min_added0205.mat'],'MG1','SCAI');



%% input is HALO:
clear MG1 SCAI
platform = 'HALO';
load([JOANNE_Dir filesep 'JNLv4_' platform '.mat']);

HALO_dates = fieldnames(JNLv4_HALO_circgrp);
num_of_HLdays = length(HALO_dates);

for ic =[7, 12] %1:num_of_HLdays   % number of days:
    
    dn = HALO_dates{ic};
    
    x0 = mean(JNLv4_HALO_circgrp.(dn).circle_lon,'omitnan');
    y0 = mean(JNLv4_HALO_circgrp.(dn).circle_lat,'omitnan');
    t0_c1 = JNLv4_HALO_circgrp.(dn).circle_time(1);
    t0_cN = JNLv4_HALO_circgrp.(dn).circle_time(end);
    
    % display circle segment id and mean time:
    circID = JNLv4_HALO_circgrp.(dn).segment_id{1};
    disp(['@' circID ':' datestr(t0_c1) ' ~ ' datestr(t0_cN)]);
    
    if floor(t0_c1-Twindow/48)<floor(t0_c1)
        tt0 = t0_c1 - 0.5/24;
    else
        tt0 = t0_c1-Twindow/48;
    end
    
    if floor(t0_cN+Twindow/48)>floor(t0_cN)
        ttN = t0_cN+0.5/24;
    else
        ttN = t0_cN+Twindow/48;
    end
    
    time_range = tt0 :tinc*tres/(24*60): ttN;
    
    % skip Feb 05 for now: (revise this in a better way later)
    if floor(t0_c1) ~= datenum(2020,2,5) && floor(t0_c1) ~= datenum(2020,2,15) && floor(t0_c1)~=datenum(2020,1,24)
        
        %% -- combine the following script into one single function script, so that for the HALO circle, the data
        %%% -- will only be read one time and the code will automatically
        %%% switch the source netCDF depending on the local time of the
        %%% snapshot of interest.
        %%%
        datasub = read_netCDF_subsets(x0,y0, radius, t0_c1, time_range, SateDataDir);
    else
        default_source = false;
        datasub = read_netCDF_subsets(x0,y0, radius, t0_c1, time_range, SateDataDir, default_source);
    end
       
%        matfilename = ['HALO_' circID '_SceneRadius_' num2str(radius) ...
%            '_dt' num2str(tinc*tres,'%2.2i') 'min.mat'];
%        load([datasvdir filesep matfilename],'datasub');
%        
        
        
        %% calculate aggregation indices with an object-based algorithm from Dias et al. (2012)
        %  refer to Janssens et al. (2021) paper, there are X object based
        %  metrics that we can use to quantify the degree of aggregation from
        %  the OLR.
        Nt = size(datasub.values,3);
        M1{Nt} = [];
        M2{Nt} = [];
        
        for it = 1:Nt
            channel = datasub.channel{it};
            range_data = max(max(datasub.values(:,:,it))) - min(min(datasub.values(:,:,it)));
            %% note: how to deal with NaN...
            % calculate area that is occupied by NaN and remove those area when
            % calculating the statistics.
            
            % scale the data values to be within range 0~1 (relative magnitude). the scaling is
            % different for VIS and IR input data;
            if strcmpi(channel, 'VIS')
                %datain = satedata(:,:,it)./max(max(satedata(:,:,it)));
                % cloudy pixels associates with larger reflectance. increasing
                % value --> clouds
                f = (datasub.values(:,:,it)-min(min(datasub.values(:,:,it))))./range_data;
                % scale data by its range:
            else
                % brightness temperature: lower --> cloudy
                % after normalization: 0=no cloud, 1, cloudy.
                f = (max(max(datasub.values(:,:,it)))-datasub.values(:,:,it))./range_data;
            end
            
            fthres = 0.1;  % values larger than 10% of the value range. Let me express this more clearly.
            geo_coord.lon = datasub.lon;
            geo_coord.lat = datasub.lat;
            
            circ_cen=[x0,y0];
            num_nan = length(find(isnan(f)==1));
            %f(isnan(f))=0;
            %ftemp=f(1:end-1, 1:end-1);
            
            [M1{it}, M2{it}, hfig0]=compute_objectbased_metrics(f, fthres, geo_coord, num_nan, circ_cen);
            % save the hfig0:
            if 1==1
            if it==1 || it==Nt || it==round(Nt/2)
                timehere = datestr(datasub.time(it),'mmddHHMM');
                figname = [JNLv4_HALO_circgrp.(dn).segment_id{1} '_scene_' timehere '.jpg' ];
                xc_savefig(hfig0,[figsvdir filesep 'cloud_scenes'], figname,[0 0 10 8]);
            end
            end
        end
        
        M1_names =fieldnames(M1{1});
        MG1{ic} = struct(M1_names{iv},0); 
        for iv = 1:length(M1_names)
            for it =1:Nt
                 
                MG1{ic}.(M1_names{iv})(it)=M1{it}.(M1_names{iv});
                
            end
        end
        
        for it =1:Nt
            
            SCAI{ic}(it) = M2{it}.SCAI;
        end
        
        
        if 1==1
        % make some plots:
        hfig=figure(1);clf;
        for iv = 1:length(M1_names)
            varn = M1_names{iv};
            subplot(3,2,iv);
            plot(datasub.time, MG1{ic}.(varn),'linewidth',1.2);
            title(varn);
            datetick('x');
            hold on;
            ylimval = get(gca,'ylim');
            plot([t0_c1, t0_c1],ylimval,'--r','linewidth',1.2);
            plot([t0_cN, t0_cN],ylimval,'--r','linewidth',1.2);
            % if t0_c1 < min(datasub.time)
            xlim([min(t0_c1, datasub.time(1)), max(datasub.time(end),t0_cN)]);
            % end
            set(gca,'fontsize',12);
        end
        
        subplot(3,2,5);
        plot(datasub.time, SCAI{ic},'linewidth',1.2);
        title('SCAI');
        datetick('x')
        hold on;
        ylimval = get(gca,'ylim');
        % plot([t0, t0],ylimval,'--r','linewidth',1.2);
        plot([t0_c1, t0_c1],ylimval,'--r','linewidth',1.2);
        plot([t0_cN, t0_cN],ylimval,'--r','linewidth',1.2);
        %if t0_c1 < min(datasub.time)
        xlim([min(t0_c1, datasub.time(1)), max(datasub.time(end),t0_cN)]);
        % end
        set(gca,'fontsize',12);
        
        %%
        % figsvdir='fig_test_v0';
        figname = [JNLv4_HALO_circgrp.(dn).segment_id{1} '_stats.jpg'];
        
        xc_savefig(hfig,[figsvdir filesep 'stats'], figname, [0 0 10 8]);
        end
        
%     else
%         MG1{ic}=NaN;
%         SCAI{ic} = NaN;
%     end
    
    %% should have saved datasub for each day to avoid re-reading data;
    %
    matfilename = ['HALO_' circID '_SceneRadius_' num2str(radius) ...
        '_dt' num2str(tinc*tres,'%2.2i') 'min.mat'];
    save([datasvdir filesep matfilename],'datasub');
    
end
%% save the data
save([testdirn filesep 'HALO_test_objective_aggregation_index_dt' num2str(tinc*tres,'%2.2i') 'min_complete.mat'],'MG1','SCAI', 'HALO_dates');
% note: HALO days: Feb05, Feb15 needs another data source located in
% "2km_10min" folder;
% P3 day: Feb 05 idem.




%% Obsolete:
if 1==0
 % determine which product to use based on the circle mean time:
    if (t0 + UTCoffset/24) <= datehere + local_sunrise_hr/24 ...
            || (t0+ UTCoffset/24)>=datehere +local_sunset_hr/24
        channel = 'IR';
    else
        channel = 'VIS';
    end
    
    % channel here should be selected based on the day (we know which are
    % the nighttime circles.)  --> (revised). 
    if strcmpi(channel, 'VIS')
        channelID='02';
        varn='C02';
    else
        channelID='13';
        varn='C13';
    end
        
    NCFN= [SateDataDir filesep datestr(datehere,format1) filesep ...
        'GOES_' channelID '_8N-18N-62W-50W_' datestr(datehere,format2) '.nc'];
    
    
    % carve out a specific area from the dataset to get the circle location for this day: (if there are more than 1
    % circle for the day, then use the averaged center of the circle.)
    
    % read data from the netCDF file above:  (the function here probably
    % needs to be combined with teh function below).
    [data, errmsg] = read_ncfile_to_struct(NCFN);
    for iv = 1:length(errmsg)
        if ~isempty(errmsg{iv})
            disp(['warning: ' errmsg{iv}]);
        end
    end
   
    % I need to make sure the product I am using covers survey time of the
    % dropsonde circles. 
    % (define sampling window and sampling space, these two can be some tunning parameters)
    %
    % i. local conditions: in a domain of similar size as the JOANNE circle. during the flight/circle period 
    %    - options: add some leads and lags arround the sampling time
    %    period: +/- 1~3 hr of the circle mean time, will take a larger
    %    window for the investigation of lag time.
    % This is taking a long time, if I am reading a 6 hour window in a 1
    % min interval.(try 5 min interval)
    datasub = take_spatiotemporal_subsets(data, x0, y0, radius, t0, Twindow, NCFN, varn, tinc);
    
   
    % 
    % ii. background conditions: enlarge the space of sampling. (What is
    % the appropriate size of the images to apply the object-based
    % algorithm? ), if we need to apply Iorg, we will need to have a space
    % that is large enough.  (read literature to understand this)
    



end
    
    
    
    
    
    
    
    
    
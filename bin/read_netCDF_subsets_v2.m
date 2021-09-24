function datasub = read_netCDF_subsets_v2(x0, y0, radius, t0, trange, datadir, data_source, channel)
% purpose:
% input: x0, y0, trange ([1xN]) has the same length, even if x0 and y0
% remains the same during the sampling time range. 
%  - trange: datenum array specifying the range of time to extract the
%  satellite data.
% output:
%%% Note: this is an updated version (that needs to be improved)

% switch nargin
%     case 7
%         default_flag = opt;
%     case 6
%         default_flag = true;
%     otherwise
%         disp('warning: not enough input arguement!')
%         return
%         
% end

%% structure:
% start loop in time:
%   1. determine the source netcdf file to read in;
%   2. take spatiotemporal subsets of the netCDF source
% end loop

%% control satellite product name:
% ProdName.default = '0.5km_01min';
% %ProdName.default = '2km_10min';
% ProdName.backup ='2km_10min';
% ProdName.backup ='2km_10min_fulldisk';

ProdName = data_source;
% if strcmp(ProdName, '0.5km_01min')
%     ProdTag = 'default';
% elseif strcmp(ProdName, '2km_10min')
%     ProdTag = 'second';
% elseif strcmp(ProdName, '2km_10min_fulldisk')
%     ProdTag = 'second_fulldisk';
% end
    

%%
UTCoffset = -4;  % barbados in Atlantic Standard Time (GMT-4)
% information obtained from web search
local_sunrise_hr = 8;   % local sunrise time at BCO during DJF (7)
local_sunset_hr = 17;   % local sunrise time at BCO during DJF (17)

format1 = 'yyyy_mm_dd';
format2 = 'yyyymmdd';


% extract data within time range trange:
localdatehere = floor(t0);
dstr_tmp = datestr(localdatehere);
Nt = length(trange);

%radius = 1.5;  % units: degree
%Twindow = 6;    % units: hour;
%tinc = 5;      % units: 1 instances = 1 min (or 1 unit of temporal resolution.)

cnt =0;
for it =1:Nt
    
    
    tnow = trange(it);
    disp(['tnow:' datestr(tnow)]);
    
%     % determine which product to use based on the circle mean time:
%     if (tnow + UTCoffset/24) <= localdatehere + local_sunrise_hr/24 ...
%             || (tnow+ UTCoffset/24)>=localdatehere +local_sunset_hr/24
%         channel = 'IR';
%     else
%         channel = 'VIS';
%     end
    
    %     if localdatehere == datenum(2020,2,2,0,0,0)
    %         channel = 'IR';
    %     end
    
    % channel here should be selected based on the day (we know which are
    % the nighttime circles.)  --> (revised).
    if strcmpi(channel, 'VIS')
        channelID='02';
        varn='C02';
    else
        channelID='13';
        varn='C13';
    end
    
    
    % the default source:
    folderN = datestr(floor(tnow),format1);
    if strcmp(ProdName, '0.5km_01min')
        NCFN.default = [datadir filesep ProdName filesep '2020' filesep  folderN filesep ...
            'GOES_' channelID '_8N-18N-62W-50W_' datestr(floor(tnow),format2) '.nc'];
    else
        url0 = [datadir filesep ProdName filesep '2020' filesep folderN filesep 'catalog.html'];
        url_thredds  = replace(url0, 'dodsC','catalog');
        
        doy = floor(tnow) - datenum(2020,1,1)+1; % day of year;
        HH = datestr(tnow,'HH');
        MM = datestr(tnow,'MM');
        
        MMr = round(str2double(MM)/10)*10;
        if MMr ==60
            HHstr = num2str(str2double(HH)+1, '%2.2i');
            MMr = 0;
        else
            HHstr = HH;
        end
        MMr_str = num2str(MMr,'%2.2i');
        
        
        if strcmp(ProdName, '2km_10min')
            fileid = ['2020' num2str(doy,'%3.3i') HHstr MMr_str]; % YYYYJJJHHMMSSZ  (4digit year, 3digit day of year, hour, minute, second, tenth second
            disp(fileid);
            url_new = ['https://observations.ipsl.fr/aeris/eurec4a-data/SATELLITES/GOES-E/2km_10min/2020/' folderN filesep];

            fileid_full = get_fileID_from_thredds_catalog(url_new, fileid);
        %         fileid_full = get_fileID_from_thredds_catalog(url_new, fileid2);

            NCFN = [datadir filesep ProdName filesep '2020' filesep folderN filesep ...
                'clavrx_OR_ABI-L1b-RadF-M6C01_G16_s' fileid_full '_BARBADOS-2KM-FD.level2.nc'];
            
            datasub_tmp = take_subsets_from_secondary_source(NCFN, x0(it), y0(it), radius, channel);
            
            
        elseif strcmp(ProdName, '2km_10min_fulldisk')
            fileid2 = [datestr(tnow,'yyyymmdd') HHstr MMr_str];    % this is for the fulldisk data.
            NCFN = [datadir filesep ProdName filesep '2020' filesep folderN filesep ...
                'Emultic2kmNC4_goes16_' fileid2 '.nc'];
            datasub_tmp = take_subsets_from_secondary_source_v2(NCFN, x0(it), y0(it), radius, channel);
            
        end
        
    end
    
    
    % if there is no default source (error in the default source)
    % then read another datasource:
    
    % carve out a specific area from the dataset to get the circle location for this day: (if there are more than 1
    % circle for the day, then use the averaged center of the circle.)
    
    
    % I need to make sure the product I am using covers survey time of the
    % dropsonde circles.
    % (define sampling window and sampling space, these two can be some tunning parameters)
    %
    % i. local conditions: in a domain of similar size as the JOANNE circle. during the flight/circle period
    %    - options: add some leads and lags arround the sampling time
    %    period: +/- 1~3 hr of the circle mean time, will take a larger
    %    window for the investigation of lag time.
    
    
    
    % if ~isnan(fileid_full)
    
    
    
    
    %% note: I might need another subroutine to take out the values from the fulldisk dataset.
    
    cnt = cnt+1;
    % if cnt ==1
    
    %end
    datasub.time(cnt) = tnow;
    
    % make sure that the shape of right hand side matches that in the
    % left hand side.
    if cnt >1
        % potential need to trim the datasub_tmp.
        [ny1,nx1,~] = size(datasub_tmp.values);
        [ny0,nx0,~] = size(datasub.values);
        
        if nx0<nx1 || ny0<ny1
            % nx0_neg + nx0_pos = nx0
            % ny0_neg + ny0_pos = ny0
            
            % nx0_neg+ ny0_pos ~= nx_neg + ny_pos
            
            
            nlon_neg0 = length(find(datasub.lon(:,cnt-1)-x0(it-1)<0));
            nlon_pos0 = length(find(datasub.lon(:,cnt-1)-x0(it-1)>=0));
            
            nlat_neg0 = length(find(datasub.lat(:,cnt-1)-y0(it-1)<0));
            nlat_pos0 = length(find(datasub.lat(:,cnt-1)-y0(it-1)>=0));
            
            
            nlon_neg = length(find(datasub_tmp.lon - x0(it) <0));
            nlon_pos = length(find(datasub_tmp.lon - x0(it) >=0));
            
            nlondif_neg = nlon_neg-nlon_neg0;
            nlondif_pos = nlon_pos - nlon_pos0;
            
            
            nlat_neg = length(find(datasub_tmp.lat - y0(it) <0));
            nlat_pos = length(find(datasub_tmp.lat - y0(it) >=0));
            
            nlatdif_neg = nlat_neg - nlat_neg0;
            nlatdif_pos = nlat_pos - nlat_pos0;
            
            if nx0~=nx1 && ny0==ny1
                datasub_tmp.lon = datasub_tmp.lon(1+nlondif_neg: end-nlondif_pos);
                datasub_tmp.values = datasub_tmp.values(:,1+nlondif_neg: end-nlondif_pos);
                
            elseif ny0~=ny1 && nx0==nx1
                datasub_tmp.lat = datasub_tmp.lat(1+nlatdif_neg: end-nlatdif_pos);
                datasub_tmp.values = datasub_tmp.values(1+nlatdif_neg:end-nlatdif_pos,:);
                
            else  % nx0~=nx1 && ny0~=ny1
                datasub_tmp.lon = datasub_tmp.lon(1+nlondif_neg: end-nlondif_pos);
                datasub_tmp.lat = datasub_tmp.lat(1+nlatdif_neg: end-nlatdif_pos);
                datasub_tmp.values = datasub_tmp.values(1+nlatdif_neg: end-nlatdif_pos, 1+nlondif_neg: end-nlondif_pos);
                
            end
            
            % sanity check:
            [ny2, nx2, ~] = size(datasub_tmp.values);
            if nx2==nx0 && ny2==ny0
                disp('out: fixed length inconsistency');
            else
                disp('out: need to fix the code and forced to stop!');
                return
            end
            
            
        elseif nx0>nx1 || ny0<ny1
            % establish a NaN extra.
            nancol = NaN(1, nx0-nx1);
            nanrow = NaN(1,ny0-ny1);
            
            if nx0>nx1
                datasub_tmp.lon = [datasub_tmp.lon nancol];
                datasub_tmp.values  =[datasub_tmp.values repmat(nancol, ny1, 1)];
            elseif ny0>ny1
                datasub_tmp.lat = [datasub_tmp.lon nanrow];
                datasub_tmp.values  =[datasub_tmp.values; repmat(nanrow', 1, nx1)];
            else
                datasub_tmp.lon = [datasub_tmp.lon nancol];
                datasub_tmp.lat = [datasub_tmp.lon nanrow];
                datasub_tmp.values  =[datasub_tmp.values repmat(nancol, ny1, 1)];
                datasub_tmp.values  =[datasub_tmp.values; repmat(nanrow', 1, nx0)];
            end
        end
        
    end
    
    
    
    
    datasub.lon(:,cnt) = squeeze(datasub_tmp.lon);
    datasub.lat(:,cnt) = squeeze(datasub_tmp.lat);
    datasub.values(:,:,cnt) = datasub_tmp.values;
    datasub.channel{cnt}=channel;
    datasub.time(cnt)=tnow;
    
    
    
    
end







return



%%%%% save: 
%     if default_flag
%         if it ==1
%             % read basic lat, lon, and time info. from the source file.
%             try
%                 [data, errmsg] = read_ncfile_to_struct(NCFN.default);
%             catch  ME
%                 disp(ME.identifier)
%                 disp(ME.message)
%                 if cnotains(ME.message, '-90')
%                     % file not found, use another data source;
%                     disp('using another source file from 2km_10min folder');
%                     default_flag = false;
%                     
%                 end
%                 %             if strcmp(channel, 'IR')
%                 %                 % try the other source:
%                 %                 NCFN0 =  [datadir filesep datestr(datehere,format1) filesep ...
%                 %                     'GOES_02_8N-18N-62W-50W_' datestr(datehere,format2) '.nc'];
%                 %                 [data, errmsg] = read_ncfile_to_struct(NCFN0);
%                 %             end
%             end
%             
%             for iv = 1:length(errmsg)
%                 if ~isempty(errmsg{iv})
%                     disp(['warning: ' errmsg{iv}]);
%                 end
%             end
%         end
%         
%     end
%     


%  if default_flag
%         % read a subset of data:
%         window=0;   % only take the time at the given snapshot.
%         try
%             datasub_tmp = take_spatiotemporal_subsets(data, x0(it), y0(it), radius, tnow, window,  NCFN, varn, 1);
%             cnt = cnt+1;
%             if cnt ==1
%                 datasub.lon = datasub_tmp.lon;
%                 datasub.lat = datasub_tmp.lat;
%             end
%             datasub.time(cnt) = datasub_tmp.time;
%             datasub.values(:,:,cnt) = datasub_tmp.values;
%             datasub.channel{cnt}=channel;
%         catch
%             disp('error in reading the source netCDF..');
%             % read it from the secondary source:
%             %datasub_tmp = take_subsets_from_secondary_source(NCFN.backup, x0, y0, radius, tnow, channel);
%             
%         end
%         clear datasub_tmp
%         
%     else  %read from the second source:
%         % here, I need a new function to read the data:
%         disp('--> reading secondary source');
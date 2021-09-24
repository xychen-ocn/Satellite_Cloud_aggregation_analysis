function [data, errmsg] = read_ncfile_to_struct(FN)
% Purpose: read all the variables in the netCDF file and store it in a
%          matlab structure
% Input: FN (absolute path to a netCDF file)
% Output: data (in a structure format)

data_info = ncinfo(FN);
num_var = length(data_info.Variables);
errmsg{num_var} = [];   % preallocate memory for the cell array;

for iv = 1:num_var
    varn = data_info.Variables(iv).Name;
    try
        data.(varn) = ncread(FN , varn);
    catch
        errmsg{iv} = [varn ': fail to read data in the file'];
    end
end

% convert time:
time_info = ncreadatt(FN, 'time', 'units');
tmp = strsplit(time_info,' ');
basedate = datenum([tmp{end-1} ' ' tmp{end}]);
if strcmp(tmp{1}, 'seconds')
    data.time = data.time/86400 + basedate;
    
elseif strcmp(tmp{1}, 'minutes')
    data.time = data.time/(24*60) + basedate;
    
elseif strcmp(tmp{1}, 'hours')
    data.time = data.time/24 + basedate;
end




end
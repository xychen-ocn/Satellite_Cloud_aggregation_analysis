function fileIDstr = get_fileID_from_thredds_catalog(url, strin)
% purpose: this function read the catalog via url, and extract out the
% the dateID string information that is used to construct the netCDF file. 

webdata = webread(url);

expression = [strin '\d+'];    % any digits of a numeric number
tmp = regexp(webdata, expression,'match');

% take the first instance from tmp str because of repeatedness;
if ~isempty(tmp)
    fileIDstr = tmp{1};
else
    fileIDstr=NaN;
end

end
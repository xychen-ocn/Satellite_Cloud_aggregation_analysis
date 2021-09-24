function datasub = take_spatiotemporal_subsets(data, x0, y0, radius, t0, window, NCFN, varn, tinc)
% purpose: carve out a spatial box center at (x0, y0) in a square with a side of
% 2xradius; take the data in center at t0, with a give time window.
% (t-window/2: t+window/2)
% input: data  (a matlab structure containing coordinate, time and var)
% output: datasub (subsampled data in space and time)
%

lon_range = [x0-radius, x0+radius];
lat_range = [y0-radius, y0+radius];
time_range = [t0-window./2/24, t0+window/2/24];

[~, lon_stid] = min(abs(data.lon - lon_range(1)));
[~, lon_edid] = min(abs(data.lon - lon_range(2)));
[~, lat_stid] = min(abs(data.lat - lat_range(1)));
[~, lat_edid] = min(abs(data.lat - lat_range(2)));

[~, time_stid] = min(abs(data.time - time_range(1)));
[~, time_edid] = min(abs(data.time - time_range(2)));


% another way is to use mask;
% lon_mask = (data.lon>=lon_range(1)) .* (data.lon <=lon_range(2));
% lat_mask = (data.lat>=lat_range(1)) .* (data.lat <=lat_range(2));
% time_mask = (data.time>=time_range(1)) .* (data.time <=time_range(2));

% use the index to extract the subset;
datasub.lon = data.lon(lon_stid:lon_edid);
datasub.lat = data.lat(lat_stid:-1:lat_edid);
datasub.time = data.time(time_stid:tinc:time_edid);

% take a subset of the C02/C13 data:
% read from the NCFN again:
nt = time_edid - time_stid +1;
disp(['--> extracting data from '  datestr(data.time(time_stid)) ' to ' ...
      datestr(data.time(time_edid))]);
disp(['--> num of cloud scenes: ' num2str(ceil(nt/tinc)) ]); 
cnt = 0;
for it =time_stid:tinc:time_edid
    cnt = cnt + 1;
    tmp = ncread(NCFN, varn, [lon_stid, lat_edid, it],[lon_edid-lon_stid+1, abs(lat_edid-lat_stid)+1, 1]);
    datasub.values(:,:,cnt) = permute(tmp(:,end:-1:1),[2,1]);   % transpose/permute to NLAT X NLON
end

return
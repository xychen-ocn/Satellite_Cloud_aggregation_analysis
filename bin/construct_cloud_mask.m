function f=construct_cloud_mask(datasub)

%% calculate aggregation indices with an object-based algorithm from Dias et al. (2012)
%  refer to Janssens et al. (2021) paper, there are X object based
%  metrics that we can use to quantify the degree of aggregation from
%  the OLR.
Nt = size(datasub.values,3);

f = zeros(size(datasub.values));
range_data = max(datasub.values(:)) - min(datasub.values(:));

for it = 1:Nt
    channel = datasub.channel{it};
   % range_data = max(max(datasub.values(:,:,it))) - min(min(datasub.values(:,:,it)));
    %% note: how to deal with NaN...
    % calculate area that is occupied by NaN and remove those area when
    % calculating the statistics.
    
    % scale the data values to be within range 0~1 (relative magnitude). the scaling is
    % different for VIS and IR input data;
    if strcmpi(channel, 'VIS')
        %datain = satedata(:,:,it)./max(max(satedata(:,:,it)));
        % cloudy pixels associates with larger reflectance. increasing
        % value --> clouds
       % f(:,:,it) = (datasub.values(:,:,it)-min(min(datasub.values(:,:,it))))./range_data;
         f(:,:,it) = (datasub.values(:,:,it)-min(datasub.values(:)))./range_data;
        % scale data by its range:
    else
        % brightness temperature: lower --> cloudy
        % after normalization: 0=no cloud, 1, cloudy.
        %f(:,:,it) = (max(max(datasub.values(:,:,it)))-datasub.values(:,:,it))./range_data;
        f(:,:,it) = abs(datasub.values(:,:,it)-max(datasub.values(:)))./range_data;

    end
    
    
    
end


return
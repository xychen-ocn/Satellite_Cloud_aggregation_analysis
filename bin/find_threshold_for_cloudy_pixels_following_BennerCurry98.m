function thres_val = find_threshold_for_cloudy_pixels_following_BennerCurry98(datasub)

% compute threshould for all the scenes in a segment:
channel = datasub.channel{1};
all_data = datasub.values(:);

if strcmp(channel,'IR')
    max_IR= mode(all_data(:));
    dx = max(all_data(:)) - max_IR;
    lower_IR = floor(max_IR - dx);
    h = histogram(all_data(all_data>lower_IR),'normalization','pdf');
    pd = fitdist(all_data(all_data>lower_IR), 'Normal');
    %pd = makedist('Normal','mu', 295 ,'sigma',0.4);
    x_values = [285:0.1:297];
    y = pdf(pd, x_values);
    hold on
    plot(x_values,y,'-r');
    hold off
    pdf_thres = 0.005*max(y);
    thres_val = max(x_values((y<=pdf_thres)&(x_values<max_IR)));

    
    
else
    
    max_VIS= mode(all_data(:));
    dx =  max_VIS - min(all_data(:));
    upper_VIS = floor(max_VIS + dx);
    window = (all_data > 2 & all_data<upper_VIS);
    h = histogram(all_data(window),'normalization','pdf');
    h = histogram(all_data,'normalization','pdf');
    pd = fitdist(all_data(window), 'Normal')
    %pd = makedist('Normal','mu', 295 ,'sigma',0.4);
    x_values = [0:0.1:upper_VIS*1.2];
    y = pdf(pd, x_values);
    hold on
    plot(x_values,y*0.1,'-r');
    hold off
    % xbins = 0.5*(h.BinEdges(1:end-1) + h.BinEdges(2:end));
    % plot(xbins, h.Values,'-k');
    pdf_thres = 0.005*max(y);
    thres_val = min(x_values((y<=pdf_thres)& (x_values>max_VIS)));
    
    
    
end

% for it = 1:Nt
% cldmask(:,:,it) = datasub.values(:,:,it) - thres_val;
% end


return
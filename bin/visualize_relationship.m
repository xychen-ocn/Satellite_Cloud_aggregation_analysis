function visualize_relationship(X, IndStruct, Xedges, svparm)
% purpose: 
% input: X
% .      Y (structure containing indices)'
%        Xedges (Xbin edges)
%        svparm.flag , svparm.figsvdir

indx_name = fieldnames(IndStruct);
%Xedges = [0:5:70];

% group data by bins and then perform box plot:
[X_sorted, sid] = sort(X);
for iv = 1:length(indx_name)
    varn = indx_name{iv};
    sorted.(varn) = IndStruct.(varn)(sid);
end

[BN,Edg]=discretize(X_sorted, Xedges(1:2:end));
BinCen = (Edg(1:end-1) + Edg(2:end))*0.5;

for ib = 1:length(BinCen)
    BinCenLabel{ib} = num2str(BinCen(ib));
end

% group by bins:
BG = nan(size(BN));
for ig = 1:length(BinCen)
    tmp= find(BN==ig);
    if ~isempty(tmp)
       BG(tmp) = BinCen(ig);
    end
end


ps = customize_subplot_size(2,1,0.1,0.05);
for iv = 1:length(indx_name)
    varn = indx_name{iv};
    figure(iv); clf;
    
    %% 2D histogram
    subplot(2,1,1);
    
    % divide the Y range into 10 sections;
    Yedges = linspace(round(min(IndStruct.(varn))), round(max(IndStruct.(varn))), 11);
    histogram2(X, IndStruct.(varn), Xedges, Yedges,'FaceColor','Flat');
    hb = colorbar;
    view(2)
    %xlabel('Concave Area')
    ylabel(varn);
    title(varn)
    %caxis([0 ]);
    ytick_ax1 = get(gca,'ytick');
    set(gca,'pos',ps{1});
    
    %% 1D boxplot:
    subplot(2,1,2);
    hold on;
    % bin averaged boxplot:
    boxplot(sorted.(varn),BN);
%     hb2 = colorbar;
%     set(hb2,'visible','off');
    %set(gca,'ytick', ytick_ax1);
    xlabel('Concave Area')
    ylabel(varn);
    %title(varn)
    set(gca,'xticklabel',BinCenLabel);
    %get(gca,'xtick')
    xlim([0 length(Edg)]);
    hold off
    set(gca,'pos',ps{2});
    % save the current figure
    if svparm.flag
        xc_savefig(gcf,svparm.figsvdir,[svparm.figprefix varn '.jpg'],[0 0 10 6]);
    end
    %caxis([0 8]);
   % pause
end



return
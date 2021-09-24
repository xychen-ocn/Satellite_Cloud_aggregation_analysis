function [blob_coord,blob_image, f1, stats_sel] = find_image_2d_updated(f,f0)

% INPUT
% f = [lon,time] (or [lat, lon])
% f0 = threshold
% datemin = [year month day hr min sec]

% OUTPUT
% blob_coord = [centr_time centr_lon width_time width_lon]
% !!!!!!!!all nondimensional (i.e. # of grid points) !!!!!!!!!!!!
% blob_image is a structure where blob_image{n} corresponds to blob_coord(n,:)


fo = f; % save original data, just for ploting at the end;
NConn = 8;   %  Janssens et al. (2021) used 6-connectivity??

[Ny,Nx]=size(f);

BI=(f>f0); %matrix with ones where f excceeds threshold
CCr = bwconncomp(BI,NConn); % defines the connected components of B
Lr = labelmatrix(CCr); % each integer corresponds to one of the connected components
num_events_orig = max(Lr(:));

fprintf(1,[int2str(num_events_orig),' - ini number of events\n']);

stats = regionprops(Lr,f,'BoundingBox','WeightedCentroid',...
          'Image','Orientation','Eccentricity','MaxFeretProperties',...
          'EquivDiameter','EulerNumber', 'Perimeter','Area', ...
           'MajorAxisLength','MinorAxisLength', 'PixelIdxList');           % add returning pixel indices (p-element vector)
      
box=[stats.BoundingBox];
cent=[stats.WeightedCentroid];
orien = [stats.Orientation];
eccen = [stats.Eccentricity];
if ~isempty(stats)
    maxFeretAng = [stats.MaxFeretAngle];
    maxFeretCoord = reshape([stats.MaxFeretCoordinates],4,[]);
    maxFeretDiam = [stats.MaxFeretDiameter];
end
EqvDiam = [stats.EquivDiameter];
EuNum = [stats.EulerNumber];
Perim = [stats.Perimeter];
area = [stats.Area];
major_axislen = [stats.MajorAxisLength];
minor_axislen = [stats.MinorAxisLength];

for i = 1:num_events_orig
    pixel_idx{i} = stats(i).PixelIdxList;
end


box_candidates=reshape(box,4,max(Lr(:))); 

% ini_t=box_candidate(1,event) ; ini_x=box_candidate(2,event) 
% width_t=box_candidate(3,event);  width_x=box_candidate(4,event) 

cent_candidates=reshape(cent,2,max(Lr(:))); 
blob_image0 = {stats.Image};

KMin= 4; % min # of grid points in time
NMin= 4; % min # of grid points in long

% KMin = prctile(major_axislen, 25);
% NMin = prctile(major_axislen, 25);
% disp(['-crit: Kmin, Nmin=', num2str(KMin), num2str(NMin)]);

% in Janssens' work, the area needs to be larger than the area of 4 pixels.
% with 4-connectivity
cond_area = area>8;
%cond_area = area>prctile(area,25);
%disp(['-crit: area=', num2str(prctile(area,25))]);


cond_size = (box_candidates(3,:)>KMin ...
            & box_candidates(4,:)>NMin);
        
%%%% ----- construct a database that stores the pixel indices for edges -----
left_edge = [1:Ny];
top_edge = 1:Ny:(Nx-1)*Ny+1;
right_edge = (Nx-1)*Ny:1:Nx*Ny;
bottom_edge = Ny:Ny:Nx*Ny;

edge_matrix = unique([left_edge, top_edge, right_edge, bottom_edge]);

cond_edge = logical(size(cond_size));
for i = 1:num_events_orig
    % if any pixel in the originally identified object contains edge
    % pixels.
    cond_edge(i) = any(ismember(pixel_idx{i}, edge_matrix  ));       % dimension: [1x num_of_events]
end

%cond_lon = (box_candidates(2,:)> 1 & (box_candidates(2,:)+box_candidates(4,:))<N); % have to deal with periodicity

 idx_events = find(cond_area  & ~cond_edge ); % keep only cases that satisfy the two conditions above
 %idx_events = find(cond_area);
% idx_events = find( ~cond_edge);
%idx_events = find(cond_area &  cond_lon); % keep only cases that satisfy the two conditions above


edge_events = find(cond_area  &cond_edge);


if ~isempty(idx_events)
    num_events = max(size(idx_events));
else
    num_events=0;
end
fprintf(1,[int2str(num_events),' -  actual number of events\n']);

blob_image1 = blob_image0(idx_events);
blob_box1 = box_candidates(3:4,idx_events)'; %(timw lonw) . (lon, lat)
blob_cent1 = cent_candidates(:,idx_events)'; %(timcent loncent) (lon,lat)
% add extra quantities here:
stats_sel.eccen = eccen(idx_events);
stats_sel.orientation  = orien(idx_events);
if ~isempty(stats)
    stats_sel.maxFeretAng = maxFeretAng(idx_events);
    stats_sel.maxFeretDia = maxFeretDiam(idx_events);
    stats_sel.maxFeretCoord = maxFeretCoord(:,idx_events);
end
stats_sel.EqvDiam = EqvDiam(idx_events);
stats_sel.EuNum = EuNum(idx_events);
stats_sel.Perimeter = Perim(idx_events);
stats_sel.CloudArea = area(idx_events);
stats_sel.Area_removed = sum(area(edge_events));
stats_sel.EqvDiam_EdgeClouds = EqvDiam(edge_events);
stats_sel.coord_EdgeClouds = [cent_candidates(:,edge_events)', box_candidates(3:4,edge_events)'];

% I need to save the area of objects I dropped out to compute the cloud
% density; 

BW1 = ismember(Lr, idx_events);
f1=(BW1>0);% for ploting only

blob_coord1 = [blob_cent1 blob_box1];

% --- Note: I don't need to deal with periodicity 
%lon_check = isempty((box_candidates(2,:)==1 | (box_candidates(2,:)+box_candidates(4,:))>=N));
% deals with blobs crossing 0^o
% if lon_check==0
%     
%     f = circshift(f,[floor(N/2) 0]);
%     lon=1:N;
%     lonc=circshift(1:N,[0 floor(N/2)]);
%     
%     BI=(f>f0); %matrix with ones where f excceeds threshold
%     CCr = bwconncomp(BI,6); % defines the connected components of B
%     Lr = labelmatrix(CCr); % each integer corresponds to one of the connected components
%     num_events_orig = max(Lr(:));
% 
%     fprintf(1,[int2str(num_events_orig),' - ini number of events crossin 0^o\n']);
% 
%     stats = regionprops(Lr,f,'BoundingBox','WeightedCentroid',...
%         'Image');
%     box=[stats.BoundingBox];
%     cent=[stats.WeightedCentroid];
%     box_candidates=reshape(box,4,max(Lr(:))); 
%     cent_candidates=reshape(cent,2,max(Lr(:))); 
%     blob_image0 = {stats.Image};
% 
%     cond_size = (box_candidates(3,:)>KMin ...
%                 & box_candidates(4,:)>NMin);
%     cond_lon = (box_candidates(2,:)< N/2 & (box_candidates(2,:)+box_candidates(4,:))>N/2);
% 
%     idx_events = find(cond_size &  cond_lon); % keep only cases that staisfy the two conditions above
% 
%     num_events = max(size(idx_events));
%     fprintf(1,[int2str(num_events),' -  actual number of events crossing 0^o\n']);
% 
%     blob_image2 = blob_image0(idx_events);
%     blob_box2 = box_candidates(3:4,idx_events)'; %(timw lonw)
%     blob_cent2 =  cent_candidates(:,idx_events)'; %(timcent loncent)
%     blob_cent2(2,:) = lonc(round(blob_cent2(2,:)));
%     blob_coord2 = [blob_cent2 blob_box2];
%     blob_image = [blob_image1 blob_image2];
%     blob_coord = [blob_coord1' blob_coord2']';
%     
% else
    blob_image = blob_image1;
    blob_coord = blob_coord1'';
    
% end

% PLOT CHECK
figure(10); 
contourf(fo','edgecolor','none');caxis([0,1]);colorbar;colormap((gray(12)));
hold on;
if num_events>0
contour(f1','b');
hold on;
plot(blob_coord(:,2),blob_coord(:,1),'r+','markersize',10);
end
xlabel('lon (# grid points)');ylabel('lat(# grid points)');
title('f, blobs, and centroids')
hold off;
end
 
 
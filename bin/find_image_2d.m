function [blob_coord,blob_image] = find_image_2d(f,f0)

% INPUT
% f = [lon,time]
% f0 = threshold
% datemin = [year month day hr min sec]

% OUTPUT
% blob_coord = [centr_time centr_lon width_time width_lon]
% !!!!!!!!all nondimensional (i.e. # of grid points) !!!!!!!!!!!!
% blob_image is a structure where blob_image{n} corresponds to blob_coord(n,:)


fo = f; % save original data, just for ploting at the end;

[N,K]=size(f);

BI=(f>f0); %matrix with ones where f excceeds threshold
CCr = bwconncomp(BI,6); % defines the connected components of B
Lr = labelmatrix(CCr); % each integer corresponds to one of the connected components
num_events_orig = max(Lr(:));

fprintf(1,[int2str(num_events_orig),' - ini number of events\n']);

stats = regionprops(Lr,f,'BoundingBox','WeightedCentroid',...
    'Image');
box=[stats.BoundingBox];
cent=[stats.WeightedCentroid];

box_candidates=reshape(box,4,max(Lr(:))); 

% ini_t=box_candidate(1,event) ; ini_x=box_candidate(2,event) 
% width_t=box_candidate(3,event);  width_x=box_candidate(4,event) 

cent_candidates=reshape(cent,2,max(Lr(:))); 
blob_image0 = {stats.Image};

KMin=4; % min # of grid points in time
NMin=4; % min # of grid points in long

cond_size = (box_candidates(3,:)>KMin ...
            & box_candidates(4,:)>NMin);
    
cond_lon = (box_candidates(2,:)> 1 & (box_candidates(2,:)+box_candidates(4,:))<N); % have to deal with periodicity

     
idx_events = find(cond_size &  cond_lon); % keep only cases that satisfy the two conditions above

num_events = max(size(idx_events));
fprintf(1,[int2str(num_events),' -  actual number of events\n']);

blob_image1 = blob_image0(idx_events);
blob_box1 = box_candidates(3:4,idx_events)'; %(timw lonw)
blob_cent1 = cent_candidates(:,idx_events)'; %(timcent loncent)

BW1 = ismember(Lr, idx_events);f1=(BW1>0);% for ploting only

blob_coord1 = [blob_cent1 blob_box1];

lon_check = isempty((box_candidates(2,:)==1 | (box_candidates(2,:)+box_candidates(4,:))>=N));
% deals with blobs crossing 0^o
if lon_check==0
    
    f = circshift(f,[floor(N/2) 0]);
    lon=1:N;
    lonc=circshift(1:N,[0 floor(N/2)]);
    
    BI=(f>f0); %matrix with ones where f excceeds threshold
    CCr = bwconncomp(BI,6); % defines the connected components of B
    Lr = labelmatrix(CCr); % each integer corresponds to one of the connected components
    num_events_orig = max(Lr(:));

    fprintf(1,[int2str(num_events_orig),' - ini number of events crossin 0^o\n']);

    stats = regionprops(Lr,f,'BoundingBox','WeightedCentroid',...
        'Image');
    box=[stats.BoundingBox];
    cent=[stats.WeightedCentroid];
    box_candidates=reshape(box,4,max(Lr(:))); 
    cent_candidates=reshape(cent,2,max(Lr(:))); 
    blob_image0 = {stats.Image};

    cond_size = (box_candidates(3,:)>KMin ...
                & box_candidates(4,:)>NMin);
    cond_lon = (box_candidates(2,:)< N/2 & (box_candidates(2,:)+box_candidates(4,:))>N/2);

    idx_events = find(cond_size &  cond_lon); % keep only cases that staisfy the two conditions above

    num_events = max(size(idx_events));
    fprintf(1,[int2str(num_events),' -  actual number of events crossing 0^o\n']);

    blob_image2 = blob_image0(idx_events);
    blob_box2 = box_candidates(3:4,idx_events)'; %(timw lonw)
    blob_cent2 =  cent_candidates(:,idx_events)'; %(timcent loncent)
    blob_cent2(2,:) = lonc(round(blob_cent2(2,:)));
    blob_coord2 = [blob_cent2 blob_box2];
    blob_image = [blob_image1 blob_image2];
    blob_coord = [blob_coord1' blob_coord2']';
    
else
    blob_image = blob_image1;
    blob_coord = blob_coord1'';
    
end

% PLOT CHECK
figure;contourf(fo','edgecolor','none');caxis(caxis);colorbar;colormap(flipud(gray(12)));
hold on;contour(f1',[0 0],'b');
hold on;plot(blob_coord(:,2),blob_coord(:,1),'r+','markersize',10);
xlabel('lon (# grid points)');ylabel('time (# grid points)');
title('f, blobs, and centroids')
end
 
 
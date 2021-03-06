function pt_nnb = find_nearest_neighbor_pairs(X)
% purpose: find nearest neighbor pairs in the input matrix. 
% input: X [N*2]: N or records (centroids, col1=X, col2=Y)
% output: nnb_pair (with field idxs and dist)

for i=1:size(X,1)
    query_point = X(i,:);
    point_cloud = X;
    point_cloud(i,:)=[];
    [idx, Dist]=knnsearch(point_cloud, query_point);
%     figure(10); %clf;
%     plot(point_cloud(:,1), point_cloud(:,2),'.k');
%     hold on
%     plot(query_point(1), query_point(2),'*r');
%     plot(point_cloud(idx,1), point_cloud(idx,2),'or');
%     plot([point_cloud(idx,1), query_point(1)], [point_cloud(idx,2), query_point(2)],'-w');
%    % quiver(query_point(1), query_point(2), 
    %hold off
    
    pt_nnb.x(i) = point_cloud(idx,1);
    pt_nnb.y(i) = point_cloud(idx,2);
    tmp = find(X==point_cloud(idx));
    pt_nnb.idx(i)=tmp(1);
    pt_nnb.dist(i)=Dist;
end


% remove repeated pair of neareast neighbors:
% cnt=0;
% clear nnb_pair;
% for i =1:size(X,1)
%     nbid = pt_nnb.idx(i);
%     
%     if ~isnan(nbid)
%         cnt=cnt+1;
%         nnb_pair.idxs(cnt,:)=[i, nbid];
%         nnb_pair.dist(cnt) = pt_nnb.dist(i);
%         
%         % see if these two points are mutual neighbors:
%         % if so, remove the repeated pair.
%         
%         if pt_nnb.idx(nbid)==i
%             % repeated:
%             % cross over: by NaN
%             pt_nnb.idx(nbid)=NaN;
%             
%         end
%     else
%         continue
%     end
%     
% end


return


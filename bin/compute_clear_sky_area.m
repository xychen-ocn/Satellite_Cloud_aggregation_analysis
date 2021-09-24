function Aopen_max = compute_clear_sky_area(cldmask, Avalid)

% Referece: Antonissen Master's thesis Page 48

% The largest clear sky region is used as an additional measure for
% organization/clustering.

% A larger clear sky region indicates more clustering.

% definition: Aclear,max = sqrt(max{Aclear})/sqrt(A);  A is the area with
% good data; range of Aclear,max = (0,1);


% 1. find all the noncloud pixel location:
non_cloud_idx = find(~cldmask);
[Ny, Nx] = size(cldmask);

for i = 1:length(non_cloud_idx)
    
    % start searching:
    id = non_cloud_idx(i);
    
    % convert back to the subscripts for the 2D matrix:
    [iy, ix] = ind2sub(size(cldmask), id);
    
    % find the cloesest cloud pixel on the same row:
    %
    % check if there is cloudy pixel to the left (west side)   
    search_result = find(cldmask(iy,1:ix));
    if  ~isempty(search_result)
        Lx_west = ix - max(search_result);
        
    else
        % compute distance to edge:
        Lx_west = ix - 1;
    end
    
    
    % check if there is cloudy pixel to the right (east side)   
    search_result = find(cldmask(iy,ix:end));
    xs = ix:Nx;
    if  ~isempty(search_result)
        Lx_east = min(xs(search_result))-ix;
        
    else
        % compute distance to edge:
        Lx_east = Nx - ix;
    end
    Lx = [Lx_west, Lx_east];
    
    
    %% repeat for the vertical direction:
    % check if there is cloudy pixel upward (north side)   
    search_result = find(cldmask(1:iy,ix));
    if  ~isempty(search_result)
        Ly_south = iy - max(search_result);
        
    else
        % compute distance to edge:
        Ly_south = iy - 1;
    end
    
    
    % check if there is cloudy pixel downward (south side)   
    search_result = find(cldmask(iy:end,ix));
    ys = iy:Ny;
    if  ~isempty(search_result)
        Ly_north = min(ys(search_result))-iy;
        
    else
        % compute distance to edge:
        Ly_north = Ny - iy;
    end
    Ly = [Ly_south, Ly_north];
    
%     if any(Lx==0) 
%         Lx = max(Lx);
%     else
%         Lx = min(Lx);
%     end
%     
%     if any(Ly==0)
%         Ly = max(Ly);
%     else
%         Ly = min(Ly);
%     end
    
        Aopen(i) = 2*min(Lx) * 2*min(Ly);
        pixel_loc(i,:) = [ix, iy];
        Lx_save(i,:) = Lx;
        Ly_save(i,:) = Ly;
        
        % plot it out:
%         figure(13);
%         pcolor(cldmask); shading flat;
%         colorbar;
%         hold on;
%         plot(ix, iy, 'sg','MarkerSize',10);
%         quiver(ix,iy, -Lx(1), 0, 'off','r','linewidth',1.2);
%         quiver(ix,iy, Lx(2), 0, 'off','r','linewidth',1.2);
%         quiver(ix,iy, 0,-Ly(1), 'off','r','linewidth',1.2);
%         quiver(ix,iy,0, Ly(2), 'off','r','linewidth',1.2);
%         hold off;
        
       % pause(0.01)
        

end

Aopen_max = sqrt(max(Aopen))./sqrt(Avalid);



%% plot the maximum open sky location:
[~, maxid] = max(Aopen);
ix = pixel_loc(maxid,1); iy =pixel_loc(maxid,2);
Lx = Lx_save(maxid,:);
Ly = Ly_save(maxid,:);
% plot it out:
figure(13);
pcolor(cldmask); shading flat; 
colormap(gray);
colorbar;
hold on;
plot(ix, iy, 'sg','MarkerSize',10);
quiver(ix,iy, -Lx(1), 0, 1,'r','linewidth',1.2);
quiver(ix,iy, Lx(2), 0, 1,'r','linewidth',1.2);
        quiver(ix,iy, 0,-Ly(1), 'off','r','linewidth',1.2);
        quiver(ix,iy,0, Ly(2), 'off','r','linewidth',1.2);
% plot a square
plot([ix-min(Lx), ix+min(Lx), ix+min(Lx), ix-min(Lx), ix-min(Lx)], ...
    [iy-min(Ly), iy-min(Ly) ,iy+min(Ly), iy+min(Ly), iy-min(Ly)],'--r', 'linewidth',1.2);
hold off;
% plot a square



return
function Iorg_i = compute_iorg_inhibition(scene_size, Centroids, R)

% Referece: Antonissen Master's thesis Page 48
% Here, periodic boundary condition is not used. 

LX = scene_size(2) ;
LY = scene_size(1);
    
% 1. sort the cloud object by its radius:
[R, sid] = sort(R, 'descend');
cx = Centroids(sid,1);
cy = Centroids(sid,2);

% 2. start generating random distribution of the cloud circles:
cnt = 0;
place_again = true;
placeable_flag = false;
occupied_x = []; occupied_y = occupied_x;
occupied_rad = [];

figure(15); clf;

while place_again
    cnt = cnt+1;
    
 
        % initial place ment:
            
        [xr, yr] = generate_random_cloud_loc;
        % place the largest cloud in the scene:
        %dist_to_edge = xr<R(1) || abs(xr-LX)<R(1) || yr<R(1) || abs(yr-LX)<R(1);
        while (xr<R(1) || abs(xr-LX)<R(1) || yr<R(1) || abs(yr-LY)<R(1))
            [xr, yr] = generate_random_cloud_loc;
            %dist_to_edge = xr<R(1) || abs(xr-LX)<R(1) || yr<R(1) || abs(yr-LX)<R(1)
        end
        occupied_x = [occupied_x, xr];
        occupied_y = [occupied_y, yr];
        
        occupied_rad = [occupied_rad, R(1)];
        
        
        ii = 2;
        while ii <= length(R)
            itr = 0;
            
            while ((itr<1000) & (~placeable_flag))
                [xr, yr] = generate_random_cloud_loc;
                
                % check if cloud exist in this location:
                if ismember(xr, occupied_x) & ismember(yr, occupied_y)
                    % cloud exist;
                    [xr, yr]=generate_random_cloud_loc();
                    
                else
                    % check if there is overlapping clouds:
                    dist = sqrt((xr - occupied_x).^2 + (yr-occupied_y).^2);
                    dist_diff=[];
                    for j = 1:size(occupied_x,2)
                        dist_diff(j) = dist(j) - (R(ii)+occupied_rad(j));
                    end
                    edge = (xr<R(ii) || abs(xr-LX)<R(ii) || yr<R(ii) || abs(yr-LY)<R(ii));
                    
                    if any(dist_diff<0) || edge
                        % overlapped cloud
                        placeable_flag = false;
                        itr = itr +1;
                        
                    else
                        occupied_x = [occupied_x, xr];
                        occupied_y = [occupied_y, yr];
                        occupied_rad = [occupied_rad, R(ii)];
                        placeable_flag = true;
                        if ii==length(R)
                            place_again = false;
                        end
                        
%                         figure(15);
%                         plot(occupied_x, occupied_y,'+k');
%                         hold on;
                        
%                         for k= 1:length(occupied_x)
%                             circle(occupied_x(k), occupied_y(k), occupied_rad(k),'r',1.2);
%                         end
                        
                        if ii==2
                            title('random placing circular clouds');
                            xlabel('NX (pixel)');
                            ylabel('NY (pixel)');
                        end
                        %pause
                        
                        
                    end
                end    
            end
            
            if placeable_flag
                ii = ii+1;
               % disp(num2str(ii));
                placeable_flag=false;
            else
                break
            end
        end
                    
            

    if cnt > 10
        place_again= false;
        %if itr>=1000 & ~placeable_flag
        Iorg_i=NaN;
        return
        %end
    end
end


% if made it:
if length(occupied_x)==length(R)
    random_centroids = [occupied_x', occupied_y'];
    
    nbins = 100000;
    % start computing the NNCDF for the randomly placed circle:
    
    maxDist =  sqrt(LX^2 + LY^2);
    bins = linspace(0, maxDist, nbins);
    
    nnb_pair_rand= find_nearest_neighbor_pairs(random_centroids);
    nnb_pair_scene= find_nearest_neighbor_pairs(Centroids);
    
    figure(12);
    hrand=histogram(nnb_pair_rand.dist,bins,'Normalization','cdf');
    
    figure(12); hold on
    hscen=histogram(nnb_pair_scene.dist,bins,'Normalization','cdf');

    
    NNCDF_rand = hrand.Values;
    NNCDF_scene = hscen.Values;
    %binEdge = h.BinEdges;
    %binave = 0.5*(bins(1:end-1) + bins(2:end));
    
    Iorg_i = trapz(NNCDF_rand, NNCDF_scene);
    figure(12);clf;
    plot(NNCDF_rand, NNCDF_scene,'-k','linewidth',1.2);
    xlabel('random NNCDF');
    ylabel('scene NNCDF');
    hold on
    plot([0:0.1:1],[0:0.1:1],'--b','linewidth',1.2);
    hold off
    axis([0 1 0 1]);
    set(gca,'fontsize',14);
    

    
    
end






% nested function:
function [xr, yr] = generate_random_cloud_loc()
    x0 = rand;
    y0 = rand;
    while (x0==0 | y0==0)
        % regenerate coordinate:
        x0 = rand;
        y0 = rand;
    end
    
    xr = round(x0*LX);
    yr = round(y0*LY);
  

end
end
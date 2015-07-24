function [L_dividedD1,L_dividedD7,L_dividedD14,V_merge_D7,V_merge_D14]  = return_L_divided(L_D1,L_D7,L_D14)
% returnL_divided returns three labeled matrix where metastases and primary tumor have been seperated (if the had merged previos to time t)

V_merge=[];

[cmD1,indexD1,cmD7,indexD7,cmD14,indexD14] = return_met_events(L_D1,L_D7,L_D14);

BWD1 = L_D1;
BWD1(BWD1>0)=1;
BWD7 = L_D7;
BWD7(BWD7>0)=1;
BWD14 = L_D14;
BWD14(BWD14>0)=1;

DMD1 = return_dist_matrix(BWD1);
DMD7 = return_dist_matrix(BWD7);
DMD14 = return_dist_matrix(BWD14);


DCD1 = 0;% returnDilationCount(L_D1);
DCD7 = 0;%returnDilationCount(L_D7);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D1
% L_D1 never needs to be seperated...
L_dividedD1 = L_D1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D7

% find out it L_D7 need to be seperated
cmlistD7 = [cmD1;cmD7];
S = size(cmlistD7);

tempD1 = L_D1;
tempD7 = L_D7;
tempD7(tempD1==0)=0;


        rL = regionprops(bwlabel(BWD7),'Centroid');   % get CM's of blobs in this time point
        rC = regionprops(L_dividedD1,'Centroid'); % get CM's of objects in previous time point

        % list of indexes of blops in labeled image
        listBlobsNotMerge = zeros(1,length(rC)); % length(rC) is number of clusters/objects at t-1

        % for each object at t-1 find the closest blob at t
        for k=1:length(rC)
            dist = 1e10;
            index = 0;
            for n=1:length(rL)
                
                cmL = rL(n).Centroid; % (Note 'Centroid' returns first x = col then y = row)
                cmC = rC(k).Centroid; % (Note 'Centroid' returns first x = col then y = row)
                
                temp_dist = sqrt((cmL(1)-cmC(1))^2 + (cmL(2)-cmC(2))^2 ); % distance between n'th blop at time+t and k'th object at time=t-1
                
                if dist>temp_dist && temp_dist<55 % save smallest distance
                    dist = temp_dist;
                    index = n;
                end
            end
            listBlobsNotMerge(k) = index; % put index of closest blop D7 in list
        end
        

L_D7DividedYesNo = 0;

% if number met events (S(1)) is 1 or if number of indexses is the same for
% L_D1 and L_D7(L_D1==0)=0 then L_D7 needs no division
if  S(1)==1 || length(unique(nonzeros(tempD1)))==length(unique(nonzeros(tempD7))) 
    L_dividedD7 = L_D7;
    
elseif  length(listBlobsNotMerge)==length(unique(listBlobsNotMerge)) && isempty(find(listBlobsNotMerge==0))==1 % if number of unique blop indexes is same as number of object no blop has has merged!
    
    [L_dividedD7,numMetEvents ] = return_clustered_C(BWD7,return_dist_matrix(BWD7),20,listBlobsNotMerge);
    numMetEvents;
    L_D7DividedYesNo = 1;
    
else
    L_D7DividedYesNo = 1;
    
    objects =[];
    objectCounter=0;
    
    % iterate through metastatic events which happened before D1
    for k=1:length(indexD1)
        
        objectCounter=objectCounter+1;
        
        L_k = return_sub_listL(L_D1,indexD1(k));
        L_k = bwmorph(L_k,'dilate',DCD1);
        [xx,yy] = return_object_voronoi_points(L_k);
        
        objects(objectCounter).voronoiPointsX = xx;
        objects(objectCounter).voronoiPointsY = yy;
        
    end
    
    for k=1:length(indexD7)
        
        objectCounter=objectCounter+1;
        
        L_k = return_sub_listL(L_D7,indexD7(k));
        
        [xx,yy]= return_object_voronoi_points(L_k);
        
        objects(objectCounter).voronoiPointsX = xx;
        objects(objectCounter).voronoiPointsY = yy;
        
    end
    
    % extract all points into one vector so they can be feed to next function
    
    x=[];
    y=[];
    for k=1:length(objects)
        x = [x; objects(k).voronoiPointsX];
        y = [y; objects(k).voronoiPointsY];
    end
    
    %%%%
    
    [rows,cols] = size(L_D7);
    
    
    
    V_labeled = return_labeled_voronoi(x,y,rows,cols);
    
    %figure
    %imshow(V_labeled,[])
    %colormap(jet)
    
    
    
    V_merge = zeros(size(V_labeled));
    
    for k=1:(length(objects)) % iterate through all objects in D7
        
        tempList=[]; % list of voronoi regions i V_labeled which should be joined
        pointsX = objects(k).voronoiPointsX;
        pointsY = objects(k).voronoiPointsY;
        for n=1:length(objects(k).voronoiPointsX) % iterate through voronoi points x,y for object k
            tempList(n)= V_labeled(round(pointsY(n)),round(pointsX(n))); % add index of region to list
        end
        
        % make empty matrix with 0's
        tempLL = zeros(size(V_labeled));
        
        % fill in 1's in the regions which need to be joined
        for g=1:length(tempList)
            tempLL(V_labeled==tempList(g))=1;
        end
        
        % dilate so regions merge
        tempLL = bwmorph(tempLL,'dilate',2);
        % Attempt to close hole by filling in white in edges
        tempLL(5:rows-5,1)=1;
        tempLL(5:rows-5,cols)=1;
        tempLL(1,2:cols-2)=1;
        tempLL(rows,2:cols-2)=1;
        
        % fill holes and erode
        tempLL = imfill(tempLL,'holes');
        tempLL = bwmorph(tempLL,'erode',2);
        
%      figure
%      imshow(tempLL,[])
%      hold on
%      plot(pointsX,pointsY,'*')
        
        % add new merged region to V_merge
        V_merge = V_merge+tempLL;
        
    end
    
    % Label regions of V_merge
    V_merge = bwlabel(V_merge);
    
    
    L_dividedD7 = V_merge;
    temp = L_D7;
    temp(V_merge==0)=0;
    
    L_dividedD7(temp==0)=0;
end

% if V_merge is empty (meaning L_D7 did not need to get divided) then funtion should just return one big region for
% entire image
if isempty(V_merge)==1
    V_merge_D7 = ones(size(L_D7));
else
    V_merge_D7 = V_merge;
end

V_merge=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D14



% find out it L_D14 need to be seperated
cmlistD14 = [cmD1;cmD7;cmD14];

S = size(cmlistD14);


        rL = regionprops(bwlabel(BWD14),'Centroid');   % get CM's of blobs in this time point
        rC = regionprops(L_dividedD7,'Centroid'); % get CM's of objects in previous time point

        % list of indexes of blops in labeled image
        listBlobsNotMerge = zeros(1,length(rC)); % length(rC) is number of clusters/objects at t-1

        % for each object at t-1 find the closest blob at t
        for k=1:length(rC)
            dist = 1e10;
            index = 0;
            for n=1:length(rL)
                
                cmL = rL(n).Centroid; % (Note 'Centroid' returns first x = col then y = row)
                cmC = rC(k).Centroid; % (Note 'Centroid' returns first x = col then y = row)
                
                temp_dist = sqrt((cmL(1)-cmC(1))^2 + (cmL(2)-cmC(2))^2 ); % distance between n'th blop at time+t and k'th object at time=t-1
                
                if dist>temp_dist && temp_dist<55% save smallest distance
                    dist = temp_dist;
                    index = n;
                end
            end
            listBlobsNotMerge(k) = index; % put index of closest blop D7 in list
        end


if  S(1)==1 % L_D14 didn't need division
    
    L_dividedD14 = L_D14;
    
    
elseif  length(listBlobsNotMerge)==length(unique(listBlobsNotMerge))  && isempty(find(listBlobsNotMerge==0))==1 % if number of unique blop indexes is same as number of object no blop has has merged!
    
    [L_dividedD14,numMetEvents ] = return_clustered_C(BWD14,return_dist_matrix(BWD14),20,listBlobsNotMerge);
    numMetEvents;
else % L_D14 needs voronoi division!!
    
    objects = []; % make data struc to hold all relevant objects/metastatic events needed for this division
    objectCounter = 0; % count the objects as they are listed
    
    % iterate through metastatic events which happened before D1
    for k=1:length(indexD1)
        
        objectCounter = objectCounter+1;
        
        % isolate object in L_D1 (where it first occured)
        L_k = return_sub_listL(L_D1,indexD1(k));
        
        % find index of this object in L_dividedD7 - we need the most recent extent of this object
        % PROBLEM: THEY MAY NOT OVERLAP!!!! AND WHEN DILATED THEY MAY OVERLAP WITH OTHER REGIONS!!!
        BWL_k = L_k; % make bw version of L_k
        BWL_k = bwmorph(BWL_k,'dilate',5); % dilate the bw version of L_k to maximize chance of overlap
        temp = L_dividedD7; % make temp version of L_dividedD7
        temp(BWL_k==0)=0; % remove everything but dilated region from day 1 (L_D1)
        
        % one way of finding the index of the object is finding the rounded of average of nonzero
        % indecis) - (hopefully the right one is more common)!
        index1 = round(mean(nonzeros(temp)));
        % another way is simply finding the max - hoping its the only one
        index2 = max(nonzeros(temp));
        
        % if the two version do not give same results or both are zero it
        % is safer just to use boundary of object at it was on day 1. (Remember we have weeded out object that disappear at later timepoint so there should be something there...)
        if index1~=index2 || index1==0 % we failed to determine correct index of object in L_dividedD7 so we use voronoi points extracted from L_D1 (using L_k)
            L_k = bwmorph(L_k,'dilate',DCD1);
            [xx,yy]= return_object_voronoi_points(L_k);
        else
            % add both points from object on D1 and D7 (otherwise we do not have enough points inside the object)
            [xxD1,yyD1]= return_object_voronoi_points(L_k);
            L_k = return_sub_listL(L_dividedD7,index1); % isolate region with index1 in L_dividedD7 and extract voronoi points
            if L_D7DividedYesNo==0
                L_k = bwmorph(L_k,'dilate',DCD7);
            end
            [xx,yy]= return_object_voronoi_points(L_k);
            xx=[xx;xxD1];
            yy=[yy;yyD1];
            
            % weed out points which are to close
            % if points are too close replace one of them with (0,0) - then later remove
            for ii=1:length(xx)
                for jj=1:length(xx)
                    dist = sqrt((xx(ii)-xx(jj))^2 +(yy(ii)-yy(jj))^2);
                    if ii~=jj && dist <10
                        xx(ii)= 0;
                        yy(ii)= 0;
                    end
                end
            end
            
            % remove zeros!
            xx = nonzeros(xx);
            yy = nonzeros(yy);
        end
        
        objects(objectCounter).voronoiPointsX = xx;
        objects(objectCounter).voronoiPointsY = yy;
        
    end
    
    % iterate through metastatic events which happened between D1 and D7
    for k=1:length(indexD7)
        
        objectCounter=objectCounter+1;
        
        L_k = return_sub_listL(L_D7,indexD7(k));
        if L_D7DividedYesNo==0
            L_k = bwmorph(L_k,'dilate',DCD7);
        end
        
        [xx,yy]= return_object_voronoi_points(L_k);
        
        objects(objectCounter).voronoiPointsX = xx;
        objects(objectCounter).voronoiPointsY = yy;
        
    end
    
    % iterate through metastatic events which happened between D7 and D14
    for k=1:length(indexD14)
        
        objectCounter=objectCounter+1;
        
        L_k = return_sub_listL(L_D14,indexD14(k));
        
        [xx,yy]= return_object_voronoi_points(L_k);
        
        objects(objectCounter).voronoiPointsX = xx;
        objects(objectCounter).voronoiPointsY = yy;
        
    end
    
    % extract all points into one vector so they can be feed to next
    % function: return_labeled_voronoi(x,y,rows,cols)
    
    x=[];
    y=[];
    for k=1:length(objects)
        x = [x; objects(k).voronoiPointsX];
        y = [y; objects(k).voronoiPointsY];
    end
    
    %%%%
    
    [rows,cols] = size(L_D14);
    
    
    
    V_labeled = return_labeled_voronoi(x,y,rows,cols);
%         
%     figure
%     imshow(V_labeled,[])
%     colormap(jet)
    
    
    V_merge = zeros(size(V_labeled));
    
    % determine which regions in V_labeled should be merged/joined
    for k=1:(length(objects)) % iterate through all objects/metastatic events 
        
        tempList=[]; % list of voronoi regions i V_labeled which should be joined
        pointsX = objects(k).voronoiPointsX;
        pointsY = objects(k).voronoiPointsY;
        for n=1:length(objects(k).voronoiPointsX) % iterate through voronoi points x,y for object k
            tempList(n)= V_labeled(round(pointsY(n)),round(pointsX(n))); % add index of region to list (remember X=cols Y=rows)
        end
        
        % make empty matrix with 0's
        tempLL = zeros(size(V_labeled));
        
        % fill in 1's in tempLL in the regions which need to be joined
        for g=1:length(tempList)
            tempLL(V_labeled==tempList(g))=1;
        end
        
        % dilate so regions merge
        tempLL = bwmorph(tempLL,'dilate',4);
        % Attempt to close holes by filling in white in edges
%         tempLL(5:rows-5,1)=1;
%         tempLL(5:rows-5,cols)=1;
%         tempLL(1,2:cols-2)=1;
%         tempLL(rows,2:cols-2)=1;
        
        % fill holes and erode
        tempLL = imfill(tempLL,'holes');
        tempLL = bwmorph(tempLL,'erode',4);
        
        % add new merged region to V_merge
        V_merge = V_merge+tempLL;

%      figure
%      imshow(tempLL,[])
%      hold on
%      plot(pointsX,pointsY,'*')
        
    end
    
%      figure
%      imshow(V_merge,[])
    
    % Label regions of V_merge
    V_merge = bwlabel(V_merge);
    
    % use regions in V_merge to make L_dividedD14
    L_dividedD14 = V_merge;
    temp = L_D14;
    temp(V_merge==0)=0;
    
    L_dividedD14(temp==0)=0;
    
    
end


if isempty(V_merge)==1
    V_merge_D14 = ones(size(L_D14));
else
    V_merge_D14 = V_merge;
end




end



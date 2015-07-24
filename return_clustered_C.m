function [C,numMetEvents ]= return_clustered_C(BW,distMatrix,threshold,doNotMergeList)
% returnClusteredC() takes black and white image with corresponding
% distmatrix and determines which objects are closer than threshold pixels
% and labels these blops with the same index in C. It does however NOT
% merge blob it their index are in the list 'doNotMergeList'


if isempty(distMatrix)==0 && isempty(BW)==0
    % distMatrix and BW was not empty so bw need to be clustered

    C  = bwlabel(BW);
    distVec = squareform(distMatrix); % returns a vector form of distmatrix (which linkage needs)
    LIN = linkage(distVec); % returns a matrix LIN that encodes a tree of hierarchical clusters of the rows of distMatrix (row/col number of object is the same as index of object in Labeled matrix)
    % LIN has size (numBlobs-1)x3. 1st and 2nd column hold index of two
    % dofferent objects and 3rd column has the pixels distance/threshold
    % which would make these two objects clustered

    
    countJoins = length(LIN(:,3)<threshold);
        
    
    % generate blobInfoMatrix 
    % 1 col: new_index
    % 2 col: old_index
    % 3 col: indexOfClosestD(t-1)Blop
    % 4 col: D(t-1)Blop YesNo
    blobInfoMatrix = zeros(max(C(:)),4);
    
    
    for i=1:max(C(:))

        blobInfoMatrix(i,1:2) = i;
        
        if isempty(nonzeros(doNotMergeList==i))==1
            blobInfoMatrix(i,4)=0; % blop is NOT on list
        else
            blobInfoMatrix(i,4)=1; % blop is on list           
        end
        
        % find closest D(t-1) blop for each blop
        smallest_dist = 1e10;
        indexOfClosestD1Blob = 0;
        for j=1:length(doNotMergeList)
            
            dist = distMatrix(i,doNotMergeList(j));%sqrt((cm_i(1)-cm_j(1))^2+(cm_i(2)-cm_j(2))^2); % dist between blop i and blop j (a D(t-1) blop)
            
            if dist<smallest_dist && dist<threshold && i~=doNotMergeList(j) && isempty(nonzeros(doNotMergeList==i))==1 
                smallest_dist = dist;
                indexOfClosestD1Blob=doNotMergeList(j);
            end
        end
        blobInfoMatrix(i,3)=indexOfClosestD1Blob;
    end
    
    
    M_merge = [];
    mergeCount = 0;
    
    for i=1:length(distMatrix)
        for j=i:length(distMatrix)        
            if distMatrix(i,j)<threshold && i~=j
                mergeCount = mergeCount+1;  
                M_merge(mergeCount,1:3) = [i j distMatrix(i,j)];
            end        
        end    
    end
    
    if isempty(M_merge)==0
    [SS I] = sort(M_merge(:,3));

    M_merge = M_merge(I,:);
        
    M_merge_copy = M_merge;

    size_Merge = size(M_merge);
    cc=0;
    numBlobs = max(C(:));
    
    for g = 1:size_Merge(1)%:-1:1
        
        
       % doNotMergeList
        index1 = M_merge(g,1);
        index2 = M_merge(g,2);
        old_index1 = M_merge_copy(g,1);
        old_index2 = M_merge_copy(g,2);
        
        if  (isempty(doNotMergeList(doNotMergeList==index1))==0 && isempty(doNotMergeList(doNotMergeList==index2))==0)
            % do nothing!!
        else

            cc=cc+1;
            
            
            % make sure blobs are not merged to clusters containing D1 blob
            % which is not the one closest to the blob
            if isempty(nonzeros(doNotMergeList==index1))==0
               % blob 1 (index1) is a D(t-1) blob
               % get old index of closets D1 blob of obhect 2
               indexOfClosestD1Blob = blobInfoMatrix(find(blobInfoMatrix(:,2)==old_index2),3); % make sure blop is merged to closest D(t-1) blob (might not be index1 blob!!)
               % use old index to find new one;
               if indexOfClosestD1Blob~=0
                   new_index2 = blobInfoMatrix(indexOfClosestD1Blob,1);
                   new_index1 = index1; 
               else 
                   new_index2 = numBlobs+cc; 
                   new_index1 = numBlobs+cc; 
               end
            elseif index2>numBlobs && isempty(nonzeros(doNotMergeList==index1))==1
               new_index2 = max(index1,index2);
               new_index1 = max(index1,index2);
            elseif isempty(nonzeros(doNotMergeList==index1))==1
               new_index2 = numBlobs+cc; 
               new_index1 = numBlobs+cc; 
            end
            
            if isempty(nonzeros(doNotMergeList==index2))==0
               % blob 2 (index2) is a D(t-1) blob
               % get old index of closets D1 blob
               indexOfClosestD1Blob = blobInfoMatrix(find(blobInfoMatrix(:,2)==old_index1),3); % make sure blop is merged to closest D(t-1) blob (might not be index1 blob!!)
               % use old index to find new one
               if indexOfClosestD1Blob~=0
                   new_index1 = blobInfoMatrix(indexOfClosestD1Blob,1);
                   new_index2 = index2;
               else 
                   new_index2 = numBlobs+cc; 
                   new_index1 = numBlobs+cc; 
               end               
            elseif index1>numBlobs && isempty(nonzeros(doNotMergeList==index2))==1
               new_index2 = max(index1,index2);
               new_index1 = max(index1,index2);
            elseif isempty(nonzeros(doNotMergeList==index2))==1
               new_index2 = numBlobs+cc; 
               new_index1 = numBlobs+cc; 
            end            

            C(C==index2) = new_index2; % replace different index in these two objects with the same number
            C(C==index1) = new_index1; % note - this number starts at numBlops + 1 so it will be different from other index already there
            
            % change both blop index in M_merge to new_index (otherwise we will fail to do remaining merges)
            M_merge(M_merge(:,1:2)==index1)= new_index1;
            M_merge(M_merge(:,1:2)==index2)= new_index2;

            blobInfoMatrix(blobInfoMatrix(:,1)==index1)=new_index1;
            blobInfoMatrix(blobInfoMatrix(:,1)==index2)=new_index2;
            
            % Change index in doNotMergeList if any of the blobs was in there
            if isempty(doNotMergeList(doNotMergeList==index1))==0 % index1 was in list
                doNotMergeList(doNotMergeList==index1) = new_index1; % put in the new index of this object in the doMotMergeList
            end
            if isempty(nonzeros(doNotMergeList==index2))==0 %  isempty(doNotMergeList(doNotMergeList==index2))==0 % index2 was in list
                doNotMergeList(doNotMergeList==index2) = new_index2; % put in the new index of this object in the doMotMergeList
            end
            
        end
    end

    % correct if some blobs have been asssigned to different D(t-1) blob
    % than the one they are closest to
    size_blobInfoMatrix = size(blobInfoMatrix); 
    for i=1:size_blobInfoMatrix(1)
        if blobInfoMatrix(i,3)~=0 % this blob was very close to D1 blob - chech that this is still the case - if not correct!
             
            indexOfClosestD1Blob = blobInfoMatrix(i,3); 
            % new index of the closest blop
            
            new_index_closest_blop = blobInfoMatrix(indexOfClosestD1Blob,1);
            new_index = blobInfoMatrix(i,1);
            
            if new_index_closest_blop ~= new_index %
                C(bwlabel(BW)==i) = new_index_closest_blop;
               % blobInfoMatrix(i,1) = new_index_closest_blop
            end
        end
    end

    numObjects = length(unique(nonzeros(C)));
    
    numMetEvents = numObjects - length(doNotMergeList);
    
   
%     We want the largest index in C to be equal to the number of
%     clustered blobs (number of objects=tumors/metastases)
%     so we need to go back and replace the large index numbers (we got above) with the
%     lower index numbers which are now freed up!

    for i=1:numObjects  % the new number of labeled clutered regions is: numBlobs-countJoins
        if isempty(find(C==i,1))==1 % if the number i is not being used as label in C then take the largest index used and replace it with i
            C(C==max(C(:))) = i;
        end    
    end
    
    else
        % M_merge was empty -  just return normal labeled matrix
        C = bwlabel(BW);
        
        numObjects = length(unique(nonzeros(C)));
    
        numMetEvents = numObjects - length(doNotMergeList);
    
        if numMetEvents<0
            numMetEvents =0;
        end
    end
    
elseif isempty(distMatrix)==1 && isempty(BW)==0 % if distMatrix is empty just return normal labeled matrix

    C = bwlabel(BW);
    
    numObjects = length(unique(nonzeros(C)));
    
    numMetEvents = numObjects - length(doNotMergeList);
    
    if numMetEvents<0
        numMetEvents =0;
    end
elseif isempty(BW)==1 % if BW is an empty matrix - return empty matrix
    
    C = [];
    
    numMetEvents = [];
     
end


end % end of function



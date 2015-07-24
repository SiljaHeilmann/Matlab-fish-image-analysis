function [cmD1,indexD1,cmD7,indexD7,cmD14,indexD14] = return_met_events(L_D1,L_D7,L_D14)
%returnMetEvents returns list with center of mass of all metastatic events which took place between two time points 
% cmD1 list of center of mass of mets/tumor which arrose between day 0 and day1 (including the injection site!!)
% cmD7 list of center of mass of mets which arrose between day 1 and day 7 
% cmD14 list of center of mass of mets which arrose between day 7 and day 14 

cD14= 0; % number of metastatic events between d7 and d14
cD7 = 0; % number of metastatic events between d1 and d7
cD1 = 0; % number of metastatic events between d0 and d1

indexD1  = [];
indexD7  = [];
indexD14 = [];

indexListD14 = nonzeros(unique(L_D14));

% go through objects on day 14
for i=1:length(indexListD14)
    
    % make image with only object i (then dilate it so change of overlap is increased)
    L_i = return_sub_listL(L_D14,indexListD14(i));
    L_i = bwmorph(L_i,'dilate',5);
    
    % make temperary version of day 7 image - se whats left when overlain with D14 object i
    temp = L_D7;
    temp(L_i==0)=0;
    
    % get list of indexes for objects at day 7 contributing to object i at day 14
    indexListD7 = nonzeros(unique(temp)); 
    
    if isempty(indexListD7)==1 % if there are no objects within the dilated object #i of D14, at day 7 it means we had a metastatic event between day 7 and day 14
        if isempty(find(indexD14==indexListD14(i),1))==1 % only add metastatic event to list if it is not there already!
            cD14 = cD14 + 1; % count metastatic event
            L_i = double(L_i);
            L_i(L_i>0)=2;
            r=regionprops(L_i,'Centroid');
            cmD14(cD14,:) = r(2).Centroid; % r(2) because you just set number>0 to 2's!!
            indexD14(cD14)=indexListD14(i);
        end
    end
    
    for j=1:length(indexListD7)
        
        L_j = return_sub_listL(L_D7,indexListD7(j));
        
        L_j = bwmorph(L_j,'dilate',5);
        
        temp = L_D1;
        
        temp(L_j==0)=0;
        
        indexListD1 = nonzeros(unique(temp));
        
        if isempty(indexListD1)==1
            if isempty(find(indexD7==indexListD7(j),1))==1 % only add metastatic event to list if it is not there already!
                cD7 = cD7 + 1;
                L_j = double(L_j);
                L_j(L_j>0)=2;
                r=regionprops(L_j,'Centroid');
                cmD7(cD7,:)=r(2).Centroid; % r(2) because you just set number>0 to 2's!!
                indexD7(cD7)=indexListD7(j);
            end
        end
        
        for k=1:length(indexListD1)
            if isempty(find(indexD1==indexListD1(k),1))==1 % only add metastatic event to list if it is not there already!
                L_k = return_sub_listL(L_D1,indexListD1(k));
                L_k = double(L_k);
                L_k(L_k>0)=2;
                r=regionprops(L_k,'Centroid');
                cD1 = cD1 + 1;
                cmD1(cD1,:)=r(2).Centroid; % r(2) because you just set number>0 to 2's!!
                indexD1(cD1)=indexListD1(k);
            end
        end
    end
end


% if no metastatic events happened make list empty!
if cD14==0
    cmD14=[];
    indexD14=[];
end
if cD7==0
    cmD7 =[];
    indexD7=[];
end
if cD1==0
    cmD1 =[];
    indexD1 =[];
end

% % some times the same cm appears twice (probably because a region on day 7 overlaps with two different dilated region on day 14... not sure) - remove!
% if isempty(cmD7)==0
%     temp1 = cmD7(:,1);
%     temp2 = cmD7(:,2);
%     cmD7 = zeros(length(unique(temp1)),2);
%     cmD7(:,1)=unique(temp1);
%     cmD7(:,2)=unique(temp2);
% end
% 
% if isempty(cmD14)==0
%     temp1 = cmD14(:,1);
%     temp2 = cmD14(:,2);
%     cmD14 = zeros(length(unique(temp1)),2);
%     cmD14(:,1)=unique(temp1);
%     cmD14(:,2)=unique(temp2);
% end

end


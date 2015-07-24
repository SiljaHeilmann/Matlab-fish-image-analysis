function [xx,yy]= return_object_voronoi_points(L)
% returnObjectVoronoiPoints takes the nonzero region found in L as returns
% sparse points from boundary and some times cm for - these points can be used
% together with other sets of object-Voronoi-points to make voronoi regions
% (using returnLabeledVoronoi)
%
% The points need to have a certain distance apart (atleast 10 pixels) so that the voronoi
% regions does not become to crowded (we need space for 1-2 pixel wide seperation lines between them)


B = bwboundaries(L);
BO=[];

% stack boundary points in one n x 2 matrix
for k=1:length(B)
    BO = [BO ; B{k}];
end

% remove many of boundary points (they are far too close and too many!!!)
temp=[];
cc=0;
for k=1:length(BO)
    if mod(k,40)==0
        cc=cc+1;
        temp(cc,:)=BO(k,:);
    end
end

BO = zeros(length(temp),2);

if isempty(temp)==0
    % interchange columns (because bwboudaries returns x and y opposite usual)
    BO(:,1) = temp(:,2);
    BO(:,2) = temp(:,1);
end

r = regionprops(L,'Centroid','Extrema');

CM=[];
EX=[];

for k=1:length(r)
    cm = r(k).Centroid;
    if isnan(cm(1))==0
        CM = [CM;cm];
    end
end
for k=1:length(r)
    ex = r(k).Extrema;
    if isnan(ex(1))==0
        EX = [EX;ex];
    end
end



% add center of masses and extrema otherwise we might loose very small mets or points on the edges of tendrils!!
BO = [BO;CM];%;EX];

S = size(BO);

% if points are too close replace one of them with (0,0) - then later remove
for ii=1:S(1)
    for jj=1:S(1)
        dist = sqrt((BO(ii,1)-BO(jj,1))^2 +(BO(ii,2)-BO(jj,2))^2);
        if ii~=jj && dist <10
            BO(ii,:)=[0 0];
        end
    end
end

% remove zeros!
temp1 = nonzeros(BO(:,1));
temp2 = nonzeros(BO(:,2));

BO = zeros(length(temp1),2);

BO(:,1) = temp1;
BO(:,2) = temp2;

xx=BO(:,1);

yy=BO(:,2);

end


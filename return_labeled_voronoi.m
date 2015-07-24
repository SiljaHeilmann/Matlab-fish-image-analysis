function V_labeled = return_labeled_voronoi(x,y,rows,cols)
% returnLabeledVoronoi takes x- and y-coordinates of N>2 points and makes a
% matrix of size (rows,cols) with the different voronoi regions (defined by
% the N points) labelled 1,2,3 ... N.
%
% each region is seperated from the other by a line of 0's of width 1-2 pixels (which is why the input points should not be too close!)
% Detailed explanation goes here


% TEST IF EVERYTHING IS IN ORDER
% test if there is enough points to do voronoi!
if length(x)<3
    display('not enough points to use voronoi - (it needs atleast 3!)')
end
% test whether some points are too close to do voronoi!
for ii=1:length(x)
    for jj=1:length(x)
        dist = sqrt((x(ii)-x(jj))^2 +(y(ii)-y(jj))^2);
        if ii~=jj && dist <10
            display('some points are too close to use returnLabeledVoronoi!)')
        end
    end
end

% MATLABs voronoi function returns the finite vertices of the Voronoi edges
% in vx and vy (coordinates of end points of lines between regions)
[vx,vy]=voronoi(x,y);

% Initialize empty matrix
V = zeros(rows,cols);

% iterate through all voronoi edges 
for h=1:length(vx)
       
    % each edge will draw a line in V - determine x range [XMIN,XMAX] this line spans
    % within V
    XMIN = max(1,min(vx(1,h),vx(2,h)));
    XMAX = min(cols,max(vx(1,h),vx(2,h)));
    
    % The MATLAB voronoi function returns arbitrary end points for edges
    % that go to infinity - sometimes these are within V - if thats the
    % case then we need to extend the x range [XMIN, XMAX] further!
    if min(x)-50>XMIN && XMIN>1
        XMIN=1;
    end
    if max(x)+50<XMAX && XMAX<cols
        XMAX=cols;
    end
    teller = vy(1,h)-vy(2,h);
    nevner = vx(1,h)-vx(2,h);
    % slope of voronoi edge
    slope = teller/nevner;
    if max(y)+50<max(vy(1,h),vy(2,h)) && max(vy(1,h),vy(2,h))<rows % the +50 is to ensure some distance from lowest voronoi point
        if slope<0
            XMIN=1;
        elseif slope>0
            XMAX=cols;
        end
    end
    if min(y)-50>min(vy(1,h),vy(2,h)) && min(vy(1,h),vy(2,h))>1 % the -50 is to ensure some distance from highest voronoi point
        if slope<0
            XMAX=cols;
        elseif slope>0
            XMIN=1;
        end
    end

    % x-coor span for drawing edge/line on V
    xx = linspace(XMIN,XMAX,20000);

    if abs(slope)>50 || slope==-Inf || slope==Inf
        xx = linspace(XMIN,XMAX,60000);
    end        
    
    % intersection with y-axis of edge/line
    b = vy(1,h)-vx(1,h).*slope;
    
    % y - values of line
    yy = xx.*slope+b;
    
    % if slope is Inf (nevner ==0, edge is vertical) we need to generate xx and yy in a
    % diferent way
    if nevner==0
        xx = xx.*0 +vx(1,h)*ones(1,length(xx));
        yy = linspace(min(vy(1,h),vy(2,h)),max(vy(1,h),vy(2,h)),20000);
        
        if abs(slope)>50 || slope==-Inf || slope==Inf
            yy = linspace(min(vy(1,h),vy(2,h)),max(vy(1,h),vy(2,h)),60000);
        end
    end
    
    % get rid of xx's and yy's that are outside V. (Set them to 1) ... some
    % of these might be overkill... hmmm
    xx(find(yy>rows))=1;
    yy(find(yy>rows))=1;
    xx(find(yy<1))=1;
    yy(find(yy<1))=1;
    xx(find(xx>cols))=1;
    yy(find(xx>cols))=1;
    xx(find(xx<1))=1;
    yy(find(xx<1))=1;
    
    % make new temp matrix with 1's (same size as V)
    temp = zeros(rows,cols);

    % draw edge/line with 1's in temp matrix of 0's   
    for i=1:length(xx)
        temp(round(yy(i)),round(xx(i)))=1;
    end
    % remove the white pixel set to 1 by the [xx,yy]=[0,0]
    temp(1,1)=0; 

    % make line fatter by dilating
    temp=bwmorph(temp,'dilate',1);
    
    % add line to V
    V = V + temp;
    
end

% set everything > 0 to 1
V(V>0)=1;
% make negative
V=V.*-1+1;

% label each region with 1,2,3,....,N      
V_labeled=bwlabel(V);

end


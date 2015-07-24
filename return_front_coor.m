function [x_coor,y_coor1,y_coor2,y_coor_base] = return_front_coor(body_bw1,body_bw2)

% find overlap of body_bw1 and body_bw2
body = body_bw1;
body(body_bw2==0)=0;

fish_length = return_fish_length(body);

% find first pixel in fish body outline
[~,firstPixel] = find(body,1,'first');

% find last pixel in fish body outline
[~,lastPixel] = find(body,1,'last');

start = firstPixel + 5; % first x-coor of boundary landmark points  
slut = firstPixel + round(fish_length*0.80);% last x-coor of boundary landmark points  

% x_coordinates for equally spaced landmark points equally spaced along the upper and
% lower fish body boundary
x_coor = round(linspace(start,slut,8));

x_coor(end-1)=[]; % skip ventral and dorsal fins! (it varies to much between right and left side images) 

% add extra points on tail (after fins)
x_coor(end) = lastPixel - 100;% slut+130; % one more point on tail boundary
x_coor(end+1) = lastPixel - 50;%slut+230; % one more point on tail boundary

% initialize
yup1_coor = zeros(1,length(x_coor));
ydown1_coor = zeros(1,length(x_coor));
yup2_coor = zeros(1,length(x_coor));
ydown2_coor = zeros(1,length(x_coor));

% find y-coordinates using fish body outline
for i=1:length(x_coor)
        yup1_coor(i) = find(body_bw1(:,x_coor(i)),1,'first');
        yup2_coor(i) = find(body_bw2(:,x_coor(i)),1,'first');
        ydown1_coor(i) = find(body_bw1(:,x_coor(i)),1,'last');
        ydown2_coor(i) = find(body_bw2(:,x_coor(i)),1,'last');        
end

% assign output values:
x_coor = [x_coor x_coor];
y_coor1 = [yup1_coor ydown1_coor];
y_coor2 = [yup2_coor ydown2_coor];
y_coor_base = round((y_coor1 + y_coor2)./2);
 
end


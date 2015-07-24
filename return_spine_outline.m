function spine_bw = return_spine_outline(bf_red,body)
% This function returns an outline/mask of the fish spine in the tail
% using a bf image of the fish and the body boundary mask


% find last pixel in fish body outline
[~,lastPixel] = find(body,1,'last');

% new name for bf_red
im = bf_red;
        
% we wish to cut out an image of fish tail with spine and muscle (but no background plate)
% at a distance temp_width (150) back from the lastPixel. This cutout will be used to
% find an optimal threshold for distinguishing spine from muscle inside
% fish body
temp_width = 150;

% find upper and lower y values of rectancle for image crop
y_corner     = max([find(body(:,round(lastPixel-temp_width)),1,'first') ,find(body(:,round(lastPixel-temp_width)+100),1,'first')]);
y_corner_low = min([find(body(:,round(lastPixel-temp_width)),1,'last') ,find(body(:,round(lastPixel-temp_width)+100),1,'last')]);

% rectangle for cropping image
spine_rect = [round(lastPixel-temp_width) y_corner+5 +100 y_corner_low-y_corner-10];

crop_bf = imcrop(im,spine_rect);

% level for threshold (using Otsu'z method) 
level = graythresh(crop_bf);
 
im = im2bw(im,level); 

% negative
im = im.*-1 +1;

% erode body and use this as mask to remove unwanted objects around spine
% outline
body = bwmorph(body,'erode',15);%6
im(body==0)=0;

% erode spine outline to get rid of small protrusions (blood vessel ect.)
im = bwmorph(im,'erode',5);

% remove objects smaller than 1000 pixels (noise/dirt)
thres = 1000;
im = bwareaopen(im,thres);    

% dilate to smoothen
im = bwmorph(im,'dilate',4);
                 
% fill in holes
im = imfill(im,'holes');

% smoothen outline
im = smoothBW(im,10);  

% assign to output value
spine_bw =im;
        
end
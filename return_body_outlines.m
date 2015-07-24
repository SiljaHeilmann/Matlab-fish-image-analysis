function [body1,body2] = return_body_outlines(bf1,R_body1,bf2,R_body2)

% erode body (find 'core' of tail - get rid of tail fin)
R_body1 = bwmorph(R_body1,'erode',20);
R_body1 = bwareaopen(R_body1,20000);
R_body2 = bwmorph(R_body2,'erode',20);
R_body2 = bwareaopen(R_body2,20000);

% find nose and tail coor of both left and right side images
[y,x] = find(R_body1,10,'first');
nose_coor1 = [round(mean(x)),round(mean(y))];
[y,x] = find(R_body1,1000,'last');
tail_coor1 = [round(mean(x)),round(mean(y))];
fish_length1 = tail_coor1(1)-nose_coor1(1);

[y,x] = find(R_body2,10,'first');
nose_coor2 = [round(mean(x)),round(mean(y))];
[y,x] = find(R_body2,1000,'last');
tail_coor2 = [round(mean(x)),round(mean(y))];
fish_length2 = tail_coor2(1)-nose_coor2(1);

% calculate mean fish length
fish_length = min(tail_coor1(1),tail_coor2(1)) - mean(nose_coor1(1),nose_coor2(1));

% rectangle for looking at only end of tail (this is where threshold for finding boundary is hardest to determine)
rect_tail = [10 + fish_length.*0.95  mean(tail_coor1(2),tail_coor2(2))-100 200 200];

% cut out end of tail in red channel of rgb image - use this image
% to find theshold
tail_im1 = imcrop(bf1(:,:,1),rect_tail);
tail_im2 = imcrop(bf2(:,:,1),rect_tail);

% find theshold using graythres (uses Otsu's method, we need to lower
%it a bit to avoid holes inside fish body so we multiply with 0.9)
T1 = graythresh(tail_im1).*0.9;
T2 = graythresh(tail_im2).*0.9;

% Make bw images using fish specific thresholds
bw1 = im2bw(bf1(:,:,1),T1).*-1+ ones(size(bf1(:,:,1)));
bw2 = im2bw(bf2(:,:,1),T2).*-1+ ones(size(bf2(:,:,1)));

% dilate rough body outline and use as mask for new body outline (remove dirt and shadows around fish)
R_body1 = bwmorph(R_body1,'dilate',40);
R_body1 = bwareaopen(R_body1,20000);
R_body2 = bwmorph(R_body2,'dilate',40);
R_body2 = bwareaopen(R_body2,20000);
bw1(R_body1==0)=0;
bw2(R_body2==0)=0;
bw1(:,round(10 + fish_length.*1.3):end)=0;
bw2(:,round(10 + fish_length.*1.3):end)=0;

% fill hole, erode and dilate to get rid of small protrusions
% (smooths boundary)
bw1 =imfill(bw1,'holes');
bw2 =imfill(bw2,'holes');
bw1 = bwmorph(bw1,'dilate',3);
bw2 = bwmorph(bw2,'dilate',3);
bw1 =imfill(bw1,'holes');
bw2 =imfill(bw2,'holes');
bw1 = bwmorph(bw1,'erode',10);
bw2 = bwmorph(bw2,'erode',10);
bw1 = bwmorph(bw1,'dilate',7);
bw2 = bwmorph(bw2,'dilate',7);

% remove regions of 1's smaller than 20000 (not part of fish body)
bw1 = bwareaopen(bw1,20000);
bw2 = bwareaopen(bw2,20000);

% make boundary of fish body smoother
bw1 = smoothBW(bw1,30);
bw2 = smoothBW(bw2,30);

body1 = bw1;
body2 = bw2;
end
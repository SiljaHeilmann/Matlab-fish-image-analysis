function bw = return_side(flourescentImage,body,eye)
% flourescentImage is type uint16

eye = logical(eye);

% make into type double with value between 0 and 1
fp = double(flourescentImage)./2^16;

% thresholds for different sehmentation processes
t_high = 0.1190;% threshold used to make pessimist
t_low  = 0.017;% threshold used to make optimist 
t_edge = 0.045;% threshold used on g (edge finder)

% initialize segmented images
pessimist = zeros(size(fp));
optimist  = zeros(size(fp));

% pessimist is a thresholded image with only very bright gpf regions
pessimist(fp>t_high)  = 1;
pessimist(fp<=t_high) = 0;
pessimist(body==0)=0;
pessimist(eye==1)=0;
pessimist=bwmorph(pessimist,'dilate',3);
pessimist=imfill(pessimist,'holes');
pessimist=bwmorph(pessimist,'erode',2);

% optimist is a thresholded image with even very dim gpf regions
optimist(fp>t_low)    = 1;
optimist(fp<=t_low)   = 0;
optimist(body==0)=0;
optimist(eye==1)=0;
optimist = bwareaopen(optimist,3);

% Edge finder g
% return_gradmag_gray gives an grayscale image which shows the magnitude of
% the local pixel intensity gradient 
g = return_gradmag_gray(fp);

% use an eroded body mask - return_gradmag_gray tends to find artifacts
% near fish edge so we want to avoid regions close to edge
body_g = bwmorph(body,'erode',8);

% make segmented image using gradient grayscale image
bw_g = g; % initialize
bw_g(g>t_edge)    = 1;
bw_g(g<=t_edge)   = 0;
bw_g(body_g==0)=0;
bw_g(eye==1)=0;
bw_g = bwmorph(bw_g,'dilate',2);
bw_g = imfill(bw_g,'holes');
bw_g = bwmorph(bw_g,'erode',2);
bw_g = bwareaopen(bw_g,10);

% Varying threshold
% make segmented image using a varying threshold base on local average pixel intensity 
delta = 50;% width of pixwl column for calculating average local intensity inside fish body
var = return_varying_thresh(fp,delta,body);
var(body_g==0)=0;
var(eye==1)=0;
var = bwareaopen(var,2);

% merge var and bw_g
bw = var;
bw(bw_g==1)=1;

% fill holes
bw = imfill(bw,'holes');

% set regions not found by optimist to 0
bw(optimist==0)=0;
% set regions found by pessimist to 1
bw(pessimist==1)=1;

% remove 1's inside eye and outside body
bw(eye==1)=0;
bw(body==0)=0;

% fill holes
bw = imfill(bw,'holes');

end
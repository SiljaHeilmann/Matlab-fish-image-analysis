function [BW] = return_varying_thresh(fp,delta,body)

% find width of image
[~,c]=size(fp);

fp(body==0) = 0; % set intensity outside fish body to zero so that it doesnt contribute to varying threshold 

% initialize
meanOfFP = zeros(1,c);
temp     = zeros(1,c);
signal   = zeros(1,c);

% iterate through image column by column and find average intensity
for nn=1:c
    signal(nn) = mean(nonzeros(fp(:,nn)));
end

% set first delta/2 places of meanOfFP equal to signal (alway in tip of fish nose so does not matter)
meanOfFP(:,1:delta/2) = signal(:,1:delta/2);

% width delta columns contribute to average
for nn=1+delta/2:c-delta/2
   meanOfFP(nn) = mean(mean(nonzeros(fp(:,nn-delta/2:nn+delta/2))));
end

BW = fp;%initialize BW

for nn=1:c
    % calculate varying threshold for each column in image 
    temp(nn) = (meanOfFP(nn)^0.6)*0.29;
    BW(:,nn) = fp(:,nn)>temp(nn); % row for row create BW image with varying thresh 'temp'
end

end
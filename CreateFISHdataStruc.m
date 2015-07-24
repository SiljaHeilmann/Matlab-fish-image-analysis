
%% Add original images to structure array "FISH"
%  As the data structure FISH is created each fish get a unique ID contaning
%  information about  original fish number, batch, implant size and implant
%  location. Original images are then imported from their folders into matlab and
%  saved in FISH  

clear all
close all


% Path to folder with original images and all matlab files - make sure
% that current matlab folder is 'Matlab_files_and_sample_images' 
baseDir = pwd; % pwd returns the current folder as a string

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set show to 0 if you DO NOT wish to visually inspect the
% progress
show = 1;

if show==1
   display('show = 1 so progress will be shown in figure(1). Set show=0 for faster performance and no visual output.') 
   display(' ' )
else
   display('show = 0 so progress will not be shown. Set show=1 for visual output.')
   display(' ' )
end
% having show = 1 will slow down the process significantly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name_for_final_data_struct = 'FISH_data_struct.mat'; % used at the very end (line ~1561) when saving final data struct with segmented and transformed image in current directory

display('Step 1: Initializing the data structure FISH and reading in original images...')

tic % start timer


timePoint = {'1' '2' '3'}; %names of time point folders

batch    = {'BATCH1'};% 'BATCH2'}; %names of batch folders

implantSize = {'1x10^6'};% '5x10^5' '1x10^5'}; % first half of group name folder names

location = {'VENTRAL'};% 'DORSAL'}; % second half of group name folder names

fish = 0; % counter for counting fish as they are added to data structure

%Initiale data structure
FISH(1).fishID = [];

for bb = 1:length(batch)
    for tt = 1:length(timePoint)
        for ii = 1:length(implantSize)
            for ll = 1:length(location)
                
                % make list of fish folders
                fileListFISH  = dir([baseDir '/'  batch{bb}  '/' timePoint{tt} '/' implantSize{ii} '' location{ll} '/*']);
                
                for f = 1:length(fileListFISH) % iterate through fish folders
                    if strcmp( fileListFISH(f).name,'.')==0 && strcmp( fileListFISH(f).name,'..')==0 && strcmp( fileListFISH(f).name,'.DS_Store')==0 % hidden files we are not interested in
                        
                        % create unique fish ID string
                        fishID = ['fish - ' fileListFISH(f).name ' - ' batch{bb} ' - ' implantSize{ii} ' - ' location{ll}];
                        
                        % return the index of the fish if it is already in 'FISH'
                        % (fish index will be empty if it is not already in
                        % data struct
                        fishIndex = find(strcmp({FISH(:).fishID}, fishID));
                        
                        % add fish to data struct if not seen before
                        if isempty(fishIndex) % fish was not in database
                            
                            fish = fish +1; % increase every time a new fish is added
                            fishIndex = fish;
                            
                            % save info about new fish
                            FISH(fishIndex).fishID       = fishID;
                            FISH(fishIndex).batch        = batch{bb};
                            FISH(fishIndex).implantSize  = implantSize{ii};
                            FISH(fishIndex).location     = location{ii};
                            
                        end
                        
                        % make list of image files inside fish folder
                        fileListIM = dir([baseDir '/'  batch{bb}  '/' timePoint{tt} '/' implantSize{ii} '' location{ll} '/' fileListFISH(f).name '/Snap-*']);
                        
                        imCount = 0; % count images from fish folder as they are added
                        
                        for ff = 1:length(fileListIM) % iterate through im files
                            imCount = imCount + 1; % order in folder determine what type image it is
                            % import images and save in FISH
                            if show==1
                                figure(1)
                            end
                            
                            switch imCount
                                case 1 % BF nose right
                                    FISH(fishIndex).times(tt).originals.bf1  = imread([baseDir '/'  batch{bb}  '/' timePoint{tt} '/' implantSize{ii} '' location{ll}  '/' fileListFISH(f).name '/' fileListIM(ff).name]);
                                    if show==1
                                    subplot(2,3,1)
                                    imshow(FISH(fishIndex).times(tt).originals.bf1,[])
                                    hold on
                                    title(['Right side. Fish number ' num2str(fishIndex) ',  time: ' num2str(tt)],'FontSize',15)
                                    end    
                                case 2 % GFP nose right
                                    FISH(fishIndex).times(tt).originals.gfp1 = imread([baseDir '/'  batch{bb}  '/' timePoint{tt} '/' implantSize{ii} '' location{ll}  '/' fileListFISH(f).name '/' fileListIM(ff).name]);
                                    if show==1
                                    subplot(2,3,2)
                                    imshow(FISH(fishIndex).times(tt).originals.gfp1,[])
                                    hold on
                                    title(['GFP'],'FontSize',15)
                                    end
                                case 3 % RFP nose right
                                    FISH(fishIndex).times(tt).originals.rfp1 = imread([baseDir '/'  batch{bb}  '/' timePoint{tt} '/' implantSize{ii} '' location{ll}  '/' fileListFISH(f).name '/' fileListIM(ff).name]);
                                    if show==1
                                    subplot(2,3,3)
                                    imshow(FISH(fishIndex).times(tt).originals.rfp1,[])
                                    hold on
                                    title(['RFP'],'FontSize',15)
                                    end
                                case 4 % RFP nose left
                                    FISH(fishIndex).times(tt).originals.rfp2 = imread([baseDir '/'  batch{bb}  '/' timePoint{tt} '/' implantSize{ii} '' location{ll}  '/' fileListFISH(f).name '/' fileListIM(ff).name]);
                                    if show==1
                                    subplot(2,3,6)
                                    imshow(FISH(fishIndex).times(tt).originals.rfp2,[])
                                    hold on
                                    title(['RFP'],'FontSize',15)
                                    end
                                case 5 % GFP nose left
                                    FISH(fishIndex).times(tt).originals.gfp2 = imread([baseDir '/'  batch{bb}  '/' timePoint{tt} '/' implantSize{ii} '' location{ll}  '/' fileListFISH(f).name '/' fileListIM(ff).name]);
                                    if show==1
                                    subplot(2,3,5)
                                    imshow(FISH(fishIndex).times(tt).originals.gfp2,[])
                                    hold on
                                    title(['GFP'],'FontSize',15)
                                    end
                                case 6 % BF nose left
                                    FISH(fishIndex).times(tt).originals.bf2  = imread([baseDir '/'  batch{bb}  '/' timePoint{tt} '/' implantSize{ii} '' location{ll}  '/' fileListFISH(f).name '/' fileListIM(ff).name]);
                                    if show==1
                                    subplot(2,3,4)
                                    imshow(FISH(fishIndex).times(tt).originals.bf2,[])
                                    hold on
                                    title(['Left side. Fish number ' num2str(fishIndex) ',  time: ' num2str(tt)],'FontSize',15)
                                    end
                            end
                            
                        end
                    end
                end
            end
        end
    end
end
toc
display(' ')
%%  Rescale overall intensity of images images based on cut out of background (BG)
%   and save intensity rescaled images in rounds(1)
%   Note: if all images are taken with the exact same settings and ligthing
%   this step should not be nessesary and may be skipped!
close all

display('Step 2: Rescaling brightness of BF, GFP and RFP images based on average background intensity...')

tic

rect = [400 1 800 110]; % [xmin ymax width height] of background cutout

% Target background intensities
meanBGintensityBF = 3.27e4;
meanBGintensityGFP = 2000;
meanBGintensityRFP = 2000;

for ff = 1:length(FISH)
    for t=1:3
        LR = 0;
        while LR<2
            LR = LR+1;
            switch LR
                % retrieave images from data struct.
                case 1                    
                    bf  = FISH(ff).times(t).originals.bf1;
                    gfp = FISH(ff).times(t).originals.gfp1;
                    rfp = FISH(ff).times(t).originals.rfp1;
                    TITLE = 'Right side';
                case 2
                    bf  = FISH(ff).times(t).originals.bf2;
                    gfp = FISH(ff).times(t).originals.gfp2;
                    rfp = FISH(ff).times(t).originals.rfp2;
                    TITLE = 'Left side';
            end
            
            % cut out piece of BG at top edge of image (contains no fish)
            cropBF = imcrop(bf(:,:,1),rect);
            cropGFP = imcrop(gfp,rect);
            cropRFP = imcrop(rfp,rect);
            
            % determine mean and std of BG in BF, GFP and RFP images
            MBF    = mean(cropBF(:));
            STDBF  = std(double(cropBF(:)));
            MGFP   = mean(cropGFP(:));
            STDGFP = std(double(cropGFP(:)));
            MRFP   = mean(cropRFP(:));
            STDRFP = std(double(cropRFP(:)));
            
            % get rid of extreme values in the cropped images - example dark
            % spot due to random dirt.
            temp_cropBF = cropBF;
            temp_cropBF(cropBF>MBF+2*STDBF) = MBF;
            temp_cropBF(cropBF<MBF-2*STDBF) = MBF;
            temp_cropGFP = cropGFP;
            temp_cropGFP(cropGFP>MGFP+2*STDGFP) = MGFP;
            temp_cropGFP(cropGFP<MGFP-2*STDGFP) = MGFP;
            temp_cropRFP = cropRFP;
            temp_cropRFP(cropRFP>MRFP+2*STDRFP) = MRFP;
            temp_cropRFP(cropRFP<MRFP-2*STDRFP) = MRFP;
            
            % calculate mean BG intensity without extreme values caused by random dirt
            MBF = mean(temp_cropBF(:));
            MGFP = mean(temp_cropGFP(:));
            MRFP = mean(temp_cropRFP(:));
            
            % correct brightness level of BF images and save in data struct FISH in rounds(1)
            fBF  = meanBGintensityBF/MBF;
            fGFP = meanBGintensityGFP/MGFP;
            fRFP = meanBGintensityRFP/MRFP;
            
            bf_new = fBF*bf;
            gfp_new = fGFP*gfp;
            rfp_new = fRFP*rfp;
            
            % save rescaled images in rounds(1)
            switch LR
                case 1
                    FISH(ff).times(t).rounds(1).bf1  = uint16(bf_new);
                    FISH(ff).times(t).rounds(1).gfp1 = uint16(gfp_new);
                    FISH(ff).times(t).rounds(1).rfp1 = uint16(rfp_new);    
                case 2
                    FISH(ff).times(t).rounds(1).bf2  = uint16(bf_new);
                    FISH(ff).times(t).rounds(1).gfp2 = uint16(gfp_new);
                    FISH(ff).times(t).rounds(1).rfp2 = uint16(rfp_new);
            end
            
            if show==1
                %visually inspect intensity normalisation/rescaling of brigthness
                figure(1)
                subplot(3,2,1)
                imshow(bf.*1.1)
                hold on
                title([ TITLE ' Original BF. Fish number ' num2str(ff) ',  time: ' num2str(t)],'FontSize',13)
                subplot(3,2,2)
                imshow(bf_new.*1.1)
                hold on
                title(['Rescaled BF. Fish number ' num2str(ff) ',  time: ' num2str(t)],'FontSize',13)
                subplot(3,2,3)
                imshow(gfp.*1.1,[0 10000])
                hold on
                title(['Original GFP. Fish number ' num2str(ff) ',  time: ' num2str(t)],'FontSize',13)
                subplot(3,2,4)
                imshow(gfp_new.*1.1,[0 10000])
                hold on
                title(['Rescaled GFP. Fish number ' num2str(ff) ',  time: ' num2str(t)],'FontSize',13)
                subplot(3,2,5)
                imshow(rfp.*1.1,[0 10000])
                hold on
                title(['Original RFP. Fish number ' num2str(ff) ',  time: ' num2str(t)],'FontSize',13)
                subplot(3,2,6)
                imshow(rfp_new.*1.1,[0 10000])
                hold on
                title(['Rescaled RFP. Fish number ' num2str(ff) ',  time: ' num2str(t)],'FontSize',13)
                
            end
        end
    end
end

toc; display(' ')
%% Flip images of fish which has nose to the right and save in rounds(1)

display('Step 3: Flipping images of fish which has nose to the right and save in rounds(1)...');

tic

close all

for ff=1:length(FISH)
    for t=1:3
        
        % retrieave images from data struct.
        gfp1 = FISH(ff).times(t).rounds(1).gfp1;
        rfp1 = FISH(ff).times(t).rounds(1).rfp1;
        bf1  = FISH(ff).times(t).rounds(1).bf1;
        gfp2 = FISH(ff).times(t).rounds(1).gfp2;
        rfp2 = FISH(ff).times(t).rounds(1).rfp2;
        bf2  = FISH(ff).times(t).rounds(1).bf2;
        
        % flip red green and blue channels of bf1 image one at a time
        temp1 = bf1;
        temp1(:,:,1) = fliplr(bf1(:,:,1));
        temp1(:,:,2) = fliplr(bf1(:,:,2));
        temp1(:,:,3) = fliplr(bf1(:,:,3));
        
        % save in rounds(1)
        FISH(ff).times(t).rounds(1).gfp1 = fliplr(gfp1);
        FISH(ff).times(t).rounds(1).rfp1 = fliplr(rfp1);
        FISH(ff).times(t).rounds(1).bf1  = temp1;
        FISH(ff).times(t).rounds(1).gfp2 = gfp2;
        FISH(ff).times(t).rounds(1).rfp2 = rfp2;
        FISH(ff).times(t).rounds(1).bf2  = bf2;
        
        if show==1
        %visually inspect that all fish face in right direction!
        figure(1)
        imshowpair(imresize(bf2,0.5),imresize(temp1,0.5))
        hold on
        title(['Flipping images where fish nose where facing right. Green right side, pink left side. (Both noses should be facing left!). Fish number: ' num2str(ff) ',  time: ' num2str(t)],'FontSize',11)
        end
    end
end

toc; display(' ')
%% Find approx. outline of fish bodies
% this rough outline/boundary of fish body will be used later for rotating the
% fish, translating the 'nose' to a specific commen point and for cropping
% the images, and for finding a more accurate outline later. 
display('Step 4: Finding approx outline of fish bodies and save in rounds(1)...');

   
tic 
close all
 
BW_threshold = 0.38; % for thresholding image (change to fit to your specific data set)
min_conn_region_size = 40000; % this number needs to be slightly less than the size of the smallest fish in the data set

% width and height of mask around fish
ww = 900; % slightly longer than length of longest fish
hh = 250; % slightly wider than width of fattest/widest (or most not horizontal) fish

for ff=1:length(FISH)
    for t=1:3
        
        % retrieave images from data struct.
        bf1  = FISH(ff).times(t).rounds(1).bf1;
        bf2  = FISH(ff).times(t).rounds(1).bf2;
        
        % merge rgb channels into on grayscale image
        bw1 = rgb2gray(bf1);
        bw2 = rgb2gray(bf2);
        
        % make black/white version of image by thresholding
        bw1 = im2bw(bw1,BW_threshold);
        bw2 = im2bw(bw2,BW_threshold);
        
        % make 0's into 1's and 1's into 0's
        bw1 = bw1.*-1 + 1;
        bw2 = bw2.*-1 + 1;
        
        % fill in holes in connected regions of 1's
        bw1 =imfill(bw1,'holes');
        bw2 =imfill(bw2,'holes');
        
        % remove connected regions of 1's with fewer than 40000 pixels
        bw1 = bwareaopen(bw1,min_conn_region_size);
        bw2 = bwareaopen(bw2,min_conn_region_size);

        % find coordinates of fish 'nose' on left and right side images
        [y1,x1] = find(bw1,5,'first');
        noseCoor1 = [round(mean(x1)) round(mean(y1))];
        [y2,x2] = find(bw2,5,'first');
        noseCoor2 = [round(mean(x2)) round(mean(y2))];
        
        % remove 1's in regions beyond the size of mask
        bw1(1:noseCoor1(2)-hh,:) = 0;
        bw2(1:noseCoor2(2)-hh,:) = 0;
        bw1(noseCoor1(2)+hh:end,:) = 0;
        bw2(noseCoor2(2)+hh:end,:) = 0;
        bw1(:,1:noseCoor1(1)-10,:) = 0;
        bw2(:,1:noseCoor2(1)-10,:) = 0;
        bw1(:,noseCoor1(1)+ww:end) = 0;
        bw2(:,noseCoor2(1)+ww:end) = 0;
        
        % remove connected regions of 1's with fewer than 40000 pixels
        bw1 = bwareaopen(bw1,min_conn_region_size);
        bw2 = bwareaopen(bw2,min_conn_region_size);
             
        % make boundary of fish body smoother
        bw1 = smoothBW(bw1,10);
        bw2 = smoothBW(bw2,10);  
        
        if show==1
        % visually inspect that approx. fish outlines look ok
        figure(1)
        imshowpair(imresize(bw1,0.5),imresize(bw2,0.5))
        hold on
        plot(noseCoor1(1).*0.5,noseCoor1(2).*0.5,'*','MarkerSize',15)
        plot(noseCoor2(1).*0.5,noseCoor2(2).*0.5,'*','MarkerSize',15)
        title(['Finding approx outline of fish bodies. Green right side, pink left side.  Fish number: ' num2str(ff) ',  time: ' num2str(t)],'FontSize',13)
        end
        
        % save approx body mask in rounds(1)
        FISH(ff).times(t).rounds(1).body1 = logical(bw1);
        FISH(ff).times(t).rounds(1).body2 = logical(bw2);

    end
end
toc; display(' ')
%% Rotate fish body to horizontal, trabslate such that nose is at fixed coor and crop images
% using rough fish body outline found above
% save rotated and cropped images in rounds(2)
display('Step 5: Rotating translating and cropping all images and save in rounds(2)...');

tic

close all

nose_coor_fix = [10,200];

ww = 1000; % width of cropped images
hh = 350;  % heigth of cropped images


for ff=1:length(FISH)
    for t=1:3
        LR = 0;
        while LR<2
            LR = LR+1;
            switch LR
                % retrieave images from data struct.
                case 1
                    bf   = FISH(ff).times(t).rounds(1).bf1;
                    gfp  = FISH(ff).times(t).rounds(1).gfp1;
                    rfp  = FISH(ff).times(t).rounds(1).rfp1;
                    body = FISH(ff).times(t).rounds(1).body1;
                case 2
                    bf   = FISH(ff).times(t).rounds(1).bf2;
                    gfp  = FISH(ff).times(t).rounds(1).gfp2;
                    rfp  = FISH(ff).times(t).rounds(1).rfp2;
                    body = FISH(ff).times(t).rounds(1).body2;
            end
            
            % find angle of fish body with respect to horizontal
            ang = find_rot_angle(body);
            
            % rotate images so that fish becomes horizontal
            bodyN = logical(imrotate(body,ang,'nearest','crop'));
            bfN   = imrotate(bf,ang,'bicubic','crop');
            gfpN  = imrotate(gfp,ang,'bicubic','crop');
            rfpN  = imrotate(rfp,ang,'bicubic','crop');
            
            % find nose coor
            [y,x] = find(bodyN,10,'first');
            nose_coor = [mean(x),mean(y)];
            
            % rectangle for cropping [upper left corner coor , width height] -
            % makes sure nose is at (10,hh/2)
            rect = [nose_coor(1)-nose_coor_fix(1) nose_coor(2)-hh/2 ww hh];%
            
            % output images (rotated and cropped)
            body_temp = imcrop(bodyN,rect);
            bf_temp   = imcrop(bfN,rect);
            gfp_temp  = imcrop(gfpN,rect);
            rfp_temp  = imcrop(rfpN,rect);
            
            if show == 1
            % visually inspect that fish where rotated and images cropped
            % correctly:
            figure(1)
            subplot(2,1,LR); imshow(bf_temp.*2,[0 4e4])
            hold on
            title(['Rotate and crop images.  Fish number: ' num2str(ff) ',  time: ' num2str(t)],'FontSize',15)
            end
            
            switch LR
                case 1
                    FISH(ff).times(t).rounds(2).bf1   = uint16(bf_temp);
                    FISH(ff).times(t).rounds(2).gfp1  = uint16(gfp_temp);
                    FISH(ff).times(t).rounds(2).rfp1  = uint16(rfp_temp);
                    FISH(ff).times(t).rounds(2).R_body1 = logical(body_temp);                    
                case 2
                    FISH(ff).times(t).rounds(2).bf2   = uint16(bf_temp);
                    FISH(ff).times(t).rounds(2).gfp2  = uint16(gfp_temp);
                    FISH(ff).times(t).rounds(2).rfp2  = uint16(rfp_temp);
                    FISH(ff).times(t).rounds(2).R_body2 = logical(body_temp);                    
            end
        end
    end
end

toc; display(' ')
%% Use the rough mask/outline (R_body) of fish body from earlier to find better more accurate outline (body)
%  save in rounds(2)

display('Step 6: Finding more accurate outline of fish bodies and save in rounds(2)...');

tic

close all

for ff=1:length(FISH)
    for t=1:3
        % retrieave images from data struct.
        bf1   = FISH(ff).times(t).rounds(2).bf1;
        R_body1 = FISH(ff).times(t).rounds(2).R_body1;        
        bf2   = FISH(ff).times(t).rounds(2).bf2;
        R_body2 = FISH(ff).times(t).rounds(2).R_body2;
        
        % the function 'return_body_outlines' uses old rough outline
        % R_body and the bf images of both sides 
        % to find new more accurate boundaries of 
        % the fish bodies
        [body1,body2] = return_body_outlines(bf1,R_body1,bf2,R_body2);
        
        % determine points along boundary for plotting
        B1 = bwboundaries(body1);
        B2 = bwboundaries(body2);
      
        if show==1
        % visually inspect the new improved outlines plotted on top of BF images
        figure(1)
        subplot(2,1,1)
        imshow(bf1.*2)
        hold on
        for k=1:length(B1)
            b=B1{k};
            plot(b(:,2),b(:,1),'r-','LineWidth',2)
        end
        title(['Find accurate body outline,  Fish number: ' num2str(ff) ' ,  Time: ' num2str(t) '  - Right side' ],'FontSize',15)
        subplot(2,1,2)
        imshow(bf2.*2)
        hold on
        for k=1:length(B2)
            b=B2{k};
            plot(b(:,2),b(:,1),'r-','LineWidth',2)
        end        
        title(['Fish number: ' num2str(ff) '  ,  Time: ' num2str(t) '  - Left side' ],'FontSize',15)
        end
        
        % Save improved body outline in rounds(2)
        FISH(ff).times(t).rounds(2).body1 = body1;
        FISH(ff).times(t).rounds(2).body2 = body2;

    end
end
toc; display(' ')
%% Find spine outline and end of spine landmark point
display('Step 7: Finding spine outlines and end of spine landmark point and save in rounds(2)...');

tic

close all

for ff=1:length(FISH)
    for t=1:3
        % retrieave images from data struct.
        bf1   = FISH(ff).times(t).rounds(2).bf1;
        body1 = FISH(ff).times(t).rounds(2).body1;
        bf2   = FISH(ff).times(t).rounds(2).bf2;
        body2 = FISH(ff).times(t).rounds(2).body2;
        
        % function 'return_spine_outline' takes red channel of bf
        % and body outline/mask and return spine outline/mask
        spine1 = return_spine_outline(bf1(:,:,1),body1);
        spine2 = return_spine_outline(bf2(:,:,1),body2);
        
        % find last pixel in outline of both spines
        [~,xx1] = find(spine1,1,'last');
        [~,xx2] = find(spine2,1,'last');
        
        % find mean x-coor of spine ends
        xx = round(mean(xx1,xx2));
        spine1(:,xx:end) = 0;
        spine2(:,xx:end) = 0;
        
        % shorten spines if they are longer than average
        if mean(xx1)>xx
            spine1 = smoothBW(spine1,10);
        elseif mean(xx2)>xx
            spine2 = smoothBW(spine2,10);
        end
        
        % find end of spine landmark points using both spine outlines
        % use mean x-coor for both and mean(yy1) and mean(yy2) for y-coor
        [yy1,xx1] = find(spine1,50,'last');
        [yy2,xx2] = find(spine2,50,'last');
        
        if show==1
        % plot spine outlines and
        % end of spine landmarks to visually inspect that they look right
        S1 = bwboundaries(spine1);
        S2 = bwboundaries(spine2); 
        figure(1)
        imshowpair(bf1,bf2)
        hold on
        plot(xx,mean(yy1),'bo','MarkerSize',15,'LineWidth',2,'MarkerFaceColor','y')
        plot(xx,mean(yy2),'ro','MarkerSize',15,'LineWidth',2,'MarkerFaceColor','y')
        for k=1:length(S1)
            b=S1{k};
            plot(b(:,2),b(:,1),'b-','LineWidth',2)
        end
        for k=1:length(S2)
            b=S2{k};
            plot(b(:,2),b(:,1),'r-','LineWidth',2)
        end
        title(['Finding spine outline and end of spine landmark   -   Fish number: ' num2str(ff) ',  Time: ' num2str(t) ' , Rigth side blue, Left side red'],'FontSize',15)
        end
                
        % save spine outlines in rounds(2)
        FISH(ff).times(t).rounds(2).spine1 = spine1;
        FISH(ff).times(t).rounds(2).spine2 = spine2;

        % save end of spine landmark coor in rounds(2)
        FISH(ff).times(t).rounds(2).spineEnd1_coor = [xx ,mean(yy1)];
        FISH(ff).times(t).rounds(2).spineEnd2_coor = [xx ,mean(yy2)];
        
    end
end
toc; display(' ')

%%  Find and save eye oulines and center of mass (landmark) in rounds(2)

display('Step 8: Finding eye outlines and eye center of mass landmark and save in rounds(2)...');

tic
close all

% put in approx coordinates of eye and approx diameter of eye
approx_eye_coor = [75 180];
approx_eye_diameter = 60;
x_shift = round(approx_eye_diameter/3);
y_shift = round(approx_eye_diameter/4);

% rectangle for cropping out part of eye, used by function
% 'return_eye_outline' inside loop
rect = [approx_eye_coor(1)-x_shift     approx_eye_coor(2)-y_shift     approx_eye_diameter*5/3     approx_eye_diameter/2]; % [xmax ymax width height]

for ff=1:length(FISH)
    for t=1:3
        LR = 0;
        
        while LR<2
            LR = LR+1;
            % retrieave images from data struct.
            switch LR
                case 1
                    bf  = FISH(ff).times(t).rounds(2).bf1;
                case 2
                    bf  = FISH(ff).times(t).rounds(2).bf2;
            end
            
            % function 'return_eye_outline' takes red channel of rgb image
            % and rectangle (based on approx. eye center and diameter) and return
            % eye outline/mask
            bw = return_eye_outline(bf(:,:,1),rect);
            
            r = regionprops(bw,'Centroid');
            
            CM_eye = [r.Centroid];
            
            if show == 1
            % visually inspect that eye outline was determined correctly
            B = bwboundaries(bw); 
            figure(1)
            imshow(bf.*4)
            hold on
            for k=1:length(B)
                b=B{k};
                plot(b(:,2),b(:,1),'r-','LineWidth',2)
            end
            plot(CM_eye(1),CM_eye(2),'g*','MarkerSize',15)
            title(['Finding eye outline and eye center of mass landmark -  Fish number: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',15)
            end
            
            % save eye outline and CM coor in rounds(2)
            switch LR
                case 1
                    FISH(ff).times(t).rounds(2).eye1 = bw;
                    FISH(ff).times(t).rounds(2).eye1_coor = CM_eye;

                case 2
                    FISH(ff).times(t).rounds(2).eye2 = bw;
                    FISH(ff).times(t).rounds(2).eye2_coor = CM_eye;

            end
                        
        end
    end
end

toc; display(' ')

%% Find boundary landmark points and save left-rigth transforms in rounds(2)
% transform images and save transformed images in rounds(3)
%
display('Step 9: Finding boundary landmark points... doing left-rigth side transforms... saving transformed images in rounds(3)...' );

tic
close all

% get width and heigth of images
[ww,hh] = size(FISH(1).times(1).rounds(2).body1);

for ff=1:length(FISH)
    for t=1:3
        
        % retrive images from data struc
        body1  = FISH(ff).times(t).rounds(2).body1;
        body2  = FISH(ff).times(t).rounds(2).body2;
        eye1  = FISH(ff).times(t).rounds(2).eye1;
        eye2  = FISH(ff).times(t).rounds(2).eye2;
        spine1  = FISH(ff).times(t).rounds(2).spine1;
        spine2  = FISH(ff).times(t).rounds(2).spine2;
        bf1  = FISH(ff).times(t).rounds(2).bf1;
        bf2  = FISH(ff).times(t).rounds(2).bf2;
        gfp1  = FISH(ff).times(t).rounds(2).gfp1;
        gfp2  = FISH(ff).times(t).rounds(2).gfp2;
        rfp1  = FISH(ff).times(t).rounds(2).rfp1;
        rfp2  = FISH(ff).times(t).rounds(2).rfp2;
        
        % retrive landmark coor from data struc
        eye_coor1 = FISH(ff).times(t).rounds(2).eye1_coor;
        eye_coor2 = FISH(ff).times(t).rounds(2).eye2_coor;
        spine_end_coor1 = FISH(ff).times(t).rounds(2).spineEnd1_coor;
        spine_end_coor2 = FISH(ff).times(t).rounds(2).spineEnd2_coor;
        
        % make average of eye and spine landmark point the base points
        eye_coor_base = round((eye_coor1+eye_coor2)./2);
        spine_end_coor_base = round((spine_end_coor1+spine_end_coor2)./2);
        
        % function 'return_front_coor' returns boundary input and base
        % points for both right and left sides images
        [x_coor,y_coor1,y_coor2,y_coor_base] = return_front_coor(body1,body2);
        
        % concatenate x- and y- coor of input and base points
        x_coor = [x_coor eye_coor_base(1) spine_end_coor_base(1)];
        y_coor1 = [y_coor1 eye_coor1(2) spine_end_coor1(2)];
        y_coor2 = [y_coor2 eye_coor2(2) spine_end_coor2(2)];
        y_coor_base = [y_coor_base eye_coor_base(2) spine_end_coor_base(2)];
        
        % initialize input point and base point vectors
        input_points1 = zeros(length(x_coor),2);
        input_points2 = zeros(length(x_coor),2);
        base_points   = zeros(length(x_coor),2);
        
        % fill in values
        input_points1(:,1) = x_coor; input_points2(:,1) = x_coor; base_points(:,1) = x_coor;
        input_points1(:,2) = y_coor1; input_points2(:,2) = y_coor2; base_points(:,2) = y_coor_base;
        
        % calculate image tranformation functions
        TFORM1 = cp2tform(input_points1, base_points,'polynomial',2);
        TFORM2 = cp2tform(input_points2, base_points,'polynomial',2);
        
        % apply tranforms to all images
        bodyT1 = imtransform(body1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        bodyT2 = imtransform(body2,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        spineT1 = imtransform(spine1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        spineT2 = imtransform(spine2,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        bfT1 = imtransform(bf1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        bfT2 = imtransform(bf2,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        eyeT1 = imtransform(eye1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfpT1 = imtransform(gfp1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfpT1 = imtransform(rfp1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        eyeT2 = imtransform(eye2,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfpT2 = imtransform(gfp2,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfpT2 = imtransform(rfp2,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        
        % save left-right side transformed images in rounds(3)
        FISH(ff).times(t).rounds(3).bf1   = uint16(bfT1);
        FISH(ff).times(t).rounds(3).bf2   = uint16(bfT2);
        FISH(ff).times(t).rounds(3).gfp1  = uint16(gfpT1);
        FISH(ff).times(t).rounds(3).gfp2  = uint16(gfpT2);
        FISH(ff).times(t).rounds(3).rfp1  = uint16(rfpT1);
        FISH(ff).times(t).rounds(3).rfp2  = uint16(rfpT2);
        FISH(ff).times(t).rounds(3).eye1  = logical(eyeT1);
        FISH(ff).times(t).rounds(3).eye2  = logical(eyeT2);
        FISH(ff).times(t).rounds(3).body1 = logical(bodyT1);
        FISH(ff).times(t).rounds(3).body2 = logical(bodyT2);
        FISH(ff).times(t).rounds(3).spine1 = logical(spineT1);
        FISH(ff).times(t).rounds(3).spine2 = logical(spineT2);
        
        % save transform functions        
        FISH(ff).times(t).transforms.left_rightTF1 = TFORM1;
        FISH(ff).times(t).transforms.left_rightTF2 = TFORM2;
        
        % visually inspect if transform
        % make left and right side fish bodies overlay nicely 
        if show==1
            figure(1)
            subplot(2,1,1); imshowpair(bf1.*2,bf2.*2);hold on;title(ff)
            plot(x_coor,y_coor1,'+r','MarkerSize',15,'LineWidth',1)
            plot(x_coor,y_coor2,'xb','MarkerSize',15,'LineWidth',1)
            plot(x_coor,y_coor_base,'og','MarkerSize',15,'LineWidth',1)
            legend('Rigth side input points','Left side input points','Base points')
            title(['Original images overlay -  Fish number: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',13)
            subplot(2,1,2); imshowpair(bfT1.*2,bfT2.*2);hold on;title(ff)
            plot(x_coor,y_coor1,'+r','MarkerSize',15,'LineWidth',1)
            plot(x_coor,y_coor2,'xb','MarkerSize',15,'LineWidth',1)
            plot(x_coor,y_coor_base,'og','MarkerSize',15,'LineWidth',1)
            title(['Transformed images overlay -  Fish number: ' num2str(ff) ',  Time: ' num2str(t)],'FontSize',13)
        end
        
    end
end

toc; display(' ')
%% Now that left and rigth side images are aligned
% Merge bf, spine, body and eye images and save in rounds(3)
display('Step 10: Making merged versions of left/right side images and saving in rounds(3)...' );

tic
close all

for ff=1:length(FISH)
   for t=1:3
        % retrive images from data struc
        body1  = FISH(ff).times(t).rounds(3).body1;
        body2  = FISH(ff).times(t).rounds(3).body2;
        eye1  = FISH(ff).times(t).rounds(3).eye1;
        eye2  = FISH(ff).times(t).rounds(3).eye2;
        spine1  = FISH(ff).times(t).rounds(3).spine1;
        spine2  = FISH(ff).times(t).rounds(3).spine2;
        bf1  = FISH(ff).times(t).rounds(3).bf1;
        bf2  = FISH(ff).times(t).rounds(3).bf2;
    
        % Merge bf1 and bf2 (pick darkest pixel from each image)
        bf = MergeBF(bf1,bf2);

        % Merge body1, body2 and spine1, spine2  and eye1, eye2 (pick all white pixels from each image)
        body = body1;
        body(body2==1)=1;        
        spine = spine1;
        spine(spine2==1)=1;        
        eye = eye1;
        eye(eye2==1)=1;        
        
        % save in rounds(3)
        FISH(ff).times(t).rounds(3).bf   = uint16(bf);
        FISH(ff).times(t).rounds(3).eye  = logical(eye);
        FISH(ff).times(t).rounds(3).body = logical(body);
        FISH(ff).times(t).rounds(3).spine = logical(spine);
        
   end     
end
toc; display(' ')

%% First scale and rotatate fish so that fish at time point 1, 2 and 3 are the same length
% save scaling transform
% use on body, spine and eye and determine landmarks for nonliniar transform
% save non liniar transform
display('Step 11: Finding transform functions for making all three time point images overlay ...' );

tic
close all

nose_coor_fix = [10,170]; % pick coordinates very close to avearge 'nose' location

for ff=1:length(FISH)
    
    % retrive images from data struc 
    % t=1
    body_t1  = FISH(ff).times(1).rounds(3).body;
    eye_t1   = FISH(ff).times(1).rounds(3).eye;
    spine_t1 = FISH(ff).times(1).rounds(3).spine;
    bf_t1    = FISH(ff).times(1).rounds(3).bf;
    
    % t=2
    body_t2  = FISH(ff).times(2).rounds(3).body;
    eye_t2   = FISH(ff).times(2).rounds(3).eye;
    spine_t2 = FISH(ff).times(2).rounds(3).spine;
    bf_t2    = FISH(ff).times(2).rounds(3).bf;
    
    % t=3
    body_t3  = FISH(ff).times(3).rounds(3).body;
    eye_t3   = FISH(ff).times(3).rounds(3).eye;
    spine_t3 = FISH(ff).times(3).rounds(3).spine;
    bf_t3    = FISH(ff).times(3).rounds(3).bf;
    
    % find end on spine for all three time points
    [y1,x1]=find(spine_t1,10,'last');
    spine_end_coor_t1 = [round(mean(x1)),round(mean(y1))];
    [y2,x2]=find(spine_t2,10,'last');
    spine_end_coor_t2 = [round(mean(x2)),round(mean(y2))];
    [y3,x3]=find(spine_t3,10,'last');
    spine_end_coor_t3 = [round(mean(x3)),round(mean(y3))];
    
    % input points for length scaling transform (make fish same length on t1, t2 and t3 image)
    % nose coor (which is already fixed) and end of spine point
    input_points_t1 = [nose_coor_fix ; spine_end_coor_t1];
    input_points_t2 = [nose_coor_fix ; spine_end_coor_t2];
    input_points_t3 = [nose_coor_fix ; spine_end_coor_t3];
    
    % base points for length scaling transform is nose coor and average
    % fish length. Y-coor of spine end base point is set to
    % nose_coor_fix(2)-20 for all fish
    base_points     = [nose_coor_fix ; [mean([spine_end_coor_t1(1),spine_end_coor_t2(1),spine_end_coor_t3(1)]) nose_coor_fix(2)-20]];
    
    % calculate image tranformation (since we use 'nonreflective similarity' type and there are only two points this transform will be only scaling and rotation)
    TFORM_t1 = cp2tform(input_points_t1, base_points,'nonreflective similarity');
    TFORM_t2 = cp2tform(input_points_t2, base_points,'nonreflective similarity');
    TFORM_t3 = cp2tform(input_points_t3, base_points,'nonreflective similarity');
    
    % save transform functions in data struct.        
    FISH(ff).times(1).transforms.lengthScaleTF = TFORM_t1;    
    FISH(ff).times(2).transforms.lengthScaleTF = TFORM_t2;
    FISH(ff).times(3).transforms.lengthScaleTF = TFORM_t3;
        
    % transform bw images (scaling and rotation)
    body_t1_T   = logical(imtransform(body_t1,TFORM_t1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
    body_t2_T   = logical(imtransform(body_t2,TFORM_t2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
    body_t3_T   = logical(imtransform(body_t3,TFORM_t3,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
    spine_t1_T  = logical(imtransform(spine_t1,TFORM_t1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
    spine_t2_T  = logical(imtransform(spine_t2,TFORM_t2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
    spine_t3_T  = logical(imtransform(spine_t3,TFORM_t3,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
    eye_t1_T    = logical(imtransform(eye_t1,TFORM_t1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
    eye_t2_T    = logical(imtransform(eye_t2,TFORM_t2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
    eye_t3_T    = logical(imtransform(eye_t3,TFORM_t3,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
    
    % return input and base points from fish body boundaries
    [x_coor,y_coor1,y_coor2,y_coor3,y_coor_base] = return_front_coor_T(body_t1_T,body_t2_T,body_t3_T);
    % return input and base points along spine in tail
    [sx_coor,sy_coor1,sy_coor2,sy_coor3,sy_coor_base] = return_spine_coor_T(spine_t1_T,spine_t2_T,spine_t3_T);
    
    r_t1 = regionprops(eye_t1_T);
    eye_coor_t1 = [r_t1.Centroid];
    r_t2 = regionprops(eye_t2_T);
    eye_coor_t2 = [r_t2.Centroid];
    r_t3 = regionprops(eye_t2_T);
    eye_coor_t3 = [r_t3.Centroid];
    
    input_points_t1_T = [[x_coor',y_coor1'];[sx_coor',sy_coor1'];eye_coor_t1];
    input_points_t2_T = [[x_coor',y_coor2'];[sx_coor',sy_coor2'];eye_coor_t2];
    input_points_t3_T = [[x_coor',y_coor3'];[sx_coor',sy_coor3'];eye_coor_t3];
    base_points_T     = [[x_coor',y_coor_base'];[sx_coor',sy_coor_base'];mean([eye_coor_t1;eye_coor_t2;eye_coor_t3])];

    % calculate image tranformation (nonliniar)
    TFORM_t1_T = cp2tform(input_points_t1_T, base_points_T,'polynomial',2);
    TFORM_t2_T = cp2tform(input_points_t2_T, base_points_T,'polynomial',2);
    TFORM_t3_T = cp2tform(input_points_t3_T, base_points_T,'polynomial',2);
    
    % save transform functions in data struct.        
    FISH(ff).times(1).transforms.time123TF = TFORM_t1_T;    
    FISH(ff).times(2).transforms.time123TF = TFORM_t2_T;
    FISH(ff).times(3).transforms.time123TF = TFORM_t3_T;
    
    if show==1       
        % outlines for plot 
        BB1 =  bwboundaries(body_t1);
        BB2 =  bwboundaries(body_t2);
        BB3 =  bwboundaries(body_t3);
        B1 =  bwboundaries(body_t1_T);
        B2 =  bwboundaries(body_t2_T);
        B3 =  bwboundaries(body_t3_T);
        
        figure(1)
        subplot(2,1,1)
        imshow(body_t1+body_t2+body_t3,[])
        hold on
        title(['Overlay of all three timepoints - BEFORE length scaling transform. Fish number: ' num2str(ff) ],'FontSize',10)
        for k=1:length(BB1)
            b=BB1{k};
            plot(b(:,2),b(:,1),'-c')
        end
        for k=1:length(BB2)
            b=BB2{k};
            plot(b(:,2),b(:,1),'-b')
        end
        for k=1:length(BB3)
            b=BB3{k};
            plot(b(:,2),b(:,1),'-k')
        end
        plot(input_points_t1(:,1),input_points_t1(:,2),'+r','MarkerSize',15)
        plot(base_points(:,1),base_points(:,2),'og','MarkerSize',15)
        plot(input_points_t2(:,1),input_points_t2(:,2),'+r','MarkerSize',15)
        plot(input_points_t3(:,1),input_points_t3(:,2),'+r','MarkerSize',15)
        legend('t = 1, outline','t = 2, outline','t = 3, outline','Input points','Base point','Location','NorthEastOutside')
        
        subplot(2,1,2)
        imshow(body_t1_T+body_t2_T+body_t3_T,[])
        hold on
        title(['Overlay of all three timepoints - AFTER length scaling transform. Fish number: ' num2str(ff) ],'FontSize',13)
        for k=1:length(B1)
            b=B1{k};
            plot(b(:,2),b(:,1),'-c')
        end
        for k=1:length(B2)
            b=B2{k};
            plot(b(:,2),b(:,1),'-b')
        end
        for k=1:length(B3)
            b=B3{k};
            plot(b(:,2),b(:,1),'-k')
        end
        plot(input_points_t1(:,1),input_points_t1(:,2),'+r','MarkerSize',15)
        plot(input_points_t2(:,1),input_points_t2(:,2),'+r','MarkerSize',15)
        plot(input_points_t3(:,1),input_points_t3(:,2),'+r','MarkerSize',15)
        plot(base_points(:,1),base_points(:,2),'og','MarkerSize',15)
        
        plot(input_points_t1_T(:,1),input_points_t1_T(:,2),'+r','MarkerSize',15)
        plot(input_points_t2_T(:,1),input_points_t2_T(:,2),'+r','MarkerSize',15)
        plot(input_points_t3_T(:,1),input_points_t3_T(:,2),'+r','MarkerSize',15)
        plot(base_points_T(:,1),base_points_T(:,2),'og','MarkerSize',15)

    end
    
end
toc; display(' ')

%% Use transforms found in previous step. Save transformed images in rounds(4)
display('Step 12: Using transform functions found in previous step for making all three time point images overlay ... Saving transformed images in rounds(4)' );

tic
close all

for ff=1:length(FISH)
    for t=1:3
        
        % retrive images from data struc
        body  = FISH(ff).times(t).rounds(3).body;
        eye  = FISH(ff).times(t).rounds(3).eye;
        spine  = FISH(ff).times(t).rounds(3).spine;
        bf  = FISH(ff).times(t).rounds(3).bf;
        gfp1  = FISH(ff).times(t).rounds(3).gfp1;
        gfp2  = FISH(ff).times(t).rounds(3).gfp2;
        rfp1  = FISH(ff).times(t).rounds(3).rfp1;
        rfp2  = FISH(ff).times(t).rounds(3).rfp2;
                
        % retriave transform functions from data struct.
        lengthScaleTF = FISH(ff).times(t).transforms.lengthScaleTF;
        time123TF     = FISH(ff).times(t).transforms.time123TF;
        
        % transform images - first lengthscale then nonliniar transformation
        body_T   = logical(imtransform(body,lengthScaleTF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        spine_T  = logical(imtransform(spine,lengthScaleTF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        eye_T    = logical(imtransform(eye,lengthScaleTF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        bf_T     = imtransform(bf,lengthScaleTF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfp1_T   = imtransform(gfp1,lengthScaleTF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfp2_T   = imtransform(gfp2,lengthScaleTF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfp1_T   = imtransform(rfp1,lengthScaleTF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfp2_T   = imtransform(rfp2,lengthScaleTF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        
        body_TT   = logical(imtransform(body_T,time123TF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        spine_TT  = logical(imtransform(spine_T,time123TF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        eye_TT    = logical(imtransform(eye_T,time123TF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        bf_TT     = imtransform(bf_T,time123TF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfp1_TT   = imtransform(gfp1_T,time123TF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfp2_TT   = imtransform(gfp2_T,time123TF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfp1_TT   = imtransform(rfp1_T,time123TF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfp2_TT   = imtransform(rfp2_T,time123TF,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        
        % save images in data struc in rounds(4)
        FISH(ff).times(t).rounds(4).body  = body_TT;
        FISH(ff).times(t).rounds(4).eye   = eye_TT;
        FISH(ff).times(t).rounds(4).spine = spine_TT;
        FISH(ff).times(t).rounds(4).bf    = bf_TT;
        FISH(ff).times(t).rounds(4).gfp1  = gfp1_TT;
        FISH(ff).times(t).rounds(4).gfp2  = gfp2_TT;
        FISH(ff).times(t).rounds(4).rfp1  = rfp1_TT;
        FISH(ff).times(t).rounds(4).rfp2  = rfp2_TT;
        
        % visually inspect transformed images
        if show==1
            figure(1)
            subplot(4,1,1)
            hold on
            title('Images AFTER length scaling and','FontSize',13)
            imshow(body_TT+spine_TT+eye_TT,[])
            subplot(4,1,2)
            hold on
            title([' nonliniar transforms. Fish #: ' num2str(ff)  ',  Time: ' num2str(t)],'FontSize',13)
            imshow(bf_TT.*2,[])
            subplot(4,1,3)
            title('RFP')
            hold on
            imshow(rfp1_TT+rfp2_TT,[0 3e4])
            subplot(4,1,4)
            title('GFP')
            hold on
            imshow(gfp1_TT+gfp2_TT,[0 5e4])            
        end
        
    end
end

toc; display(' ')


%% Subtract autoflourescence using RFP image and segment gfp expressing regions in fish 
display('Step 13: Subtracting autoflourescence, segmenting images and saving in rounds(4)...' );

tic
close all

for ff=1:length(FISH)
    for t=1:3
        
        % retrive images from data struc
        body  = FISH(ff).times(t).rounds(4).body;
        eye  = FISH(ff).times(t).rounds(4).eye;
        gfp1  = FISH(ff).times(t).rounds(4).gfp1;
        gfp2  = FISH(ff).times(t).rounds(4).gfp2;
        rfp1  = FISH(ff).times(t).rounds(4).rfp1;
        rfp2  = FISH(ff).times(t).rounds(4).rfp2;
        
        %subtract autoflourescence signal
        im1 = gfp1-rfp1;
        im2 = gfp2-rfp2;
        
        % find boundary points around eye for plotting
        E = bwboundaries(eye);

        % erode eye as to not remove metastasis around eye
        eye = bwmorph(eye,'erode',10);
        
        % function 'return_side()' uses 'return_gradmag_gray()',
        % 'return_varying_thresh()'
        bw1 = return_side(im1,body,eye);
        bw2 = return_side(im2,body,eye);
        
        % merge left and right segmented images
        bw = bw1;
        bw(bw2==1)=1;
        
        if t==3 % be more critical at last time point very small objects (<15 pixels) could be dirt/fish poop and we do not have another time point to look at to see if they stick around...
            bw = bwmorph(bw,'erode',1);
            bw = bwareaopen(bw,15);
            bw = imfill(bw,'holes');
        else     % remove very very small objects - they are most likely dirt
            bw = bwareaopen(bw,5);
            bw = imfill(bw,'holes'); 
        end
        
        % take brigthtest pixel from the left and right pure signal images
        % and make merged image
        im = Merge(im1,im2);

        % save merged pure signal and merged segmented image in rounds(4) 
        FISH(ff).times(t).rounds(4).bw = logical(bw);
        FISH(ff).times(t).rounds(4).pure = uint16(im);

        % Show images with autoflourescence subtracted
        if show==1                      
            B = bwboundaries(body);            
            figure
            imshow(im.*7)
            hold on
            title(['GFP signal left/rigth merged, minus autoflourescence. Fish #: ' num2str(ff)  ',  Time: ' num2str(t)],'FontSize',13)
            b=B{1};
            plot(b(:,2),b(:,1),'b-')
            b=E{1};
            plot(b(:,2),b(:,1),'r-')
            
        end
    end
end

toc; display(' ')

%% Mask day 1 using day 7 bw image and mask day 7 bw image using day 14 bw image
%we assume that object not segmented a week later where noise when they
%where first detected. Save in rounds(4) (NB thus overwriting previous image)

display('Step 14: Removing objects which are no longer there in the next timepoint...' );

tic
close all

for ff=1:length(FISH)
    
    % retrive images from data struc
    bw_t1  = FISH(ff).times(1).rounds(4).bw;
    bw_t2  = FISH(ff).times(2).rounds(4).bw;
    bw_t3  = FISH(ff).times(3).rounds(4).bw;    
    im_t1  = FISH(ff).times(1).rounds(4).pure;
    im_t2  = FISH(ff).times(2).rounds(4).pure;
    im_t3  = FISH(ff).times(3).rounds(4).pure;
    
    % dilate masks to allow for inclusion of regions which where not
    % perfectly mapped on top of each other by the transformations
    mask_t3 = bwmorph(bw_t3,'dilate',1);
    mask_t2 = bwmorph(bw_t2,'dilate',1);
    
    % mask_t2 should not have objects which are not present in mask_t3
    mask_t2(mask_t3==0)=0;
    
    clean_bw_t1 = bw_t1;
    clean_bw_t2 = bw_t2;
    
    % clean up the past
    clean_bw_t1(mask_t2==0)=0;
    clean_bw_t2(mask_t3==0)=0;
    
    % show final segmentation result on top of pure=gfp-rfp images
    if show==1
        C1 = bwboundaries(clean_bw_t1);
        C2 = bwboundaries(clean_bw_t2);
        C3 = bwboundaries(bw_t3);             
        figure
        subplot(3,1,1)
        imshow(im_t1.*7)
        hold on
        title(['Green line: segmentation. Fish #: ' num2str(ff)  ',  Time: ' num2str(1)],'FontSize',13)
        for k=1:length(C1)
            b=C1{k};
            plot(b(:,2),b(:,1),'-g')
        end        
        subplot(3,1,2)
        imshow(im_t2.*7)
        hold on
        title(['Fish #: ' num2str(ff)  ',  Time: ' num2str(2)],'FontSize',13)
        for k=1:length(C2)
            b=C2{k};
            plot(b(:,2),b(:,1),'-g')
        end        
        subplot(3,1,3)
        imshow(im_t3.*7)
        hold on
        title(['Fish #: ' num2str(ff)  ',  Time: ' num2str(3)],'FontSize',13)
        for k=1:length(C3)
            b=C3{k};
            plot(b(:,2),b(:,1),'-g')
        end
    end
    
    % save cleaned up segmentation from t=1 and t=2 images in rounds(4)
    FISH(ff).times(1).rounds(4).bw = logical(clean_bw_t1);
    FISH(ff).times(2).rounds(4).bw = logical(clean_bw_t2);
    
end

toc; display(' ')

%%  Construct DistMatrix form BW segmented images (pairwise distance
% matrix holding the shortest distances between all seperate regions with
% 1's. (distMatrix is a symmetric matrix))
% Distances are shortest possible distances between points on the
% boundary of both objects
% use this matrix to make labeled segmented image where each cluster of
% close objects have same value 1,2,3 ect and the rest is 0
close all
display('Step 15: Finding pairwise distance matrix for all connected regions in segmented images and use this to cluster close objects...' );

tic

% threshold (in pixels) for clustering close objects
threshold = 10;

for ff=1:length(FISH)
    for t=1:3
        
        % retrive images from data struc
        bw    = FISH(ff).times(t).rounds(4).bw;
        body  = FISH(ff).times(t).rounds(4).body;
        
        % function 'return_dist_matrix' finds matrix holding the shortest distances between all seperate regions with
        % 1's. Distances are shortest possible distances between points on the
        % boundary of both objects. If there is only one object in fish (only
        % primary) distmatrix is empty []
        distMatrix = return_dist_matrix(bw);        
        
        L_clustered = return_clustered_L(bw,distMatrix,threshold,t);
        
        % save matrix and L_clustered in data struct
        FISH(ff).times(t).rounds(4).distMatrix = distMatrix;
        FISH(ff).times(t).rounds(4).L = L_clustered;
       
        if show==1
            figure(1)
            subplot(3,1,t)
            imshow(L_clustered,[])
            hold on
            title(['Labeled and clustered segmented image. Fish #: ' num2str(ff)  ',  Time: ' num2str(t)],'FontSize',13)
            colormap('jet')
            
        end
        
    end
end

toc; display(' ')
%%
close all
display('Step 16: Separating objects which has merged due to melanoma growth since last timepoint...' );
tic

for ff=1:length(FISH)
    
    % retrive images from data struc
    L_t1 = FISH(ff).times(1).rounds(4).L;
    L_t2 = FISH(ff).times(2).rounds(4).L;
    L_t3 = FISH(ff).times(3).rounds(4).L;
    
    % function return_L_divided (uses the functions: return_met_events, return_clustered_C, return_object_voronoi_points and return_labeled_voronoi)
    % The function seperates objects in labeled matrix L_t which started out as
    % seperated object in L_t-1 but has merged over time. When two objects have fully merged (no 0's between then) seperation line is found using
    % dilation of ojects in the previous timepoints and using the
    % boundary points of these objects to determine voronoi regions
    % and assign these regions to either of the two objects
    % (see more information in the individual functions listed above)
    [L_divided_t1,L_divided_t2,L_divided_t3,V_merge_t2,V_merge_t3]  = return_L_divided(L_t1,L_t2,L_t3);
    
    if show==1
        figure
        subplot(3,1,1)
        imshow(L_divided_t1,[])
        hold on
        title(['Labeled and clustered segmented image. Fish #: ' num2str(ff)  ',  Time: ' num2str(1)],'FontSize',13)
        colormap('jet')
        subplot(3,1,2)
        imshow(L_divided_t2+V_merge_t2,[])
        hold on
        title(['Labeled, clustered and divided segmented image. Fish #: ' num2str(ff)  ',  Time: ' num2str(2)],'FontSize',13)
        colormap('jet')
        subplot(3,1,3)
        imshow(L_divided_t3+V_merge_t3,[])
        hold on
        title(['Labeled, clustered and divided segmented image. Fish #: ' num2str(ff)  ',  Time: ' num2str(3)],'FontSize',13)
        colormap('jet')
    end
    
    % save labeled segmented image in data struct.
    FISH(ff).times(1).rounds(4).L_divided = int8(L_divided_t1);
    FISH(ff).times(2).rounds(4).L_divided = int8(L_divided_t2);
    FISH(ff).times(3).rounds(4).L_divided = int8(L_divided_t3);
      
end


toc; display(' ')

%% Transform to standard fish shape
% First transform all fish to have same length
% the determine average fish shape and transform all images to this
% shape. If the boundary of a fish bulges out over the standard fish shape
% close to the injection site we do not use standard shape as base coor in
% this region - this way we avoid transforming a bulging tumor smaller 

close all
display('Step 17: Finding average fish shape of entire data set and transforming to this shape...' );


tic

% find average fish length
FL=[];% vector for storing fish lengths
cc=0; % index for fish length vector

for ff=1:length(FISH)
    for t=1:3
        body  = FISH(ff).times(t).rounds(4).body;        
        cc=cc+1;
        FL(cc)= return_fish_length(body);
    end
end

% find average fish length
ave_fish_length = round(mean(FL));
nose_coor_fix = [10,170]; % pick coordinates very close to avearge 'nose' location
end_spine_coor_fix = [ave_fish_length, 150];

% initialize sum images
bodySum  = zeros(size(FISH(1).times(1).rounds(4).body));
eyeSum   = zeros(size(FISH(1).times(1).rounds(4).body));
spineSum = zeros(size(FISH(1).times(1).rounds(4).body));

for ff=1:length(FISH)
    for t=1:3
        % retrive images from data struc
        body   = FISH(ff).times(t).rounds(4).body;
        spine  = FISH(ff).times(t).rounds(4).spine;
        eye    = FISH(ff).times(t).rounds(4).eye;
        bf     = FISH(ff).times(t).rounds(4).bf;
        
        % find coordinates of end of spine
        [y,x]=find(spine,20,'last');
        end_spine_coor = [round(mean(x)),round(mean(y))];
        
        % calculate image tranformation (nonliniar)
        TFORM = cp2tform([nose_coor_fix;end_spine_coor], [nose_coor_fix;end_spine_coor_fix],'nonreflective similarity');
        
        % save liniar transform function in data struct.
        FISH(ff).times(t).transforms.lengthScaleTF2 = TFORM;
        
        % transform all fish bodies to be the same length
        body_T  = logical(imtransform(body,TFORM,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        eye_T   = logical(imtransform(eye,TFORM,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        spine_T = logical(imtransform(spine,TFORM,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        bf_T    = uint16(imtransform(bf,TFORM,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        
        % add all length-rescaled fish bodies, eye and spine in bodySum,
        % eyeSum and spineSum images
        bodySum  = bodySum + body_T - eye_T;
        eyeSum   = eyeSum + eye_T;
        spineSum = spineSum + spine_T;
        
        if show==1
            figure(1)
            imshow(bodySum,[])
            hold on
            plot(end_spine_coor(1),end_spine_coor(2),'b*','MarkerSize',20)
            plot(end_spine_coor_fix(1),end_spine_coor_fix(2),'r*','MarkerSize',20)
            legend('Input points for lengthScaleTF2','Base points for lengthScaleTF2')
            plot(nose_coor_fix(1),nose_coor_fix(2),'b*','MarkerSize',20)
            plot(nose_coor_fix(1),nose_coor_fix(2),'r*','MarkerSize',20)
            title('Finding average fish shape and average eye shape and eye location','FontSize',15);
        end
    end
end

% Using the image bodySum find region where 50% or more of the rescaled fish bodies in the dataset overlap
% this will be the standard fish shape which we map all bodies to
bodySum = bodySum./max(bodySum(:));
eyeSum = eyeSum./max(eyeSum(:));
spineSum = spineSum./max(spineSum(:));
aveFish = smoothBW(im2bw(bodySum,0.5),30);
aveEye = smoothBW(im2bw(eyeSum,0.5),10);
aveSpine = smoothBW(im2bw(spineSum,0.9),10);

% find center of mass of average eye
r = regionprops(aveEye,'Centroid');
aveEyeCoor = [r.Centroid];

if show==1 % plot stadard fish shape
    A = bwboundaries(aveFish);
    E = bwboundaries(aveEye);
    S = bwboundaries(aveSpine);
    hold on
    b=A{1};
    plot(b(:,2),b(:,1),'-y','LineWidth',2)
    b=E{1};
    plot(b(:,2),b(:,1),'-y','LineWidth',2)
    plot(aveEyeCoor(1),aveEyeCoor(2),'oy','LineWidth',2)
    hold on
    b=S{1};
    plot(b(:,2),b(:,1),'-y','LineWidth',2)
end

for ff=1:length(FISH)
    for t=1:3
        % retrive images from data struct.
        body   = FISH(ff).times(t).rounds(4).body;
        spine   = FISH(ff).times(t).rounds(4).spine;
        eye   = FISH(ff).times(t).rounds(4).eye;
        bf   = FISH(ff).times(t).rounds(4).bf;
        
        % retrive liniar transform function in data struct.
        TFORM1 = FISH(ff).times(t).transforms.lengthScaleTF2;
        
        % transform all fish bodies to be the same length
        body_T = logical(imtransform(body,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        spine_T = logical(imtransform(spine,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        eye_T = logical(imtransform(eye,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        bf_T = imtransform(bf,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        
        % find center of mass of eye
        r = regionprops(eye_T,'Centroid');
        eye_T_Coor = [r.Centroid];
        
        % find input points for transform from fish body and spine in tail
        [input_points,base_points] = return_front_coor_Std(body_T,aveFish);
        [x_coor,y_coor,y_coor_base] = return_spine_coor_Std(spine_T,aveSpine); % uses function return_Y_coor()
        
        % concatenate points
        input_points = [input_points ;[x_coor; y_coor]'     ;eye_T_Coor];
        base_points  = [base_points  ;[x_coor; y_coor_base]';aveEyeCoor];
        
        % calculate image tranformation (nonliniar)
        TFORM2 = cp2tform(input_points, base_points,'polynomial',2);
        
        % save non liniar transform function in data struct.
        FISH(ff).times(t).transforms.standardFishTF = TFORM2;

        if show==1
            % test non liniar transform on bf image and plot
            bf_TT = imtransform(bf_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);

            figure(1)
            imshowpair(bf_T.*2,bf_TT.*2)
            hold on
            b=A{1};
            plot(b(:,2),b(:,1),'-y','LineWidth',2)
            plot(input_points(:,1),input_points(:,2),'b*','MarkerSize',20)
            plot(base_points(:,1),base_points(:,2),'r*','MarkerSize',20)
            legend('Input points','Base points')
            title('Transforming to standard fish shape. (Avoid transforming fish shape close to primary tumor)','FontSize',15);
        end
    end
end

toc; display(' ')
            
%% apply standard fish transform to segmented images

close all
display('Step 18: Applying transform functions found in previos step to all images and saving in finals...' );
tic

% boundary points of standard fish shape
A = bwboundaries(bwmorph(aveFish,'dilate',2));

for ff=1:length(FISH)
    for t=1:3
        
        % retriavae images from data struct.
        body   = FISH(ff).times(t).rounds(4).body;
        spine   = FISH(ff).times(t).rounds(4).spine;
        eye   = FISH(ff).times(t).rounds(4).eye;
        bw   = FISH(ff).times(t).rounds(4).bw;
        bf   = FISH(ff).times(t).rounds(4).bf;
        gfp1   = FISH(ff).times(t).rounds(4).gfp1;
        gfp2   = FISH(ff).times(t).rounds(4).gfp2;
        rfp1   = FISH(ff).times(t).rounds(4).rfp1;
        rfp2   = FISH(ff).times(t).rounds(4).rfp2;
        pure   = FISH(ff).times(t).rounds(4).pure;
        L      = FISH(ff).times(t).rounds(4).L;
        L_divided      = FISH(ff).times(t).rounds(4).L_divided;
        
        % retriavae transform functions from data struct.
        TFORM1 = FISH(ff).times(t).transforms.lengthScaleTF2;
        TFORM2 = FISH(ff).times(t).transforms.standardFishTF;
        
        % transform images
        bf_T    = imtransform(bf,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        bf_TT   = imtransform(bf_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfp1_T  = imtransform(gfp1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfp1_TT = imtransform(gfp1_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfp2_T  = imtransform(gfp2,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        gfp2_TT = imtransform(gfp2_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfp1_T  = imtransform(rfp1,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfp1_TT = imtransform(rfp1_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfp2_T  = imtransform(rfp2,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        rfp2_TT = imtransform(rfp2_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        pure_T  = imtransform(pure,TFORM1,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        pure_TT = imtransform(pure_T,TFORM2,'XData',[1 hh], 'YData',[1 ww],'FillValues',0);
        body_T  = logical(imtransform(body,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        body_TT = logical(imtransform(body_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        spine_T  = logical(imtransform(spine,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        spine_TT = logical(imtransform(spine_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        eye_T  = logical(imtransform(eye,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        eye_TT = logical(imtransform(eye_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        bw_T  = logical(imtransform(bw,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        bw_TT = logical(imtransform(bw_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        L_T  = uint8(imtransform(L,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        L_TT = uint8(imtransform(L_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        L_divided_T  = uint8(imtransform(L_divided,TFORM1,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        L_divided_TT = uint8(imtransform(L_divided_T,TFORM2,'nearest','XData',[1 hh], 'YData',[1 ww],'FillValues',0));
        
        % save transformed images in finals
        FISH(ff).times(t).finals.body  = body_TT;
        FISH(ff).times(t).finals.spine = spine_TT;
        FISH(ff).times(t).finals.eye   = eye_TT;
        FISH(ff).times(t).finals.bw    = bw_TT;
        FISH(ff).times(t).finals.bf    = bf_TT;
        FISH(ff).times(t).finals.gfp1  = gfp1_TT;
        FISH(ff).times(t).finals.gfp2  = gfp2_TT;
        FISH(ff).times(t).finals.rfp1  = rfp1_TT;
        FISH(ff).times(t).finals.rfp2  = rfp2_TT;
        FISH(ff).times(t).finals.pure  = pure_TT;
        FISH(ff).times(t).finals.L          = L_TT;
        FISH(ff).times(t).finals.L_divided  = L_divided_TT;
        
        % show L_divided_TT in standard fish shape
        if show==1
            figure(1)
            imshow(L_divided_TT,[])
            hold on
            colormap('jet')
            title(['Labeled, clustered and divided segmented image transformed to standard fish shape. Fish #: ' num2str(ff)  ',  Time: ' num2str(t)],'FontSize',13);
            b=A{1};
            plot(b(:,2),b(:,1),'-y','LineWidth',2)
        end
    end
end
toc; display(' ')

%%
display([ 'Step 19: Saving data struct. FISH in current directory as ' file_name_for_final_data_struct]);
tic

% save data struct FISH as a .mat file
save(file_name_for_final_data_struct,'FISH','-v7.3')
toc


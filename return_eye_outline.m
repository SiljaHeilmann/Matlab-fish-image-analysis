function bw = return_eye_outline(bf_red,rect)
% takes rect and red channel of rgb bf image and returns eye outline
            
            crop = imcrop(bf_red,rect);       
                   
            bw = im2bw(bf_red,graythresh(crop));
       
            bw = bw.*-1 +1;
            
            bw(:,1:round(rect(1)/2)) = 0;
            bw(:,rect(3)+round(rect(4)/2):end) = 0;
            bw(1:rect(3)+rect(4),:) = 0;
            bw(end-rect(3)+rect(4):end,:) = 0;
            
            bw = bwmorph(bw,'dilate',1);
            
            bw = imfill(bw,'holes');
         
            bw = bwareaopen(bw,300);

            bw = bwmorph(bw,'erode',2);
            bw = bwareaopen(bw,300);
            
            bw = bwmorph(bw,'dilate',2);
            
            bw = smoothBW(bw,10);



end


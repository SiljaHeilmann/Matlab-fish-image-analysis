function ang = find_rot_angle(bodyTail_bw) 
% input is rough bw mask of fish body (with or without tail fin - doesnt matter)
% output is ang - how much images needs to be rotated so that fish is approx aligned with horizontal 

% erode body (find 'core' of body and tail - get rid of tail fin)
bodyTail_bw = bwmorph(bodyTail_bw,'erode',20); % not strictly nessesary
bodyTail_bw = bwareaopen(bodyTail_bw,20000);   % not strictly nessesary
 
% find nose and tail coor (for finding length of fish and later angle)
[y,x]=find(bodyTail_bw,10,'first');
nose_coor = [round(mean(x)),round(mean(y))];
[y,x]=find(bodyTail_bw,1000,'last');
tail_coor = [round(mean(x)),round(mean(y))];

fish_length = tail_coor(1)-nose_coor(1);

% point in the middle of the fish body at 0.35 of the full fish
% length from the 'nose'
midbody_coor(1) = nose_coor(1) + round(0.35*fish_length);

First = find(bodyTail_bw(:,midbody_coor(1)),1,'first');
Last  = find(bodyTail_bw(:,midbody_coor(1)),1,'last');

midbody_coor(2) = First+round(abs(First-Last)/2);


% find angle of fish with respect to horizontal using midbody coor
a = midbody_coor(1)-nose_coor(1);
b = midbody_coor(2)-nose_coor(2);
c = sqrt(a^2+b^2);
ang1 = 0;
if b>=0
    ang1 = acos(a/c)*(360/(2*pi));
elseif b<0
    ang1 = -acos(a/c)*(360/(2*pi));
end

% find angle of fish with respect to horizontal using tail coor
a = tail_coor(1)-nose_coor(1);
b = tail_coor(2)-nose_coor(2);
c = sqrt(a^2+b^2);
ang2 = 0;
if b>=0
    ang2 = acos(a/c)*(360/(2*pi));
elseif b<0
    ang2 = -acos(a/c)*(360/(2*pi));
end

% use mean of angle found using tail and midbody coor
ang = mean([ang1 ang2]);


% plot midbody and tail coor
% figure(1)
% imshow(bodyTail_bw)
% hold on
% plot(midbody_coor(1),midbody_coor(2),'r*')
% plot(tail_coor(1),tail_coor(2),'g*')


end
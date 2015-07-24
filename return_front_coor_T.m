function [x_coor,y_coor1,y_coor2,y_coor3,y_coor_base] = return_front_coor_T(body_bw1,body_bw2,body_bw3)

fish_length1 = return_fish_length(body_bw1);
fish_length2 = return_fish_length(body_bw2);
fish_length3 = return_fish_length(body_bw3);

fish_length = min([fish_length1 fish_length2 fish_length3]);

start = 40;
slut = round(fish_length*0.99); %

x_coor = round(linspace(start+1,slut-1,11)); %5

% dont put points near tumor
x_coor(7:8)=[];


yup1_coor = zeros(1,length(x_coor));

ydown1_coor = zeros(1,length(x_coor));

yup2_coor = zeros(1,length(x_coor));

ydown2_coor = zeros(1,length(x_coor));

yup3_coor = zeros(1,length(x_coor));

ydown3_coor = zeros(1,length(x_coor));

for i=1:length(x_coor)
        yup1_coor(i) = find(body_bw1(:,x_coor(i)),1,'first');
        yup2_coor(i) = find(body_bw2(:,x_coor(i)),1,'first');
        yup3_coor(i) = find(body_bw3(:,x_coor(i)),1,'first');
        ydown1_coor(i) = find(body_bw1(:,x_coor(i)),1,'last');
        ydown2_coor(i) = find(body_bw2(:,x_coor(i)),1,'last');    
        ydown3_coor(i) = find(body_bw3(:,x_coor(i)),1,'last');    
end

x_coor = [x_coor x_coor];

y_coor1 = [yup1_coor ydown1_coor];

y_coor2 = [yup2_coor ydown2_coor];

y_coor3 = [yup3_coor ydown3_coor];


y_coor_base = round((y_coor1 + y_coor2+ y_coor3)./3);
 
end


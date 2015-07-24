function [x_coor,y_coor1,y_coor2,y_coor3,y_coor_base] = return_spine_coor_T(spine_bw1,spine_bw2,spine_bw3)

     spine1 = sum(spine_bw1(:,200:end));
     spine2 = sum(spine_bw2(:,200:end));
     spine3 = sum(spine_bw3(:,200:end));

     index1=find(spine1);
     spine1 = spine1(index1);
     index2=find(spine2);
     spine2 = spine2(index2);
     index3=find(spine3);
     spine3 = spine3(index3);
              
         
     width1 = spine1(round(0.82*length(spine1)))*1.5;
     width2 = spine2(round(0.82*length(spine2)))*1.5;
     width3 = spine3(round(0.82*length(spine3)))*1.5;
   
     start1 = find(spine1<width1,1) +200;
     start2 = find(spine2<width2,1) +200;
     start3 = find(spine3<width3,1) +200;

     slut1 = find(spine1,1,'last')+200-20;
     slut2 = find(spine2,1,'last')+200-20;
     slut3 = find(spine3,1,'last')+200-20;
     
     start = max([start1 start2 start3]) ;
     slut = min([slut1 slut2 slut3]);
     
     x_coor = round(linspace(start+1,slut-1,5));
     
     y_coor1 = x_coor.*0;
     y_coor2 = y_coor1;
     y_coor3 = y_coor1;


         for n=1:length(x_coor)
             y_coor1(n) = return_Y_coor(x_coor(n),spine_bw1);             
             y_coor2(n) = return_Y_coor(x_coor(n),spine_bw2);             
             y_coor3(n) = return_Y_coor(x_coor(n),spine_bw3);             
         end
         
         
     y_coor_base = round((y_coor1 + y_coor2 + y_coor3)./3);
         

end


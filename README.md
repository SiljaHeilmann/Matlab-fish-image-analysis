# Matlab-fish-image-analysis

Run only this file: CreateFISHdataStruc.m - it calls all other funcIons.

Make sure that current matlab contains folder with image.

Important: the script CreateFISHdataStruc.m require the following matlab toolboxes to run:
• Matlab statistics toolbox
• Matlab imaging processing toolbox


Functions called by 'CreateFISHdataStruc.m':
find_rot_angle.m
Merge.m
MergeBF.m
return_side.m
return_body_outlines.m
return_clustered_L.m
return_dist_matrix.m
return_eye_outline.m
return_fish_length.m
return_front_coor_Std.m
return_front_coor_T.m
return_front_coor.m
return_L_divided.m
return_spine_coor_Std.m
return_spine_coor_T.m
return_spine_outline.m
smoothBW.m
      
      Functions called only inside other functions:
      return_clustered_C.m
      return_gradmag_gray.m
      return_labeled_voronoi.m
      return_met_events.m
      return_object_voronoi_points.m
      return_varying_thresh.m
      return_Y_coor.m
      return_sub_listL.m
      smooth_boundary.m

function threshold=thres_calculator(base_threshold,dist_threshold,cur_dis)

threshold=((1-base_threshold)/dist_threshold*(45-cur_dis))+base_threshold;


end
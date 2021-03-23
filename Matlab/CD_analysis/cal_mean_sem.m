function [lap_range CA1_mean_fy CA1_f_sem]=cal_mean_sem(CA1_f_x, CA1_f_y)
 lap_range=1:1:max(CA1_f_x);
[count id]=histc(CA1_f_x,lap_range);
CA1_mean_fy=[];
CA1_f_sem=[];
for i=lap_range
    cur_fy=CA1_f_y(id==i);
    CA1_mean_fy=[CA1_mean_fy mean(cur_fy)];
    CA1_f_sem= [CA1_f_sem std(cur_fy)/sqrt(length(cur_fy))];
end

end
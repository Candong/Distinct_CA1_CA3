function [x y]=find_center( ROI )
rol_factor=(1:size(ROI,1));
col_factor=(1:size(ROI,2));
mean_ROI=mean(ROI(:));
% norm_ROI=ROI/max(max(ROI));
% % [rol col]=find(norm_ROI>0);
% pos=(norm_ROI>0).*norm_ROI;
% 
% w_x=rol_factor.*norm_ROI;
% w_y=col_factor.*norm_ROI;
% 
% x_pos=w_x(find(w_x>0));
% y_pos=w_y(find(w_y>0));
[X Y]=meshgrid(col_factor,rol_factor);
x=mean(mean(ROI.*X))/mean_ROI;
y=mean(mean(ROI.*Y))/mean_ROI;

% y=mean(rol);
% x=mean(col);


end
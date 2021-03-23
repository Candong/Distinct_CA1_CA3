clear all;close all;

[PC_filepaths, temp]=(uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on'));

if ~isa(PC_filepaths,'cell') 
  PC_filepaths={PC_filepaths};  
end 
window=5;
binM=1:50;
total_lap_width=[];
total_COM=[];
com=[];
start_bin=[];
end_bin=[];

count=0;
for f=1:size(PC_filepaths,2)
load([temp PC_filepaths{f}]);
for i=1:size(sig_PFs,2)  
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean=sig_PFs{j,i};
            width_lap=find_pfwidth(binmean,window);
            total_lap_width=[total_lap_width; width_lap];
            
            
            COM_lap=find_pfCOM(binmean,binM);
            total_COM=[total_COM COM_lap];

            [meanCOM SP]=meanCOMandSP(binmean',1,50,50);
            com=[com meanCOM];
            start_bin=[start_bin PF_start_bins(j,i)];
            end_bin=[end_bin PF_end_bins(j,i)];
            if PF_start_bins(j,i)==1
                if (PF_end_bins(j,i)-meanCOM)>(meanCOM-PF_start_bins(j,i))
                    count=count+1;
                end
            elseif PF_start_bins(j,i)==50
                if (PF_end_bins(j,i)-meanCOM)<(meanCOM-PF_start_bins(j,i))
                    count=count+1;
                end
                
            end
                
               

        end
    end
end
end

figure;

plot(1:size(total_lap_width,2),(total_lap_width-total_lap_width(:,15))')
% figure;
% 
% plot(1:34,mean((total_lap_width-total_lap_width(:,15))',2))

figure;
total_COM=total_COM';
plot(1:size(total_COM,2),(total_COM-total_COM(:,15))')

% figure;
% total_COM=total_COM';
% plot(1:35,mean(total_COM-total_COM(:,15),1)')
figure;
norm_lap=15;
total_lap_width(total_lap_width==0)=NaN;
for i = 1:size(total_lap_width,1)
    cur_y=total_lap_width(i,:);
    cur_x=1:size(total_lap_width,2);
    remove_id=cur_y==NaN;
    cur_y=cur_y-cur_y(norm_lap);
    cur_y(remove_id)=[];
    cur_x(remove_id)=[];
    
 model=fitlm(cur_x,cur_y);
hold on;
a=model.('Coefficients').('Estimate')(2);
plot(cur_x,model.('Coefficients').('Estimate')(2).*cur_x+ model.('Coefficients').('Estimate')(1));
title(['slope=' num2str(a)]);
% ylim([-50 50]);
end

figure;
norm_lap=1;
total_rsquare=[];
for i = 1:size(total_COM,1)
    cur_y=total_COM(i,:);
    cur_x=1:size(total_COM,2);
    remove_id=isnan(cur_y);
    cur_y=cur_y-cur_y(norm_lap);
    cur_y(remove_id)=[];
    cur_x(remove_id)=[];
    cur_x(isnan(cur_y))=[];
    cur_y(isnan(cur_y))=[];
 if ~isempty(cur_y)
 model=fitlm(cur_x,cur_y,'Intercept',false);
hold on;
%a=model.('Coefficients').('Estimate')(1);
%plot(cur_x,model.('Coefficients').('Estimate')(2).*cur_x+ model.('Coefficients').('Estimate')(1));
plot(cur_x,model.('Coefficients').('Estimate')(1).*cur_x);
%title(['slope=' num2str(a)]);
% ylim([-50 50]);
total_rsquare=[total_rsquare model.Rsquared.Ordinary];
end


end

function COM_lap=find_pfCOM(binmean,binM)


COM_lap=sum(binmean'.*binM,2)./sum(binmean',2);

end
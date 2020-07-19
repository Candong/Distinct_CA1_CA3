clear all; %close all;

[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose same threshold files to load:','MultiSelect','on');
range=-1:50;

plane_name={'plain1','plain2','plain3'};
figure
for p=1:3
f_distribute=zeros(size(behavior_filepaths,2),length(range));
f_distribute_norm=zeros(size(behavior_filepaths,2),length(range));
n_distribute=zeros(size(behavior_filepaths,2),length(range));
n_distribute_norm=zeros(size(behavior_filepaths,2),length(range));
combine_distribute=zeros(size(behavior_filepaths,2),length(range));
combine_distribute_norm=zeros(size(behavior_filepaths,2),length(range));

for f=1:size(behavior_filepaths,2)
 load([temp behavior_filepaths{f}]);  
 for i=1:size(PF_info,2)
     if contains(PF_info(i).filepaths,'_f_')
      
        if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,plane_name{p}) & ~isempty(PF_info(i).delay_lap_distribute)
           f_distribute(f,:)=PF_info(i).delay_lap_distribute;
           f_distribute_norm(f,:)=PF_info(i).delay_lap_distribute/sum(PF_info(i).delay_lap_distribute);
              
          end
        end
        
     elseif contains(PF_info(i).filepaths,'_n_')
      
        if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,plane_name{p}) & ~isempty(PF_info(i).delay_lap_distribute)
           n_distribute(f,:)=PF_info(i).delay_lap_distribute;
           n_distribute_norm(f,:)=PF_info(i).delay_lap_distribute/sum(PF_info(i).delay_lap_distribute);
              
          end
          
        elseif contains(PF_info(i).filepaths,'2nd25lap')
          
          if contains(PF_info(i).filepaths,plane_name{p}) & ~isempty(PF_info(i).delay_lap_distribute)
           combine_distribute(f,:)=PF_info(i).delay_lap_distribute;
           combine_distribute_norm(f,:)=PF_info(i).delay_lap_distribute/sum(PF_info(i).delay_lap_distribute);
              
          end
 
        end
        
     end 
     
     
 end
    
end
F_d{p}=f_distribute;
F_d_norm{p}=f_distribute_norm;
N_d{p}=n_distribute;
N_d_norm{p}=n_distribute_norm;
Comb_d{p}=combine_distribute;
Comb_d_norm{p}=combine_distribute_norm;

subplot(3,3,p)
bar(f_distribute')
title(['familiar plane' num2str(p)]);
ylim([0 25]);

subplot(3,3,3+p)
bar(n_distribute')
title(['novel plane' num2str(p)]);
ylim([0 25]);

subplot(3,3,6+p)
bar(combine_distribute')
title(['combine plane' num2str(p)]);
ylim([0 25]);

end


F=F_d{1,1}+F_d{1,2}+F_d{1,3};
figure
bar(F')
ylim([0 40])
title('F')

figure
bar(sum(F,1)')
ylim([0 150])
title('F')

N=N_d{1,1}+N_d{1,2}+N_d{1,3};
figure
bar(N')
ylim([0 40])
title('N')

figure
bar(sum(N,1)')
ylim([0 150])
title('N')

F_norm=F./sum(F,2);
N_norm=N./sum(N,2);
figure;
bar(mean(F_norm,1));
hold on; plot(F_norm','o')
ylim([0 1])

figure;
bar(mean(N_norm,1));
hold on; plot(N_norm','o')
ylim([0 1])

figure;
cdfplot(mean(F_norm(:,2:end),1));
hold on;
cdfplot(mean(N_norm(:,2:end),1));
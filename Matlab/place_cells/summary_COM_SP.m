clear all; close all;

[filepaths, temp]=uigetfile('*.mat', 'Chose COM files to load:','MultiSelect','on');
   p1_f_SP=[];
   p1_n_SP=[];
   p2_f_SP=[];
   p2_n_SP=[];
   p3_f_SP=[];
   p3_n_SP=[];
for f=1:size(filepaths,2)
load([temp filepaths{f}]);  

   p1_f_SP=[p1_f_SP f_SP{1}];
   p1_n_SP=[p1_n_SP n_SP{1}];
   
   p2_f_SP=[p2_f_SP f_SP{2}];
   p2_n_SP=[p2_n_SP n_SP{2}];
   
   p3_f_SP=[p3_f_SP f_SP{3}];
   p3_n_SP= [p3_n_SP n_SP{3}];

   
end

p1_f_SP(find(isnan(p1_f_SP)))=[];
p1_f_SP(find(isinf(p1_f_SP)))=[];

p1_n_SP(find(isnan(p1_n_SP)))=[];
p1_n_SP(find(isinf(p1_n_SP)))=[];

p2_f_SP(find(isnan(p2_f_SP)))=[];
p2_f_SP(find(isinf(p2_f_SP)))=[];

p2_n_SP(find(isnan(p2_n_SP)))=[];
p2_n_SP(find(isinf(p2_n_SP)))=[];

p3_f_SP(find(isnan(p3_f_SP)))=[];
p3_f_SP(find(isinf(p3_f_SP)))=[];

p3_n_SP(find(isnan(p3_n_SP)))=[];
p3_n_SP(find(isinf(p3_n_SP)))=[];

figure;hold on;
plot(1,p1_f_SP','o')
plot(2,p1_n_SP','o')
plot(4,p2_f_SP','o')
plot(5,p2_n_SP','o')
plot(7,p3_f_SP','o')
plot(8,p3_n_SP','o')

figure
bar([mean(p1_f_SP) mean(p1_n_SP) mean(p2_f_SP) mean(p2_n_SP) mean(p3_f_SP) mean(p3_n_SP)])



%plane_name={'plain1','plain2','plain3'};
% for p=1:3
% f_distribute=zeros(size(behavior_filepaths,2),length(range));
% f_distribute_norm=zeros(size(behavior_filepaths,2),length(range));
% n_distribute=zeros(size(behavior_filepaths,2),length(range));
% n_distribute_norm=zeros(size(behavior_filepaths,2),length(range));
% combine_distribute=zeros(size(behavior_filepaths,2),length(range));
% combine_distribute_norm=zeros(size(behavior_filepaths,2),length(range));
% 
% for f=1:size(behavior_filepaths,2)
%  load([temp behavior_filepaths{f}]);  
%  for i=1:size(PF_info,2)
%      if contains(PF_info(i).filepaths,'_f_')
%       
%         if contains(PF_info(i).filepaths,'late15lap')
%           
%           if contains(PF_info(i).filepaths,plane_name{p}) & ~isempty(PF_info(i).delay_lap_distribute)
%            f_distribute(f,:)=PF_info(i).delay_lap_distribute;
%            f_distribute_norm(f,:)=PF_info(i).delay_lap_distribute/sum(PF_info(i).delay_lap_distribute);
%               
%           end
%         end
%         
%      elseif contains(PF_info(i).filepaths,'_n_')
%       
%         if contains(PF_info(i).filepaths,'late15lap')
%           
%           if contains(PF_info(i).filepaths,plane_name{p}) & ~isempty(PF_info(i).delay_lap_distribute)
%            n_distribute(f,:)=PF_info(i).delay_lap_distribute;
%            n_distribute_norm(f,:)=PF_info(i).delay_lap_distribute/sum(PF_info(i).delay_lap_distribute);
%               
%           end
%           
%         elseif contains(PF_info(i).filepaths,'2nd25lap')
%           
%           if contains(PF_info(i).filepaths,plane_name{p}) & ~isempty(PF_info(i).delay_lap_distribute)
%            combine_distribute(f,:)=PF_info(i).delay_lap_distribute;
%            combine_distribute_norm(f,:)=PF_info(i).delay_lap_distribute/sum(PF_info(i).delay_lap_distribute);
%               
%           end
%  
%         end
%         
%      end 
%      
%      
%  end
%     
% end
% F_d{p}=f_distribute;
% F_d_norm{p}=f_distribute_norm;
% N_d{p}=n_distribute;
% N_d_norm{p}=n_distribute_norm;
% Comb_d{p}=combine_distribute;
% Comb_d_norm{p}=combine_distribute_norm;
% 
% subplot(3,3,p)
% bar(f_distribute')
% title(['familiar plane' num2str(p)]);
% ylim([0 25]);
% 
% subplot(3,3,3+p)
% bar(n_distribute')
% title(['novel plane' num2str(p)]);
% ylim([0 25]);
% 
% subplot(3,3,6+p)
% bar(combine_distribute')
% title(['combine plane' num2str(p)]);
% ylim([0 25]);
% 
% end



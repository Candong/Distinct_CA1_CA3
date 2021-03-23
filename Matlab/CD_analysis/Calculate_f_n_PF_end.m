clear all; %close all;

%%
[PC_filepaths, temp]=sort(uigetfile('*.mat', 'Chose f & n files to load:','MultiSelect','on'));
planenum = input(['What are the planes you would like to analyze  >> ']);
summary_cor=[];
inst_summary=[];

f_count=0;
n_count=0;
window=6;
threshold=3;
max_lap=50;
unpass_count=0;
f_multi_PF=0;
n_multi_PF=0;

for f=1:size(PC_filepaths,2)
 if contains(PC_filepaths{f},'_f_')
        f_count=f_count+1;
    for p=planenum
        cur_delay_PF=[];
        cur_inst_PF=[];
    if contains(PC_filepaths{f},['plain' num2str(p)]) 

        load(PC_filepaths{f});        
        f_sig_PFs{p,f_count}=sig_PFs;
        f_sig_PFs_with_noise{p,f_count}=sig_PFs_with_noise;
        f_meantrans_withnoise{p,f_count}=mean_trans;
        f_PF_start_bins{p,f_count}=PF_start_bins;
        f_PF_end_bins{p,f_count}=PF_end_bins;
        
        
        [pf_num pf_id]=find(cellfun(@isempty,sig_PFs)==0);  
        if any(pf_num>1)
            f_multi_PF=f_multi_PF+1;
        end
            
            
            
        pf_id=unique(pf_id);

        
            
    end
    end
 elseif contains(PC_filepaths{f},'_n_')
        n_count=n_count+1;
    for p=planenum
        cur_delay_PF=[];
        cur_inst_PF=[];
    if contains(PC_filepaths{f},['plain' num2str(p)]) 
        load(PC_filepaths{f});
        
        n_sig_PFs{p,n_count}=sig_PFs;
        n_sig_PFs_with_noise{p,n_count}=sig_PFs_with_noise;
        n_meantrans_withnoise{p,n_count}=mean_trans;
        n_PF_start_bins{p,n_count}=PF_start_bins;
        n_PF_end_bins{p,n_count}=PF_end_bins;
        
        [pf_num pf_id]=find(cellfun(@isempty,sig_PFs)==0);
        
        
        if pf_num>1
            n_multi_PF=n_multi_PF+1;
        end
        
        pf_id=unique(pf_id);

            
    end  
  end  
    
 end

end
%% 
 
for i=1:f_count
    
for p=planenum
f_COM_by_plane=[]; 
n_COM_by_plane=[];
%   for i=1:size(f_sig_PFs{p},2)
if ~isempty(f_sig_PFs{p,i})
[f_pf_num{p,i} f_pf_id{p,i}]=find(cellfun(@isempty,f_sig_PFs{p,i})==0);
[n_pf_num{p,i} n_pf_id{p,i}]=find(cellfun(@isempty,n_sig_PFs{p,i})==0);




if ~isempty(find(f_pf_num{p,i}>1))
    multi_id=find(f_pf_num{p,i}>1);
    f_multipf_id{p,i}=f_pf_id{p,i}(multi_id);
    f_pf_id{p,i}([multi_id; multi_id-1])=[];
    
      
end

if ~isempty(find(n_pf_num{p,i}>1))
    multi_id=find(n_pf_num{p,i}>1);
    n_multipf_id{p,i}=n_pf_id{p,i}(multi_id);
    n_pf_id{p,i}([multi_id; multi_id-1])=[];
      
end


end

    
end

end


%%
range=-1.01:0.02:1.01;


%%
% NOTICE only used the first and last 5 laps that have activities. amount
% of shifting might be different.
numbins=50;
all_bins=false;
ave_bin=5;

 
%%
stepsize=5;
standardlap=12;
stoplap=25;
onsetlap=0;
find_start=0;

[f_x f_y f_x_M f_y_M]=linreg_for_bylap_PF_end(f_sig_PFs,f_pf_id,stepsize,standardlap,onsetlap,stoplap,find_start);
[n_x n_y n_x_M n_y_M]=linreg_for_bylap_PF_end(n_sig_PFs,n_pf_id,stepsize,standardlap,onsetlap,stoplap,find_start);

%[f_x f_y f_jump_y f_jump_x]=linreg_for_bylap_groupedCOM(f_sig_PFs,f_pf_id,f_norm_act_COMshift,f_startlap,stepsize,standardlap,stoplap);
%[n_x n_y n_jump_y n_jump_x]=linreg_for_bylap_groupedCOM(n_sig_PFs,n_pf_id,n_norm_act_COMshift,n_startlap,stepsize,standardlap,stoplap);
 lap_range=1:1:max(f_x);
[count id]=histc(f_x,lap_range);
mean_fy=[];
f_sem=[];
for i=lap_range
    cur_fy=f_y(id==i);
    mean_fy=[mean_fy mean(cur_fy)];
    f_sem= [f_sem std(cur_fy)/sqrt(length(cur_fy))];
end
figure;
%plot(range,mean_y,'o')
errorbar(lap_range,mean_fy,f_sem,'o')

 lap_range=1:1:max(n_x);
[count id]=histc(n_x,lap_range);
mean_ny=[];
n_sem=[];
for i=lap_range
    cur_ny=n_y(id==i);
    mean_ny=[mean_ny mean(cur_ny)];
    n_sem= [n_sem std(cur_ny)/sqrt(length(cur_ny))];
end
hold on;
%plot(range,mean_y,'o')
errorbar(lap_range,mean_ny,n_sem,'o')
legend({'f','n'})
%legend({'nday1','nday2'})
%ylim([-2 2])
title('CA1 f and n')
%title('CA3 n day1 and day2')
save(['CA3_novel2day_delta_end_step_' num2str(stepsize) '_onsset_' num2str(onsetlap)],'f_x','f_y','n_x','n_y','f_x_M','f_y_M','n_x_M','n_y_M');
 %save(['CA3_fnovel_deltaCOM_' num2str(stoplap)],'f_x','f_y','n_x','n_y','f_jump_y','n_jump_y',...
 % 'f_jump_x','n_jump_x');

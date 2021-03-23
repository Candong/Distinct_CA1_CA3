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
        for i=1:length(pf_id)
        start=PF_start_bins(1,pf_id(i));
        PFend=PF_end_bins(1,pf_id(i));
        binmean=sig_PFs{1,pf_id(i)}; 
        max_lap=size(binmean,2);
        lap=find_delaylap(binmean,start,PFend,window,threshold,max_lap) ;
        if lap>1
            cur_delay_PF=[cur_delay_PF pf_id(i)];
        elseif lap==1
            cur_inst_PF=[cur_inst_PF pf_id(i)];
        elseif lap<1
            unpass_count=unpass_count+1;
        end
        end
        f_delay_id{p,f_count}=cur_delay_PF;
        f_inst_id{p,f_count}=cur_inst_PF;
        
            
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
%         for i=1:length(pf_id)
%         start=PF_start_bins(1,pf_id(i));
%         PFend=PF_end_bins(1,pf_id(i));
%         binmean=sig_PFs{1,pf_id(i)};           
%         lap=find_delaylap(binmean,start,PFend,window,threshold,max_lap) ;
%         if lap>1
%             cur_delay_PF=[cur_delay_PF pf_id(i)];
%         elseif lap==1
%             cur_inst_PF=[cur_inst_PF pf_id(i)];
%         end
%         end
%         n_delay_id{p,n_count}=cur_delay_PF;
%         n_inst_id{p,n_count}=cur_inst_PF;
            
    end  
  end  
    
 end

end
%% 
total_cor=[];
exclude_total_cor=[];
common_total_cor=[];
fn_total_cor=[];
inst_total_cor=[];
n_delay_total_cor=[];
f_inst_total_cor=[];
f_delay_total_cor=[];
total_f_mean=[];
total_n_mean=[];
f_total_COM=[];
n_total_COM=[];
f_total_SP=[];
n_total_SP=[];
binM=1:50;

f_total_neuron=[];
%n_total_neuron=[];
f_pf_num_summary=[];
n_pf_num_summary=[];
f_multipf_num_summary=[];
n_multipf_num_summary=[];
common_PC_num=[];

f_start_mean=[];
f_end_mean=[];
n_start_mean=[];
n_end_mean=[];
    


for i=1:f_count
    
for p=planenum
f_COM_by_plane=[]; 
n_COM_by_plane=[];
%   for i=1:size(f_sig_PFs{p},2)
if ~isempty(f_sig_PFs{p,i})
[f_pf_num{p,i} f_pf_id{p,i}]=find(cellfun(@isempty,f_sig_PFs{p,i})==0);
[n_pf_num{p,i} n_pf_id{p,i}]=find(cellfun(@isempty,n_sig_PFs{p,i})==0);

f_total_neuron=[f_total_neuron size(f_sig_PFs{p,i},2)];
%n_total_neuron=[n_total_neuron size(f_sig_PFs{p,i},2)];
f_multipf_num_summary=[f_multipf_num_summary sum(f_pf_num{p,i}>1)];
n_multipf_num_summary=[n_multipf_num_summary sum(n_pf_num{p,i}>1)];



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

f_pf_num_summary=[f_pf_num_summary size(f_pf_id{p,i},1)];
n_pf_num_summary=[n_pf_num_summary size(n_pf_id{p,i},1)];



cur_f_mean=f_meantrans_withnoise{p,i};
cur_n_mean=n_meantrans_withnoise{p,i};
total_f_mean=[total_f_mean; cur_f_mean];
total_n_mean=[total_n_mean; cur_n_mean];

cur_f_pfid=f_pf_id{p,i};
cur_n_pfid=n_pf_id{p,i};

common_PC_id{p,i}=intersect(f_pf_id{p,i},n_pf_id{p,i});
common_PC_num=[common_PC_num size(common_PC_id{p,i},1)];
if ~isempty(cur_f_pfid)
cur_f_COM=[];
    
    for j=1:length(cur_f_pfid)
        f_cur_binmean=f_sig_PFs{p,i}{1,cur_f_pfid(j)};
        %f_cur_binmean2=f_cur_binmean(:,21:25);
        %n_cur_binmean=n_sig_PFs{p,i}{1,common_PC_id{p,i}(j)};
        [f_meanCOM f_SP]=meanCOMandSP(f_cur_binmean',1,50,50);
        %[n_meanCOM n_SP]=meanCOMandSP(n_cur_binmean',1,50,50);
        f_total_COM=[f_total_COM f_meanCOM];
        cur_f_COM=[cur_f_COM f_meanCOM];
        %n_total_COM=[n_total_COM n_meanCOM];
        f_total_SP=[f_total_SP f_SP];
        %n_total_SP=[n_total_SP n_SP];
        
        f_cur_com=sum(f_cur_binmean'.*binM,2)./sum(f_cur_binmean',2);
        %n_cur_com=sum(n_cur_binmean'.*binM,2)./sum(n_cur_binmean',2);
        
        f_COM_by_plane=[f_COM_by_plane f_cur_com];
        %n_COM_by_plane=[n_COM_by_plane n_cur_com];
        cur_binmean_noise=f_sig_PFs_with_noise{p,i}{1,cur_f_pfid(j)};
        cur_start_mean=mean(cur_binmean_noise(:,1:10),2);
        cur_end_mean=mean(cur_binmean_noise(:,end-9:end),2);
        f_start_mean=[f_start_mean;cur_start_mean'];
        f_end_mean=[f_end_mean;cur_end_mean'];
        
        
    end
f_COM{p,i}=cur_f_COM;
end

if ~isempty(cur_n_pfid)
    cur_n_COM=[];
    
    for j=1:length(cur_n_pfid)
        %f_cur_binmean=f_sig_PFs{p,i}{1,cur_f_pfid(j)};
        n_cur_binmean=n_sig_PFs{p,i}{1,cur_n_pfid(j)};
        %n_cur_binmean2=n_cur_binmean(:,21:25);
        %[f_meanCOM f_SP]=meanCOMandSP(f_cur_binmean',1,50,50);
        [n_meanCOM n_SP]=meanCOMandSP(n_cur_binmean',1,50,50);
        %f_total_COM=[f_total_COM f_meanCOM];
        n_total_COM=[n_total_COM n_meanCOM];
        cur_n_COM=[cur_n_COM n_meanCOM];
        %f_total_SP=[f_total_SP f_SP];
        n_total_SP=[n_total_SP n_SP];
        
        %f_cur_com=sum(f_cur_binmean'.*binM,2)./sum(f_cur_binmean',2);
        n_cur_com=sum(n_cur_binmean'.*binM,2)./sum(n_cur_binmean',2);
        
        %f_COM_by_plane=[f_COM_by_plane f_cur_com];
        n_COM_by_plane=[n_COM_by_plane n_cur_com];
        cur_binmean_noise=n_sig_PFs_with_noise{p,i}{1,cur_n_pfid(j)};
        cur_start_mean=mean(cur_binmean_noise(:,1:10),2);
        cur_end_mean=mean(cur_binmean_noise(:,end-9:end),2);
        n_start_mean=[n_start_mean;cur_start_mean'];
        n_end_mean=[n_end_mean;cur_end_mean'];        
        
    end
  n_COM{p,i}=cur_n_COM;
end



    
end
f_animal_COM{p,i}=f_COM_by_plane;
n_animal_COM{p,i}=n_COM_by_plane;
end
%summary_cor=[summary_cor fn_total_cor];
%inst_summary=[inst_summary inst_total_cor];
end



remove_id=find(mean(f_start_mean,2)>5);
f_start_mean(remove_id,:)=[];
f_end_mean(remove_id,:)=[];

pop_cor_bybin=corr(f_start_mean,f_end_mean);


%%
range=-1.01:0.02:1.01;


%%
% NOTICE only used the first and last 5 laps that have activities. amount
% of shifting might be different.
numbins=50;
all_bins=false;
ave_bin=5;
[f_meanCOM_first5 f_meanCOM_last5  f_Act_lapnum f_norm_act_COMshift f_norm_from_startCOMshift f_startlap]=find_firstlastCOM(f_sig_PFs,f_pf_id,f_PF_start_bins,f_PF_end_bins,numbins,all_bins,ave_bin);
[n_meanCOM_first5 n_meanCOM_last5  n_Act_lapnum n_norm_act_COMshift n_norm_from_startCOMshift n_startlap]=find_firstlastCOM(n_sig_PFs,n_pf_id,n_PF_start_bins,n_PF_end_bins,numbins,all_bins,ave_bin);

[day1_meanCOM_first5 day1_meanCOM_last5  day1_Act_lapnum day1_norm_act_COMshift day1_norm_from_startCOMshift day1_startlap]=find_firstlastCOM(f_sig_PFs,common_PC_id,f_PF_start_bins,f_PF_end_bins,numbins,all_bins,ave_bin);
[day2_meanCOM_first5 day2_meanCOM_last5  day2_Act_lapnum day2_norm_act_COMshift day2_norm_from_startCOMshift day2_startlap]=find_firstlastCOM(n_sig_PFs,common_PC_id,n_PF_start_bins,n_PF_end_bins,numbins,all_bins,ave_bin);



total_f_shift=[];
total_n_shift=[];
total_twoday_shift=[];
total_f_norm_fromstart_shift=[];
total_n_norm_fromstart_shift=[];
total_f_act_shift=[];
total_n_act_shift=[];

common_day1_shift=[];
common_day2_shift=[];

total_f_start_lap=[];
total_n_start_lap=[];
total_day1_startlap=[];
total_day2_startlap=[];

for i =1:size(f_sig_PFs,2)
    for p=1:size(f_sig_PFs,1)
        f_shift=f_meanCOM_first5{p,i}-f_meanCOM_last5{p,i};
        n_shift=n_meanCOM_first5{p,i}-n_meanCOM_last5{p,i};
%         f_norm_fs=f_norm_from_startCOMshift{p,i};
%         n_norm_fs=n_norm_from_startCOMshift{p,i};
        f_act_shift=f_norm_act_COMshift{p,i};
        n_act_shift=n_norm_act_COMshift{p,i};
        twoday_shift=day1_meanCOM_last5{p,i}-day2_meanCOM_first5{p,i};
        day1_shift= day1_norm_act_COMshift{p,i};
        day2_shift= day2_norm_act_COMshift{p,i};
%         if ~isempty(abs(f_act_shift)>10)
%             remap_id=(abs(f_act_shift)>10);
%             true_id=f_pf_id{p,i}(remap_id);
%             for j=1:size(true_id)
%                 f_part=f_sig_PFs_with_noise{p,i}(1,true_id(j));
%                 n_part=n_sig_PFs_with_noise{p,i}(1,true_id(j));
%                 fd_part=f_sig_PFs{p,i}(1,true_id(j));
%                 nd_part=n_sig_PFs{p,i}(1,true_id(j));
%                 figure;
%                 subplot(1,2,1)
%                 imagesc([f_part{:}']);
%                 %imagesc([f_part{:}'; n_part{:}']);
%                 title([num2str(p) ' ' num2str(i) ' ' num2str(true_id(j))])
%                 subplot(1,2,2)
%                 imagesc([fd_part{:}']);
%                 %imagesc([fd_part{:}'; nd_part{:}']);
%                 title([num2str(sum(cellfun(@isempty,f_sig_PFs{p,i}(:,true_id(j)))==0)) ' ' num2str(sum(cellfun(@isempty,n_sig_PFs{p,i}(:,true_id(j)))==0))] );
%             end
%             
%             
%         end
        
        
%         if ~isempty(abs(n_act_shift)>10)
%             remap_id=(abs(n_act_shift)>10);
%             true_id=n_pf_id{p,i}(remap_id);
%             for j=1:size(true_id)
%                 f_part=f_sig_PFs_with_noise{p,i}(1,true_id(j));
%                 n_part=n_sig_PFs_with_noise{p,i}(1,true_id(j));
%                 fd_part=f_sig_PFs{p,i}(1,true_id(j));
%                 nd_part=n_sig_PFs{p,i}(1,true_id(j));
%                 figure;
%                 subplot(1,2,1)
%                 imagesc([n_part{:}']);
%                 %imagesc([f_part{:}'; n_part{:}']);
%                 title([num2str(p) ' ' num2str(i) ' ' num2str(true_id(j))])
%                 subplot(1,2,2)
%                 imagesc([nd_part{:}']);
%                 %imagesc([fd_part{:}'; nd_part{:}']);
%                 title([num2str(sum(cellfun(@isempty,f_sig_PFs{p,i}(:,true_id(j)))==0)) ' ' num2str(sum(cellfun(@isempty,n_sig_PFs{p,i}(:,true_id(j)))==0))] );
%             end
%             
%             
%         end
        
        
        total_f_shift=[total_f_shift f_shift];
        total_n_shift=[total_n_shift n_shift];
        total_twoday_shift=[total_twoday_shift twoday_shift];
%         total_f_norm_fromstart_shift=[total_f_norm_fromstart_shift f_norm_fs];
%         total_n_norm_fromstart_shift=[total_n_norm_fromstart_shift n_norm_fs];
        total_f_act_shift=[total_f_act_shift f_act_shift];
        total_n_act_shift=[total_n_act_shift n_act_shift];
        common_day1_shift=[common_day1_shift day1_shift];
        common_day2_shift=[common_day2_shift day2_shift];
        total_f_start_lap=[total_f_start_lap f_startlap{p,i}];
        total_n_start_lap=[total_n_start_lap n_startlap{p,i}];
        total_day1_startlap=[total_day1_startlap day1_startlap{p,i}];
        total_day2_startlap=[total_day2_startlap day2_startlap{p,i}];
    end
end

[f_bincounts,f_ind] = histc(total_f_shift,range);
[n_bincounts,n_ind] = histc(total_n_shift,range);
%[twoday_bincounts,f_ind] = histc(total_twoday_shift,range);

%bar(binranges,bincounts,'histc')
f_notout=find(~isoutlier(total_f_act_shift,'mean'));
n_notout=find(~isoutlier(total_n_act_shift,'mean'));
%% Shuffle the data
shuffle_num=1;
[f_shift_dist_M f_start5 f_end5]=find_shuffle_dist(shuffle_num,f_sig_PFs,f_pf_id,f_PF_start_bins,f_PF_end_bins,numbins,all_bins,ave_bin);
[day1_shift_dist_M day1_start5 day1_end5]=find_shuffle_dist(shuffle_num,f_sig_PFs,common_PC_id,f_PF_start_bins,f_PF_end_bins,numbins,all_bins,ave_bin);


f_total_h=[];f_total_p=[];
for i=1:shuffle_num
  [f_h,f_p]=kstest2(f_shift_dist_M(i,:),total_f_act_shift,'Tail','larger') ; 
  f_total_h=[f_total_h f_h];
  f_total_p=[f_total_p f_p];
  
end

%%

[n_shift_dist_M n_start5 n_end5]=find_shuffle_dist(shuffle_num,n_sig_PFs,n_pf_id,n_PF_start_bins,n_PF_end_bins,numbins,all_bins,ave_bin);
[day2_shift_dist_M day2_start5 day2_end5]=find_shuffle_dist(shuffle_num,n_sig_PFs,common_PC_id,n_PF_start_bins,n_PF_end_bins,numbins,all_bins,ave_bin);


n_total_h=[];n_total_p=[];
for i=1:shuffle_num
  [n_h,n_p]=kstest2(n_shift_dist_M(i,:),total_n_act_shift,'Tail','larger') ; 
  n_total_h=[n_total_h n_h];
  n_total_p=[n_total_p n_p];
  
end

%%
shuffle_2day_shift=day1_end5-day2_start5;
[shuffle_h shuffle_p]=kstest2(total_twoday_shift,shuffle_2day_shift)
 figure;%histogram(shuffle_2day_shift,-50:1:50,'Normalization','Probability')
 hold on;histogram(total_twoday_shift,-50:1:50,'Normalization','Probability')
 [sym_p sym_h]=signrank(total_twoday_shift)
 %legend('shuffle','across days')
 xlabel('bin')
 title(['across day shift.sign rank test p=' num2str(sym_p)])
[sym_p sym_h]=signrank(total_twoday_shift)



%%
f_instant_id=find(total_f_start_lap==1);
f_delay_id=find(total_f_start_lap>1);
n_instant_id=find(total_n_start_lap==1);
n_delay_id=find(total_n_start_lap>1);



%%
pc_label={'non-PC','f-PC','common-PC','n-PC'};
pf_label={'f-PF', 'f-multiPF','n-pf','n-multipf','common-pf'};
% figure;
% %subplot(1,2,1)
% title('PC summary')
% pie([sum(f_total_neuron)-sum(f_pf_num_summary)-sum(n_pf_num_summary)+sum(common_PC_num)...
%     sum(f_pf_num_summary)-sum(common_PC_num)+sum(f_multipf_num_summary)...
%      sum(common_PC_num) sum(n_pf_num_summary)-sum(common_PC_num)+sum(n_multipf_num_summary)])
%  legend(pc_label)
 
 delay_range=1:50;
 figure;
 subplot(1,2,1)
 histogram(total_day1_startlap,delay_range);
 title('f common pc')
 subplot(1,2,2)
 common_day1_dist=histc(total_day1_startlap,delay_range);
 f_dist=histc(total_f_start_lap,delay_range);
bar(f_dist-common_day1_dist);
title('f uncommon pc')
  figure;
  subplot(1,2,1)
 histogram(total_day2_startlap,delay_range);
 title('n common pc')
  subplot(1,2,2)
 common_day2_dist=histc(total_day2_startlap,delay_range);
 n_dist=histc(total_n_start_lap,delay_range);
bar(n_dist-common_day2_dist);
title('n uncommon pc')
 
figure;histogram(total_twoday_shift,-50:1:50) 
 
%  subplot(1,2,2)
% title('PF summary')
% num_uncommon=sum(abs(total_twoday_shift)>=10);
% pie([sum(f_pf_num_summary)+sum(common_PC_num)-num_uncommon sum(f_multipf_num_summary)*2 ... 
%      sum(n_pf_num_summary)+sum(common_PC_num)-num_uncommon sum(n_multipf_num_summary)*2 ...
%      sum(common_PC_num)-num_uncommon ])
%  legend(pf_label)
 
%%
stepsize=5;
standardlap=12;
stoplap=25;
[f_x f_y  f_x_M f_y_M]=linreg_for_bylap_groupedCOM_onsetlap(f_sig_PFs,f_pf_id,f_norm_act_COMshift,f_startlap,stepsize,standardlap,stoplap);
[n_x n_y  n_x_M n_y_M]=linreg_for_bylap_groupedCOM_onsetlap(n_sig_PFs,n_pf_id,n_norm_act_COMshift,n_startlap,stepsize,standardlap,stoplap);
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
ylim([-2 2])
title('CA1 nand f')

%   save(['new_CA3_novel2day_onsetlap_deltaCOM_' num2str(stoplap)],'f_x','f_y','n_x','n_y',...
%    'f_x_M','f_y_M','n_x_M','n_y_M');

%  lap_range=unique(f_jump_x);
% [count id]=histc(f_jump_x,lap_range);
% mean_fjumpy=[];
% f_jump_sem=[];
% for i=lap_range
%     cur_fy=f_jump_y(f_jump_x==i);
%     mean_fjumpy=[mean_fjumpy mean(cur_fy)];
%     f_jump_sem= [f_jump_sem std(cur_fy)/sqrt(length(cur_fy))];
% end
% figure;
% %plot(range,mean_y,'o')
% errorbar(lap_range,mean_fjumpy,f_jump_sem,'o')
% 
%  lap_range=unique(n_jump_x);
% [count id]=histc(n_jump_x,lap_range);
% mean_njumpy=[];
% n_jump_sem=[];
% for i=lap_range
%     cur_ny=n_jump_y(n_jump_x==i);
%     mean_njumpy=[mean_njumpy mean(cur_ny)];
%     n_jump_sem= [n_jump_sem std(cur_ny)/sqrt(length(cur_ny))];
% end
% %plot(range,mean_y,'o')
% hold on;
% errorbar(lap_range,mean_njumpy,n_jump_sem,'o')


%%
% [f_model f_x f_y]=linreg_for_allCOM(f_animal_COM,f_norm_act_COMshift,f_COM,f_startlap);
% [n_model n_x n_y]=linreg_for_allCOM(n_animal_COM,n_norm_act_COMshift,n_COM,n_startlap);
% 
%  lap_range=1:1:max(f_x);
% [count id]=histc(f_x,lap_range);
% mean_fy=[];
% f_sem=[];
% for i=lap_range
%     cur_fy=f_y(f_x==i);
%     mean_fy=[mean_fy mean(cur_fy)];
%     f_sem= [f_sem std(cur_fy)/sqrt(length(cur_fy))];
% end
% figure;
% %plot(range,mean_y,'o')
% errorbar(lap_range,mean_fy,f_sem,'o')
% 
%  lap_range=1:1:max(n_x);
% [count id]=histc(n_x,lap_range);
% mean_ny=[];
% n_sem=[];
% for i=lap_range
%     cur_ny=n_y(n_x==i);
%     mean_ny=[mean_ny mean(cur_ny)];
%     n_sem= [n_sem std(cur_ny)/sqrt(length(cur_ny))];
% end
% hold on;
% %plot(range,mean_y,'o')
% errorbar(lap_range,mean_ny,n_sem,'o')
% legend({'f','n'})
%ylim([-3 2])
%%
% figure;
% subplot(2,2,1); 
% bar(range,f_bincounts/length(common_total_cor),'histc');title('f shift')
% xlim([-50 50]);
% hold on; line([0 0],[0 1],'Color','r');
% subplot(2,2,2);
% bar(range,n_bincounts/length(common_total_cor),'histc');title('n shift')
% xlim([-50 50]);
% hold on; line([0 0],[0 1],'Color','r');
% subplot(2,2,3);
% bar(range,twoday_bincounts/length(common_total_cor),'histc');title('twaoday shift')
% xlim([-50 50]);
% hold on;line([0 0],[0 1],'Color','r');
% subplot(2,2,4);
% bar(range,f_bincounts/length(common_total_cor),'histc');
% hold on;bar(range,n_bincounts/length(common_total_cor),'histc');
% bar(range,twoday_bincounts/length(common_total_cor),'histc');
% xlim([-50 50]);
% legend({'f shift','n shift','twoday shift'})
% line([0 0],[0 1],'Color','r');
% set(gcf, 'Position',  [100, 100, 1000, 800])

% f_p=f_bincounts/length(common_total_cor);
% n_p=n_bincounts/length(common_total_cor);
% %twoday_p=twoday_bincounts/length(common_total_cor);
% figure;plot(range,cumsum(f_p));
% hold on;plot(range,cumsum(n_p));
%plot(range,cumsum(twoday_p));
%legend({'f shift','n shift','twoday shift'})
%xlim([-50 50])


% 
% figure;
% subplot(2,2,1);histogram(total_f_shift,range); hold on; %line([0 0],[0 150],'Color','r');
% % histfit(total_f_shift);
% % hold on;
% % histfit(total_n_shift);
% % histfit(total_twoday_shift);
% subplot(2,2,2);histogram(total_n_shift,range,'FaceColor','m');hold on; %line([0 0],[0 150],'Color','r');
% % histfit(total_n_shift);
% %subplot(2,2,3);histogram(total_twoday_shift,range,'FaceColor','y');hold on; %line([0 0],[0 150],'Color','r');
% % histfit(total_twoday_shift);
% subplot(2,2,4);
% histogram(total_f_shift,range)
% hold on;histogram(total_n_shift,range)
% %histogram(total_twoday_shift,range)
% legend({'f shift','n shift','twoday shift'})
% line([0 0],[0 150],'Color','r');
% set(gcf, 'Position',  [100, 100, 1000, 800])

% f = fit(range.',f_bincounts.','gauss1')
% figure;plot(f,range,f_bincounts)
% 
% pd = fitdist(total_f_shift','Normal');
% x_values = range;
% fit_f = pdf(pd,x_values);
% figure;plot(x_values,fit_f,'LineWidth',2)
% hold on;bar(range,f_bincounts/length(common_total_cor),'histc');title('f shift')

%%
%save('CA3_f_shift_new','f_shift_dist_M','n_shift_dist_M','total_f_act_shift','total_n_act_shift')



%%
% function [meanCOM_first5 meanCOM_last5 total_lapstart Act_lapnum norm_act_COMshift norm_from_startCOMshift]=find_firstlastCOM(sig_PFs,common_PC_id,PF_start_bins,PF_end_bins,numbins,all_bins)
% start_bin=1;
% end_bin=50;
% window=6;
% threshold=3;
% for i =1:size(sig_PFs,2)
%     for p=1:size(sig_PFs,1)
%         beginCOM=[];
%         endCOM=[];
%         lap_Start=[];
%         active_lap_num=[];
%         last_actlap=[];
%         if ~isempty(common_PC_id{p,i})
%             for j=1:size(common_PC_id{p,i},1)
%                 binM=1:numbins;
%                 cur_binmean=sig_PFs{p,i}{1,common_PC_id{p,i}(j)}';
%                 max_lap=size(cur_binmean,2);
%                 cur_act=find(~isnan(sum(cur_binmean,1)));
%                 last_act=cur_act(end);
%                 start_lap=find_delaylap(cur_binmean',start_bin,end_bin,window,threshold,max_lap);
%                 if start_lap>1
%                 cur_binmean(1:start_lap-1,:)=[];
%                 end
%                 cur_binmean(sum(cur_binmean,2)==0,:)=[];
%                 active_lap=size(cur_binmean,1);
%                 if all_bins
%                     start_b=1;
%                     end_b=numbins;
%                     [begin_meanCOM SP]=meanCOMandSP(cur_binmean(1:5,:),start_b,end_b,numbins);
%                     [end_meanCOM SP]=meanCOMandSP(cur_binmean(end-5:end,:),start_b,end_b,numbins);
%                 else
%                     start_b=PF_start_bins{p,i}(1,common_PC_id{p,i}(j));
%                     end_b=PF_end_bins{p,i}(1,common_PC_id{p,i}(j));
%                     [begin_meanCOM SP]=meanCOMandSP(cur_binmean(1:5,:),start_b,end_b,numbins);
%                     [end_meanCOM SP]=meanCOMandSP(cur_binmean(end-5:end,:),start_b,end_b,numbins);
%                 end
%                 beginCOM=[beginCOM begin_meanCOM];
%                 
%                 endCOM=[endCOM end_meanCOM];
%                 lap_Start=[lap_Start start_lap];
%                 active_lap_num=[active_lap_num active_lap];
%                 last_actlap=[last_actlap last_act];
%             end
%             meanCOM_first5{p,i}=beginCOM;
%             meanCOM_last5{p,i}=endCOM;
%             norm_COMshift{p,i}=(beginCOM-endCOM)/size(sig_PFs{p,i}{1,common_PC_id{p,i}(j)},2);
%             norm_from_startCOMshift{p,i}=(beginCOM-endCOM)./(size(sig_PFs{p,i}{1,common_PC_id{p,i}(j)},2)-lap_Start).*50;
%             norm_act_COMshift{p,i}=(beginCOM-endCOM)./(last_actlap-lap_Start).*50;
%             %norm_act_COMshift{p,i}=(beginCOM-endCOM)./active_lap_num.*50;
%             total_lapstart{p,i}=lap_Start;
%             Act_lapnum{p,i}=active_lap_num;
%         end
%     end
% end
% 
%  end 
clear all; %close all;

%%
[PC_filepaths, temp]=sort(uigetfile('*.mat', 'Chose f & n files to load:','MultiSelect','on'));
planenum = input(['What are the planes you would like to analysis  >> ']);
summary_cor=[];
inst_summary=[];

f_count=0;
n_count=0;
window=6;
threshold=3;
max_lap=25;
unpass_count=0;
f_lap_onset=[];
n_lap_onset=[];
for f=1:size(PC_filepaths,2)
 if contains(PC_filepaths{f},'_f_')
        f_count=f_count+1;
    for p=planenum
        cur_delay_PF=[];
        cur_inst_PF=[];
        f_lap=[];
    if contains(PC_filepaths{f},['plain' num2str(p)]) 

        load(PC_filepaths{f});  
        load([PC_filepaths{f}(1:end-26) 'binMatrix']); 
        f_sig_PFs{p,f_count}=sig_PFs;
        f_sig_PFs_with_noise{p,f_count}=sig_PFs_with_noise;
        f_meantrans_withnoise{p,f_count}=mean_trans;
        f_PF_start_bins{p,f_count}=PF_start_bins;
        f_PF_end_bins{p,f_count}=PF_end_bins;
        binMatrix_session1(isnan(binMatrix_session1))=0;
        f_binM{p,f_count}=binMatrix_session1;
        
        
        
        [pf_num pf_id]=find(cellfun(@isempty,sig_PFs)==0);
        pf_id=unique(pf_id);
        
        for i=1:length(pf_id)
        start=PF_start_bins(1,pf_id(i));
        PFend=PF_end_bins(1,pf_id(i));
        binmean=sig_PFs{1,pf_id(i)}; 
        max_lap=size(binmean,2);
        lap=find_delaylap(binmean,start,PFend,window,threshold,max_lap) ;
        f_lap=[f_lap lap];
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
    f_lap_onset=[f_lap_onset f_lap];
    end
 elseif contains(PC_filepaths{f},'_n_')
        n_count=n_count+1;
    for p=planenum
        cur_delay_PF=[];
        cur_inst_PF=[];
    if contains(PC_filepaths{f},['plain' num2str(p)]) 
        load(PC_filepaths{f});
        load([PC_filepaths{f}(1:end-26) 'binMatrix']); 
        
        n_sig_PFs{p,n_count}=sig_PFs;
        n_sig_PFs_with_noise{p,n_count}=sig_PFs_with_noise;
        n_meantrans_withnoise{p,n_count}=mean_trans;
        n_PF_start_bins{p,n_count}=PF_start_bins;
        n_PF_end_bins{p,n_count}=PF_end_bins;
        binMatrix_session1(isnan(binMatrix_session1))=0;
        n_binM{p,n_count}=binMatrix_session1;
        
        [pf_num pf_id]=find(cellfun(@isempty,sig_PFs)==0);
        pf_id=unique(pf_id);
        n_lap=[];
        for i=1:length(pf_id)
        start=PF_start_bins(1,pf_id(i));
        PFend=PF_end_bins(1,pf_id(i));
        binmean=sig_PFs{1,pf_id(i)};           
        lap=find_delaylap(binmean,start,PFend,window,threshold,max_lap) ;
        n_lap=[n_lap lap];
        if lap>1
            cur_delay_PF=[cur_delay_PF pf_id(i)];
        elseif lap==1
            cur_inst_PF=[cur_inst_PF pf_id(i)];
        end
        end
        n_delay_id{p,n_count}=cur_delay_PF;
        n_inst_id{p,n_count}=cur_inst_PF;
            
    end  
  n_lap_onset=[n_lap_onset n_lap];  
  end  
    
 end

end

%%
f_mean_total_cor=[];
f_last_total_cor=[];
f_n_total_cor=[];
total_COM=[];
total_SP=[];
total_f_mean=[];
total_n_mean=[];
for i=1:f_count
    
for p=planenum
  
%   for i=1:size(f_sig_PFs{p},2)
if ~isempty(f_sig_PFs{p,i})
%[f_pf_num{p,i} f_pf_id{p,i}]=unique(find(cellfun(@isempty,f_sig_PFs{p,i})==0));
%[n_pf_num{p,i} n_pf_id{p,i}]=unique(find(cellfun(@isempty,n_sig_PFs{p,i})==0));
[f_pf_num{p,i} f_pf_id{p,i}]=find(cellfun(@isempty,f_sig_PFs{p,i})==0);
[n_pf_num{p,i} n_pf_id{p,i}]=find(cellfun(@isempty,n_sig_PFs{p,i})==0);

f_pf_id{p,i}=unique(f_pf_id{p,i});
n_pf_id{p,i}=unique(n_pf_id{p,i});


cur_f_mean=f_meantrans_withnoise{p,i};
cur_n_mean=n_meantrans_withnoise{p,i};
cur_n_binM=n_binM{p,i};
cur_f_binM=f_binM{p,i};
f_n_PCid{p,i}=union(f_pf_id{p,i},n_pf_id{p,i});





 
% common_PC_id{p,i}=intersect(f_pf_id{p,i},n_pf_id{p,i});
% if ~isempty(common_PC_id{p,i})
% common_cor{p,i}=diag(corr(cur_f_mean(common_PC_id{p,i},:)',cur_n_mean(common_PC_id{p,i},:)'));
% common_total_cor=[common_total_cor common_cor{p,i}']; 
% end
if ~isempty(f_pf_id{p,i})
    f_last_lapM=cur_f_binM(:,end,f_pf_id{p,i});
    n_first_lapM=cur_n_binM(:,1,f_pf_id{p,i});
    first_lap_Fc3=reshape(n_first_lapM,[size(n_first_lapM,1), size(n_first_lapM,3)]);
    last_lap_Fc3=reshape(f_last_lapM,[size(f_last_lapM,1), size(f_last_lapM,3)]);
    f_pfmean_cor{p,i}=diag(corr(cur_f_mean(f_pf_id{p,i},:)',first_lap_Fc3));
    f_pflast_cor{p,i}=diag(corr(last_lap_Fc3,first_lap_Fc3));
    f_n_mean_cor{p,i}=diag(corr(cur_f_mean(f_pf_id{p,i},:)',cur_n_mean(f_pf_id{p,i},:)'));
    
    f_mean_total_cor=[f_mean_total_cor f_pfmean_cor{p,i}'];
    f_last_total_cor=[f_last_total_cor f_pflast_cor{p,i}'];
    f_n_total_cor=[f_n_total_cor f_n_mean_cor{p,i}'];
    
    %f_mean_cor_M{p,i}(:,:,
    total_f_mean=[total_f_mean; cur_f_mean(f_pf_id{p,i},:)];
    total_n_mean=[total_n_mean; cur_n_mean(f_pf_id{p,i},:)];
    
    cur_sig_PFs=f_sig_PFs{p,i};
    cur_COM=[];
    cur_SP=[];
    for n=f_pf_id{p,i}'
        
        cur_binmean=cur_sig_PFs{1,n};
        [meanCOM SP]=meanCOMandSP(cur_binmean',1,50,50);
        cur_COM=[cur_COM meanCOM];
        cur_SP=[cur_SP SP];   
    end
    total_COM=[total_COM cur_COM];
    total_SP=[total_SP cur_SP];
    
    
    
end    
end
end
%summary_cor=[summary_cor fn_total_cor];
%inst_summary=[inst_summary inst_total_cor];
end
f_mean_total_cor(isnan(f_mean_total_cor))=-1.05;
f_last_total_cor(isnan(f_last_total_cor))=-1.05;
f_n_total_cor(isnan(f_n_total_cor))=-1.05;
range=-1.1:0.1:1;
f_mean_dist=histc(f_mean_total_cor,range);
figure;histogram(f_mean_total_cor,range);
title('cor of day1 mean and day2 1st lap')


f_last_dist=histc(f_last_total_cor,range);
figure;histogram(f_last_total_cor,range);
title('cor of day1 last and day2 1st lap')


f_n_dist=histc(f_n_total_cor,range);
figure;histogram(f_n_total_cor,range);
title('cor of day1 mean and day2 mean lap')

%%
[sorted_COM, COM_id]=sort(total_COM);
sorted_f_n_total_cor=f_n_total_cor(COM_id);
sorted_SP=total_SP(COM_id);
bin_COM=ceil(sorted_COM);

mean_cor=[];
mean_SP=[];
for i=min(bin_COM):max(bin_COM)
    cur_meancor=mean(sorted_f_n_total_cor(bin_COM==i));
    mean_cor=[mean_cor cur_meancor];
    
    cur_meanSP=mean(sorted_SP(bin_COM==i));
    mean_SP=[mean_SP cur_meanSP];
    
end
figure;
plot(bin_COM,sorted_f_n_total_cor,'o')
hold on;plot(min(bin_COM):max(bin_COM),mean_cor);
title('cor by bin')

figure;
plot(bin_COM,sorted_SP,'o')
hold on;plot(min(bin_COM):max(bin_COM),mean_SP);
title('SP by bin')

high_fn_cor=sorted_f_n_total_cor(sorted_f_n_total_cor>=0.7);
high_bin_COM=bin_COM(sorted_f_n_total_cor>=0.7);
unique_COM=unique(high_bin_COM);
cor_dist=[];
for i=1:50
    cur_cor_num=sum(high_bin_COM==i);
    cor_dist=[cor_dist cur_cor_num];
  
end
figure;
bar(cor_dist)
%plot(high_bin_COM,high_fn_cor,'o')
%hold on;plot(min(bin_COM):max(bin_COM),mean_cor);
title('high_cor by bin')
%%
figure;
lap_range=-1:50;
high_cor_lap=f_lap_onset(f_n_total_cor>0.7);
f_dist=histc(f_lap_onset,lap_range);
subplot(1,2,1)
bar(f_dist);
title('all distribution')

highcor_f_dist=histc(high_cor_lap,lap_range);
subplot(1,2,2)
bar(highcor_f_dist)
title('high cor distribution');

n_dist=histc(n_lap_onset,lap_range);
figure;bar(n_dist)
%%
remove_id=find(total_f_mean(:,1)>10);
total_f_mean(remove_id,:)=[];
total_n_mean(remove_id,:)=[];
pop_cor_bybin=corr(total_f_mean,total_n_mean);
figure;
imagesc(flip(pop_cor_bybin));
title('pop corelation by bin')
hold on;plot(1:50,50:-1:1,'LineWidth',2)


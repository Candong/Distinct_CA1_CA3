clear all; close all;

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
f_trans_max=[];
n_trans_max=[];
f_frq=[];
n_frq=[];
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
        f_trans_max=[f_trans_max, max(binmean)];
        f_frq=[f_frq sum(max(binmean)~=0)/size(binmean,2)];
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
        for i=1:length(pf_id)
        start=PF_start_bins(1,pf_id(i));
        PFend=PF_end_bins(1,pf_id(i));
        binmean=sig_PFs{1,pf_id(i)};           
        lap=find_delaylap(binmean,start,PFend,window,threshold,max_lap) ;
        if lap>1
            cur_delay_PF=[cur_delay_PF pf_id(i)];
        elseif lap==1
            cur_inst_PF=[cur_inst_PF pf_id(i)];
        end
        n_trans_max=[n_trans_max, max(binmean)];
        n_frq=[n_frq sum(max(binmean)~=0)/size(binmean,2)];
        end
        n_delay_id{p,n_count}=cur_delay_PF;
        n_inst_id{p,n_count}=cur_inst_PF;
            
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



%%
range=-1.01:0.02:1.01;
% high_cor_id=find(common_total_cor>0.6);
% figure;histogram(f_total_COM(high_cor_id)-n_total_COM(high_cor_id),range);

%%
% NOTICE only used the first and last 5 laps that have activities. amount
% of shifting might be different.
numbins=50;
all_bins=false;
[f_meanCOM_first5 f_meanCOM_last5  f_Act_lapnum f_norm_act_COMshift f_norm_from_startCOMshift f_startlap]=find_firstlastCOM_character(f_sig_PFs,f_pf_id,f_PF_start_bins,f_PF_end_bins,numbins,all_bins,5);
[n_meanCOM_first5 n_meanCOM_last5  n_Act_lapnum n_norm_act_COMshift n_norm_from_startCOMshift n_startlap]=find_firstlastCOM_character(n_sig_PFs,n_pf_id,n_PF_start_bins,n_PF_end_bins,numbins,all_bins,5);

[day1_meanCOM_first5 day1_meanCOM_last5  day1_Act_lapnum day1_norm_act_COMshift day1_norm_from_startCOMshift day1_startlap]=find_firstlastCOM_character(f_sig_PFs,common_PC_id,f_PF_start_bins,f_PF_end_bins,numbins,all_bins,5);
[day2_meanCOM_first5 day2_meanCOM_last5  day2_Act_lapnum day2_norm_act_COMshift day2_norm_from_startCOMshift day2_startlap]=find_firstlastCOM_character(n_sig_PFs,common_PC_id,n_PF_start_bins,n_PF_end_bins,numbins,all_bins,5);






% end

%%
f_trans_max(f_trans_max==0)=[];
n_trans_max(n_trans_max==0)=[];
range=0:0.2:7;
figure;
histogram(f_trans_max,range,'Normalization','probability')
hold on
histogram(n_trans_max,range,'Normalization','probability')
legend({'f','n'})
title('f n trans peak')
xlabel('dF/F')
% legend({'nday1','nday2'})
% title('n 2 day trans peak')
box off
figure
cdfplot(f_trans_max);
hold on; cdfplot(n_trans_max);
%xlim([0 30])
[p h stats]=ranksum(f_trans_max,n_trans_max);
title(['fn transient peak wilcoxon test,p= ' num2str(p)]);
xlabel('dF/F')
legend({'f','n'})
xlim([0 7])
box off;grid off;

figure;
histogram(f_frq,0.15:0.05:1,'Normalization','probability')
hold on
histogram(n_frq,0.15:0.05:1,'Normalization','probability')
legend({'f','n'})
title('f n pf frq')
% legend({'nday1','nday2'})
% title('n 2 day frq')
box off
figure
cdfplot(f_frq);
hold on; cdfplot(n_frq);
%xlim([0 30])
[p h stats]=ranksum(f_frq,n_frq);
title(['fn pf frq wilcoxon test,p= ' num2str(p)]);
xlabel('percentage')
legend({'f','n'})
%xlim([0 10])
box off;grid off;

 %%
 function [meanCOM_first5 meanCOM_last5 Act_lapnum norm_act_COMshift norm_from_startCOMshift total_lapstart]=find_firstlastCOM_character(f_sig_PFs,common_PC_id,PF_start_bins,PF_end_bins,numbins,all_bins,ave_bin)
start_bin=1;
end_bin=50;
window=6;
threshold=3;
for i =1:size(f_sig_PFs,2)
    for p=1:size(f_sig_PFs,1)
        beginCOM=[];
        endCOM=[];
        lap_Start=[];
        active_lap_num=[];
        last_actlap=[];
        if ~isempty(common_PC_id{p,i})
            for j=1:size(common_PC_id{p,i},1)
                binM=1:numbins;
                cur_binmean=f_sig_PFs{p,i}{1,common_PC_id{p,i}(j)}';
                %n_cur_binmean=n_sig_PFs{p,i}{1,common_PC_id{p,i}(j)}';
                %if ~isclipped(cur_binmean') & ~isclipped(n_cur_binmean')
                %cur_binmean2=sig_PFs{p,i}{1,common_PC_id{p,i}(j)}';
                max_lap=size(cur_binmean,2);
                cur_act=find(~isnan(sum(cur_binmean,1)));
                last_act=cur_act(end);
                start_lap=find_delaylap(cur_binmean',start_bin,end_bin,window,threshold,max_lap);
                if start_lap>1
                cur_binmean(1:start_lap-1,:)=[];
                end
                cur_binmean(sum(cur_binmean,2)==0,:)=[];
                active_lap=size(cur_binmean,1);
                if all_bins
                    start_b=1;
                    end_b=numbins;
                    [begin_meanCOM SP]=meanCOMandSP(cur_binmean(1:ave_bin,:),start_b,end_b,numbins);
                    [end_meanCOM SP]=meanCOMandSP(cur_binmean(end-ave_bin+1:end,:),start_b,end_b,numbins);
                else
                    start_b=PF_start_bins{p,i}(1,common_PC_id{p,i}(j));
                    end_b=PF_end_bins{p,i}(1,common_PC_id{p,i}(j));
                    [begin_meanCOM SP]=meanCOMandSP(cur_binmean(1:ave_bin,:),start_b,end_b,numbins);
                    [end_meanCOM SP]=meanCOMandSP(cur_binmean(end-ave_bin+1:end,:),start_b,end_b,numbins);
                end
                beginCOM=[beginCOM begin_meanCOM];
                
                endCOM=[endCOM end_meanCOM];
                lap_Start=[lap_Start start_lap];
                active_lap_num=[active_lap_num active_lap];
                last_actlap=[last_actlap last_act];
                %end
            end
            meanCOM_first5{p,i}=beginCOM;
            meanCOM_last5{p,i}=endCOM;
            norm_COMshift{p,i}=(beginCOM-endCOM)/size(f_sig_PFs{p,i}{1,common_PC_id{p,i}(j)},2);
            norm_from_startCOMshift{p,i}=(beginCOM-endCOM)./(size(f_sig_PFs{p,i}{1,common_PC_id{p,i}(j)},2)-lap_Start);
            norm_act_COMshift{p,i}=(beginCOM-endCOM)./(last_actlap-lap_Start);
            %norm_act_COMshift{p,i}=(beginCOM-endCOM)./active_lap_num.*50;
            total_lapstart{p,i}=lap_Start;
            Act_lapnum{p,i}=active_lap_num;
            
        end
    end
end

 end 
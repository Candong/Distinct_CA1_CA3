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

if ~isa(PC_filepaths,'cell') 
  PC_filepaths={PC_filepaths};  
    
end 

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
        if pf_num>1
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
for i=1:f_count
    
for p=planenum
f_COM_by_plane=[]; 
n_COM_by_plane=[];
%   for i=1:size(f_sig_PFs{p},2)
if ~isempty(f_sig_PFs{p,i})
[f_pf_num{p,i} f_pf_id{p,i}]=find(cellfun(@isempty,f_sig_PFs{p,i})==0);
[n_pf_num{p,i} n_pf_id{p,i}]=find(cellfun(@isempty,n_sig_PFs{p,i})==0);

cur_f_mean=f_meantrans_withnoise{p,i};
cur_n_mean=n_meantrans_withnoise{p,i};
total_f_mean=[total_f_mean; cur_f_mean];
total_n_mean=[total_n_mean; cur_n_mean];
% if sum(mean(cur_f_mean,2)>5)
%     p
%     i
%     find(mean(cur_f_mean,2)>5)
% end
% f_n_PCid{p,i}=union(f_pf_id{p,i},n_pf_id{p,i});
% 
% transient_cor{p,i}=diag(corr(cur_f_mean(f_n_PCid{p,i},:)',cur_n_mean(f_n_PCid{p,i},:)'));
%       
% total_cor=[total_cor transient_cor{p,i}'];   
% 
% exclude_PC_id{p,i}=setxor(f_pf_id{p,i},n_pf_id{p,i});
% exclude_cor{p,i}=diag(corr(cur_f_mean(exclude_PC_id{p,i},:)',cur_n_mean(exclude_PC_id{p,i},:)'));
% exclude_total_cor=[exclude_total_cor exclude_cor{p,i}'];    


common_PC_id{p,i}=intersect(f_pf_id{p,i},n_pf_id{p,i});
if ~isempty(common_PC_id{p,i})
common_cor{p,i}=diag(corr(cur_f_mean(common_PC_id{p,i},:)',cur_n_mean(common_PC_id{p,i},:)'));
common_total_cor=[common_total_cor common_cor{p,i}']; 
    
    for j=1:length(common_PC_id{p,i})
        f_cur_binmean=f_sig_PFs{p,i}{1,common_PC_id{p,i}(j)};
        n_cur_binmean=n_sig_PFs{p,i}{1,common_PC_id{p,i}(j)};
        [f_meanCOM f_SP]=meanCOMandSP(f_cur_binmean',1,50,50);
        [n_meanCOM n_SP]=meanCOMandSP(n_cur_binmean',1,50,50);
        f_total_COM=[f_total_COM f_meanCOM];
        n_total_COM=[n_total_COM n_meanCOM];
        f_total_SP=[f_total_SP f_SP];
        n_total_SP=[n_total_SP n_SP];
        
        f_cur_com=sum(f_cur_binmean'.*binM,2)./sum(f_cur_binmean',2);
        n_cur_com=sum(n_cur_binmean'.*binM,2)./sum(n_cur_binmean',2);
        
        f_COM_by_plane=[f_COM_by_plane f_cur_com];
        n_COM_by_plane=[n_COM_by_plane n_cur_com];
        
        
    end

end



    
end
f_animal_COM{p,i}=f_COM_by_plane;
n_animal_COM{p,i}=n_COM_by_plane;
end
%summary_cor=[summary_cor fn_total_cor];
%inst_summary=[inst_summary inst_total_cor];
end
figure;
boxplot([f_total_SP; n_total_SP]')
%%
range=-50.5:1:50.5;
% high_cor_id=find(common_total_cor>0.6);
% figure;histogram(f_total_COM(high_cor_id)-n_total_COM(high_cor_id),range);
%%
% f_COMshift=nan(150,200);
% n_COMshift=nan(150,200);
% x_max=1;
% y=[];
% x=[];
% for i =1:size(f_animal_COM,2)
%     for p=1:size(f_animal_COM,1)
%         if ~isempty(f_animal_COM{p,i})
%         %x_max=max(x_max,size(f_animal_COM{p,i},1)-1);
%         cur_COM_M=f_animal_COM{p,i};
% %         cur_x=ones(size(cur_y)).*(1:size(cur_y,1))';
%         for j=1:size(f_animal_COM{p,i},2)
%             cur_y=cur_COM_M(:,j)';
%             cur_x=1:length(cur_y)';
%             remove_id=isnan(cur_y);
%             cur_y(remove_id)=[];
%             cur_y=cur_y-cur_y(1);           
%             cur_x(remove_id)=[];
%             cur_x=cur_x-cur_x(1);
%             y=[y cur_y];
%             x=[x cur_x];
%         end
%             
% %         remove_id=isnan(cur_y);
% %         cur_y(remove_id)=[];
% %         cur_y=cur_y-cur_y(1);
% %         cur_x(remove_id)=[];
% %         cur_x=cur_x-curx
% %         
% % %         for j=1:size(f_aniamal_COM{p,i},2)
% % %         cur_x=
% % %         end
% % %         
% % %         
% %          y=[y cur_y];
% %          x=[x cur_x];
%         end
%         
%         
%     end
% end
[f_model f_x f_y]=linreg_for_allCOM(f_animal_COM);
[n_model n_x n_y]=linreg_for_allCOM(n_animal_COM);
% figure;plot(f_x,f_y,'o')
% 
% %mdl=fitlm(x,y);
% hold on;
% a=f_model.('Coefficients').('Estimate')(2);
% plot(f_x,f_model.('Coefficients').('Estimate')(2).*f_x+ f_model.('Coefficients').('Estimate')(1));
% title(['slope=' num2str(a)]);

%%
%%
% NOTICE only used the first and last 5 laps that have activities. amount
% of shifting might be different.
numbins=50;
all_bins=false;
[f_meanCOM_first5 f_meanCOM_last5 f_total_lap_start f_Act_lapnum f_norm_act_COMshift f_norm_from_startCOMshift]=find_firstlastCOM(f_sig_PFs,common_PC_id,f_PF_start_bins,f_PF_end_bins,numbins,all_bins);
[n_meanCOM_first5 n_meanCOM_last5 n_total_lap_start n_Act_lapnum n_norm_act_COMshift n_norm_from_startCOMshift]=find_firstlastCOM(n_sig_PFs,common_PC_id,n_PF_start_bins,n_PF_end_bins,numbins,all_bins);

total_f_shift=[];
total_n_shift=[];
total_twoday_shift=[];
total_f_norm_fromstart_shift=[];
total_n_norm_fromstart_shift=[];
total_f_act_shift=[];
total_n_act_shift=[];

for i =1:size(f_sig_PFs,2)
    for p=1:size(f_sig_PFs,1)
        f_shift=f_meanCOM_first5{p,i}-f_meanCOM_last5{p,i};
        n_shift=n_meanCOM_first5{p,i}-n_meanCOM_last5{p,i};
        f_norm_fs=f_norm_from_startCOMshift{p,i};
        n_norm_fs=n_norm_from_startCOMshift{p,i};
        f_act_shift=f_norm_act_COMshift{p,i};
        n_act_shift=n_norm_act_COMshift{p,i};
        twoday_shift=f_meanCOM_last5{p,i}-n_meanCOM_first5{p,i};
        
        total_f_shift=[total_f_shift f_shift];
        total_n_shift=[total_n_shift n_shift];
        total_twoday_shift=[total_twoday_shift twoday_shift];
        total_f_norm_fromstart_shift=[total_f_norm_fromstart_shift f_norm_fs];
        total_n_norm_fromstart_shift=[total_n_norm_fromstart_shift n_norm_fs];
        total_f_act_shift=[total_f_act_shift f_act_shift];
        total_n_act_shift=[total_n_act_shift n_act_shift];
    
    end
end

[f_bincounts,f_ind] = histc(total_f_shift,range);
[n_bincounts,n_ind] = histc(total_n_shift,range);
[twoday_bincounts,f_ind] = histc(total_twoday_shift,range);

%bar(binranges,bincounts,'histc')
%%
figure;
subplot(2,2,1)
histogram(total_f_norm_fromstart_shift,range);
title('from start f')
subplot(2,2,2)
histogram(total_n_norm_fromstart_shift,range);
title('from start n')
subplot(2,2,3)
histogram(total_f_act_shift,range);
title('act f')
subplot(2,2,4)
histogram(total_n_act_shift,range);
title('act n')
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

f_p=f_bincounts/length(common_total_cor);
n_p=n_bincounts/length(common_total_cor);
twoday_p=twoday_bincounts/length(common_total_cor);
figure;plot(range,cumsum(f_p));
hold on;plot(range,cumsum(n_p));
plot(range,cumsum(twoday_p));
legend({'f shift','n shift','twoday shift'})
%xlim([-50 50])



figure;
subplot(2,2,1);histogram(total_f_shift,range); hold on; %line([0 0],[0 150],'Color','r');
% histfit(total_f_shift);
% hold on;
% histfit(total_n_shift);
% histfit(total_twoday_shift);
subplot(2,2,2);histogram(total_n_shift,range,'FaceColor','m');hold on; %line([0 0],[0 150],'Color','r');
% histfit(total_n_shift);
subplot(2,2,3);histogram(total_twoday_shift,range,'FaceColor','y');hold on; %line([0 0],[0 150],'Color','r');
% histfit(total_twoday_shift);
subplot(2,2,4);
histogram(total_f_shift,range)
hold on;histogram(total_n_shift,range)
histogram(total_twoday_shift,range)
legend({'f shift','n shift','twoday shift'})
line([0 0],[0 150],'Color','r');
set(gcf, 'Position',  [100, 100, 1000, 800])

% f = fit(range.',f_bincounts.','gauss1')
% figure;plot(f,range,f_bincounts)
% 
% pd = fitdist(total_f_shift','Normal');
% x_values = range;
% fit_f = pdf(pd,x_values);
% figure;plot(x_values,fit_f,'LineWidth',2)
% hold on;bar(range,f_bincounts/length(common_total_cor),'histc');title('f shift')

%%
save('CA1_shift','f_bincounts','n_bincounts','twoday_bincounts')


%%
% function [meanCOM_first5 meanCOM_last5 norm_COMshift]=find_firstlastCOM(sig_PFs,common_PC_id,PF_start_bins,PF_end_bins,numbins,all_bins)
% for i =1:size(sig_PFs,2)
%     for p=1:size(sig_PFs,1)
%         beginCOM=[];
%         endCOM=[];
%         if ~isempty(common_PC_id{p,i})
%             for j=1:size(common_PC_id{p,i},1)
%                 binM=1:numbins;
%                 cur_binmean=sig_PFs{p,i}{1,common_PC_id{p,i}(j)}';
%                 cur_binmean(sum(cur_binmean,2)==0,:)=[];
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
%                 endCOM=[endCOM end_meanCOM];
%             end
%             meanCOM_first5{p,i}=beginCOM;
%             meanCOM_last5{p,i}=endCOM;
%             norm_COMshift{p,i}=(beginCOM-endCOM)/size(sig_PFs{p,i}{1,common_PC_id{p,i}(j)},2);
%         end
%     end
% end
% 
%  end 
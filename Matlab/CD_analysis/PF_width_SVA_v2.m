clear all

[PC_filepaths, temp]=(uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on'));

if ~isa(PC_filepaths,'cell') 
  PC_filepaths={PC_filepaths};  
end 
%%
f_count = 0;
f_start_bin = [];
f_end_bin = [];
f_COM=[];
f_SP=[];
f_all_SP=[];
f_all_remove=[];
for f=1:size(PC_filepaths,2)
    
    if contains(PC_filepaths{f},'_f_')
        f_count = f_count + 1;
        if contains(PC_filepaths{f},['plain']) 
            load(PC_filepaths{f});      
        end
    f_start_bin = [f_start_bin,PF_start_bins];
    f_end_bin = [f_end_bin,PF_end_bins];
    COM_M=ones(size(PF_start_bins)).*(-1);
    SP_M=ones(size(PF_start_bins)).*(-1);
    for m=1:size(sig_PFs,1)
       for n=1:size(sig_PFs,2)
          if ~isempty(sig_PFs{m, n})
             cur_remove=0;
             [COM SP]=meanCOMandSP(sig_PFs{m, n}',1,50,50); 
             COM_M(m,n)=COM;
             SP_M(m,n)=SP;
             f_all_SP=[f_all_SP SP];
             meanPF=mean(sig_PFs{m, n}',1);
             if meanPF(1)~=0
                cur_remove=meanPF(1)<=(0.1*max(meanPF));
             elseif meanPF(end)~=0
                cur_remove=meanPF(end)<=(0.1*max(meanPF)); 

             end
             f_all_remove=[f_all_remove cur_remove];
                           
          end
           
       end
                
    end
    f_COM=[f_COM,COM_M];
    f_SP=[f_SP,SP_M];
    
    end


 
end

f_all_SP(logical(f_all_remove))=[];

for i = 1:length(f_start_bin(:,1))
    for j = 1:length(f_end_bin(1,:))
        PF_width_fam(i,j) = f_end_bin(i,j) - f_start_bin(i,j);
    end
end

%PF_width_fam(PF_width_fam == 0) = nan; 
% PF_w_f=reshape(PF_width_fam,1,[]);
% PF_w_f(isnan(PF_w_f))=[];
% f_SP_reshape=reshape(f_SP,1,[]);
% f_SP_reshape(f_SP_reshape==-1)=[];
% f_COM_reshape=reshape(f_COM,1,[]);
% f_COM_reshape(f_COM_reshape==-1)=[];
% f_COM_reshape=round(f_COM_reshape);
% f_start_reshape=reshape(f_start_bin,1,[]);
% f_start_reshape(isnan(f_start_reshape))=[];
% f_end_reshape=reshape(f_end_bin,1,[]);
% f_end_reshape(isnan(f_end_reshape))=[];
% 
% begin_id=f_start_reshape==1;
% end_id=f_end_reshape==1;
% start_PFs=(f_COM_reshape(begin_id)-1-(f_end_reshape(begin_id)-f_COM_reshape(begin_id)))<0;
% end_PFs=(50-f_COM_reshape(end_id)-(f_COM_reshape(end_id)-f_start_reshape(end_id)))<0;
% f_remove_id=[];
% if ~isempty(start_PFs)
%     id=find(begin_id);
%     id(start_PFs==0)=[];
%     f_remove_id=[f_remove_id id];
% end
%     
% if ~isempty(end_PFs)
%     id=find(end_id);
%     id(end_PFs==0)=[];
%     f_remove_id=[f_remove_id id];
% end
% 
% new_pf_w_f=PF_w_f;
% new_pf_w_f(f_remove_id)=[];
% f_SP_reshape(f_remove_id)=[];



figure;
% 
% range = 0:1:25;
% histogram(PF_w_f,range,'Normalization','probability')%,'FaceColor',[0.5 0.25 0.7],'EdgeColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.8);
% 
hold on
histogram(PF_width_fam(1,:),range,'FaceColor','#7E2F8E','FaceAlpha',0.4);


% xlabel('PF width');
% title('Fam PF Widths');
%hold on

%%


n_count = 0;
n_start_bin = [];
n_end_bin = [];
nPF_count=0;
n_COM=[];
n_SP=[];
n_all_SP=[];
n_all_remove=[];
for f=1:size(PC_filepaths,2)
    
    if contains(PC_filepaths{f},'_n_')
        n_count = n_count + 1;
        if contains(PC_filepaths{f},['plain']) 
            load(PC_filepaths{f});      
        end
        nPF_count=nPF_count+sum(number_of_PFs(~isnan(number_of_PFs)));
        n_start_bin = [n_start_bin,PF_start_bins];
        n_end_bin = [n_end_bin,PF_end_bins]; 
        COM_M=ones(size(PF_start_bins)).*(-1);
        SP_M=ones(size(PF_start_bins)).*(-1);
        for m=1:size(sig_PFs,1)
           for n=1:size(sig_PFs,2)
              if ~isempty(sig_PFs{m, n})
                 [COM SP]=meanCOMandSP(sig_PFs{m, n}',1,50,50); 
                 COM_M(m,n)=COM;
                 SP_M(m,n)=SP;
                 n_all_SP=[n_all_SP SP];
                 meanPF=mean(sig_PFs{m, n}',1);
                 if meanPF(1)~=0
                    cur_remove=meanPF(1)<=(0.1*max(meanPF));
                 elseif meanPF(end)~=0
                    cur_remove=meanPF(end)<=(0.1*max(meanPF)); 

                 end
                 n_all_remove=[n_all_remove cur_remove];
              end

           end

        end
        n_COM=[n_COM,COM_M];
        n_SP=[n_SP,SP_M];
    end
   
 
    
   
end

n_all_SP(logical(n_all_remove))=[];

for i = 1:length(n_start_bin(:,1))
    for j = 1:length(n_end_bin(1,:))
        PF_width_nov(i,j) = n_end_bin(i,j) - n_start_bin(i,j);
    end
end

% PF_width_nov(PF_width_nov == 0) = nan; 
% PF_w_n=reshape(PF_width_nov,1,[]);
% PF_w_n(isnan(PF_w_n))=[];
% n_SP_reshape=reshape(n_SP,1,[]);
% n_SP_reshape(n_SP_reshape==-1)=[];
% n_COM_reshape=reshape(n_COM,1,[]);
% n_COM_reshape(n_COM_reshape==-1)=[];
% n_COM_reshape=round(n_COM_reshape);
% n_start_reshape=reshape(n_start_bin,1,[]);
% n_start_reshape(isnan(n_start_reshape))=[];
% n_end_reshape=reshape(n_end_bin,1,[]);
% n_end_reshape(isnan(n_end_reshape))=[];
% 
% begin_id=n_start_reshape==1;
% end_id=n_end_reshape==1;
% start_PFs=(n_COM_reshape(begin_id)-1-(n_end_reshape(begin_id)-n_COM_reshape(begin_id)))<0;
% end_PFs=(50-n_COM_reshape(end_id)-(n_COM_reshape(end_id)-n_start_reshape(end_id)))<0;
% n_remove_id=[];
% if ~isempty(start_PFs)
%     id=find(begin_id);
%     id(start_PFs==0)=[];
%     n_remove_id=[n_remove_id id];
% end
%     
% if ~isempty(end_PFs)
%     id=find(end_id);
%     id(end_PFs==0)=[];
%     n_remove_id=[n_remove_id id];
% end
% 
% new_pf_w_n=PF_w_n;
% new_pf_w_n(n_remove_id)=[];
% n_SP_reshape(n_remove_id)=[];


% hold on;
% %range = 0:2.5:25;
% histogram(PF_w_n,range,'Normalization','probability')%,'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor',[0.5 0.25 0.7],'FaceAlpha',0.6);
% 
% % hold on
% %histogram(PF_width_nov(1,:),range,'FaceColor','#7E2F8E','FaceAlpha',0.4);
% 
% 
% xlabel('PF width /bin');
% title('CA1 f n PF Widths');
% legend({'familiar','novel'})
% 
% 
% figure
% cdfplot(PF_w_f);,
% hold on; cdfplot(PF_w_n);
% %xlim([0 30])
% p=ranksum(PF_w_f,PF_w_n);
% title(['CA1 PF width wilcoxon test,p= ' num2str(p)]);
% xlabel('width /bin')
% legend({'familiar','novel'})
% grid off
%%
% range = 0:1:25;
% figure;histogram(new_pf_w_f,range,'Normalization','probability')
% hold on;histogram(new_pf_w_n,range,'Normalization','probability')
% figure
% cdfplot(new_pf_w_f);,
% hold on; cdfplot(new_pf_w_n);
% %xlim([0 30])
% p=ranksum(new_pf_w_f,new_pf_w_n);
% title(['CA1 PF width wilcoxon test,p= ' num2str(p)]);
% xlabel('width /bin')
% legend({'familiar','novel'})
% grid off
%%

save('CA1_SP_fn','n_all_SP','f_all_SP')

%%

%save('CA3_PF_width_n2day', 'PF_w_f','PF_w_n');
%save('male PF_width_nov','PF_width_nov');



            
            
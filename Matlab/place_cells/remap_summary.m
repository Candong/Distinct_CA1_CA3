%clear all; close all;

[remap_filepaths, temp]=uigetfile('*.mat', 'Chose remap info files to load:','MultiSelect','on');
remap_filepaths=sort(remap_filepaths);
%%
f_summary_meantrans=[];
n_summary_meantrans=[];
f_in_n_summary_meantrans=[];
n_fandn_meantrans=cell(14,1);
f_fandn_meantrans=cell(14,1);
fandn_delaylap=cell(14,1);
f_PC=cell(14,1);
n_PC=cell(14,1);
f_neuron=[];
n_neuron=[];

PC_fandn=cell(14,1);
unique_PC_fandn=cell(14,1);
num_of_neurons=0;

f_COM=[];
n_COM=[];
plane_number=3;
for f=1:size(remap_filepaths,2)
    load([temp remap_filepaths{f}]);
   for i=1:2*plane_number
       if i<=plane_number
       PC_id=remap_info(i).place_cellid;
       f_summary_meantrans=[f_summary_meantrans;remap_info(i).meantrans(PC_id,:)];    
       f_COM=[f_COM remap_info(i).COM];   
       f_PC{f}=[f_PC{f} remap_info(i).place_cellid];
       f_neuron=[f_neuron size(remap_info(i).meantrans,1)];   
         
       else
       PC_id=remap_info(i).place_cellid;
       n_summary_meantrans=[n_summary_meantrans;remap_info(i).meantrans(PC_id,:)];         
       n_COM=[n_COM remap_info(i).COM];
       n_PC{f}=[n_PC{f} remap_info(i).place_cellid];
       n_neuron=[n_neuron size(remap_info(i).meantrans,1)];   
       
       f_in_n_PC_id =remap_info(i-plane_number).place_cellid;  
       f_in_n_summary_meantrans=[f_in_n_summary_meantrans;remap_info(i).meantrans(f_in_n_PC_id,:)]; 
       
       unique_PC_fandn{f}=intersect(remap_info(i-plane_number).place_cellid,remap_info(i).place_cellid);
       PC_fandn{f}=remap_info(i).place_cellid(find(ismember(remap_info(i).place_cellid,unique_PC_fandn{f})));
       n_fandn_meantrans{f}=remap_info(i).meantrans(unique_PC_fandn{f},:);
       f_fandn_meantrans{f}= remap_info(i-plane_number).meantrans(unique_PC_fandn{f},:);
       fandn_delaylap{f}=remap_info(i).delay_lap(ismember(PC_id,PC_fandn{f}));
       end
             
   end
   
end

[f_PF_sorted,sort_id]=sort(f_COM);
f_meantrans_sorted=f_summary_meantrans(sort_id,:);
familiar_sort=sort_id;
figure;imagesc(f_meantrans_sorted,[0 5])%,[0 2]);
title('familiar') 
colorbar;

f_in_n_meantrans_sorted=f_in_n_summary_meantrans(familiar_sort,:);
figure;imagesc(f_in_n_meantrans_sorted,[0 5])%,[0 2]);
title('familiar in novel') 
colorbar;

[n_PF_sorted,n_sort_id]=sort(n_COM);
n_meantrans_sorted=n_summary_meantrans(n_sort_id,:);
novel_sort=n_sort_id;
figure;imagesc(n_meantrans_sorted,[0 2.5])%,[0 2]);
title('novel')
colorbar;

%f_meantrans_sorted(isnan(f_meantrans_sorted))=0;
%f_in_n_meantrans_sorted(isnan(f_in_n_meantrans_sorted))=0;
ffn_trans=[];
nfn_trans=[];
inst_corr=[];
total_corr=[];
total_delaylap=[ ];
for f=1:size(remap_filepaths,2)
  cur_corr=[];
  ffn_trans=[ffn_trans; f_fandn_meantrans{f}];
  nfn_trans=[nfn_trans; n_fandn_meantrans{f}];
  n_instcell_id=find(fandn_delaylap{f}==1);
  if ~isempty(f_fandn_meantrans{f}')&~isempty(n_fandn_meantrans{f}')
  cur_corr=diag(corr(f_fandn_meantrans{f}',n_fandn_meantrans{f}'));
  end
  n_instcell=ismember(unique_PC_fandn{f},unique(PC_fandn{f}(n_instcell_id)));
  total_corr=[total_corr; cur_corr];
  inst_corr=[inst_corr; cur_corr(n_instcell)];
  total_delaylap=[total_delaylap fandn_delaylap{f}];
end

%f_in_n_corr=diag(corr(f_fandn_meantrans',n_fandn_meantrans'));
figure;plot(total_corr);
figure;plot(inst_corr);

edges=1:25;
figure;
histogram(total_delaylap,edges)


% instant_cor=f_in_n_corr(fandn_delaylap==1);
% figure;plot(instant_corr);


f_PCnum=[];
n_PCnum=[];
norm_f_PCnum=[];
norm_n_PCnum=[];
f_cell=reshape(f_neuron,[],size(remap_filepaths,2));
for f=1:size(remap_filepaths,2)
   curfnum=length(f_PC{f});
   f_PCnum=[f_PCnum curfnum];
   norm_f_PCnum=[norm_f_PCnum curfnum/sum(f_cell(:,f))];
   curnnum=length(n_PC{f});
   n_PCnum=[n_PCnum curnnum];
   norm_n_PCnum=[norm_n_PCnum curnnum/sum(f_cell(:,f))];
end

figure;
bar([f_PCnum mean(f_PCnum);n_PCnum mean(n_PCnum)]')

figure;
bar([norm_f_PCnum mean(norm_f_PCnum);norm_n_PCnum mean(norm_n_PCnum)]')

[h,p,ci,stats]= ttest(norm_f_PCnum,norm_n_PCnum)

figure;
plot([1 2],[norm_f_PCnum; norm_n_PCnum;],'-o')
ylim([0 1])


figure;
plot([1], total_corr','o')
hold on; plot([1],mean(total_corr(~isnan(total_corr))),'*');

figure;
plot([1], inst_corr','o')
hold on; plot([1],mean(inst_corr(~isnan(inst_corr))),'*');
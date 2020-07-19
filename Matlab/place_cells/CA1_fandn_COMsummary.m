clear all; close all;

[COM_filepaths, temp]=uigetfile('*.mat', 'Chose same threshold files to load:','MultiSelect','on');
plane_name={'plain1','plain2','plain3'};
plane_num=2;

f_COM=cell(1,plane_num);
f_SP=cell(1,plane_num);
n_COM=cell(1,3);
n_SP=cell(1,3);
f_empty=cell(1,3);
n_empty=cell(1,3);

for f=1:size(COM_filepaths,2)
 load([temp COM_filepaths{f}]);  
 
 for i=1:size(PF_info,2)
     for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_f_')
      
     if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
%           f_combine_pfid{1}=PF_info(i).combinePC;
          f_late15_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
%           f_combine_pfid{2}=PF_info(i).combinePC;
          f_late15_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
%           f_combine_pfid{3}=PF_info(i).combinePC;
          f_late15_id{3}=i;      
          end
     elseif  contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          f_1st25_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          f_1st25_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          f_1st25_id{3}=i;      
          end
     end
  end
  
    if contains(PF_info(i).filepaths,'_n_')
       if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
%           f_combine_pfid{1}=PF_info(i).combinePC;
          n_late15_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
%           f_combine_pfid{2}=PF_info(i).combinePC;
          n_late15_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
%           f_combine_pfid{3}=PF_info(i).combinePC;
          n_late15_id{3}=i;      
          end
     elseif  contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_1st25_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_1st25_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          n_1st25_id{3}=i;      
          end
     end
  end
    
end
end
    

 
 for p=1:3
     f_COM{p}=[f_COM{p} PF_info(f_1st25_id{p}).COM];
     f_SP{p}=[f_SP{p} PF_info(f_1st25_id{p}).SP];
     
     n_COM{p}=[n_COM{p} PF_info(n_1st25_id{p}).COM];
     n_SP{p}=[n_SP{p} PF_info(n_1st25_id{p}).SP];
     
     f_delay_lap=PF_info(f_1st25_id{p}).delay_lap;
     f_PF_empty=PF_info(f_1st25_id{p}).PF_empty;
     f_empty{p}=[f_empty{p}  num_emp(f_delay_lap, f_PF_empty)];
     
     n_delay_lap=PF_info(n_1st25_id{p}).delay_lap;
     n_PF_empty=PF_info(n_1st25_id{p}).PF_empty;
     n_empty{p}=[n_empty{p} num_emp(n_delay_lap, n_PF_empty)];
     
%      if length(PF_info(f_1st25_id{p}).SP)~=size(PF_info(f_1st25_id{p}).PF_empty,2)
%      p
%      PF_info(f_1st25_id{p}).filepaths
%      end
 end
 
 

end
f_COM_mean=cellfun(@mean,f_COM);
f_SP_mean=cellfun(@mean,f_SP);

n_COM_mean=cellfun(@mean,n_COM);
n_SP_mean=cellfun(@mean,n_SP);

SP_f_n=[f_SP_mean;n_SP_mean];
figure;
bar(SP_f_n');

f_stderror=cellfun(@std,f_SP)./sqrt(cellfun(@length,f_SP));
n_stderror=cellfun(@std,n_SP)./sqrt(cellfun(@length,f_SP));
hold on
er = errorbar([0.8 1.8 2.8; 1.2 2.2 3.2]',SP_f_n',[f_stderror;n_stderror]');
er(1).LineStyle = 'none';  

for p=1:3
 [f_sort_emp{p} f_I]=sort(f_empty{p});   
 f_sort_SP{p}= f_SP{p}(f_I);  
  
 [n_sort_emp{p} n_I]=sort(n_empty{p});   
 n_sort_SP{p}= n_SP{p}(n_I); 
    
end 

% figure
% plot(f_sort_emp{1}, f_sort_SP{1},'o')
% mdl=fitlm(f_sort_emp{1},f_sort_SP{1});
% r_square=mdl.Rsquared.Ordinary;
% 
% P=polyfit(f_sort_emp{1},f_sort_SP{1},1);
% yfit=polyval(p,f_sort_emp{1});
% yresid=f_sort_SP{1}-yfit;
% SSresid=sum(yresid.^2);
% SStotal=(length(f_sort_SP{1})-1)*var(f_sort_SP{1});
% R2=1-SSresid/SStotal;
% 
% figure
% plot(n_sort_emp{1}, n_sort_SP{1},'o')


%% 
function emptylapnum=num_emp(delay_lap, PF_empty) 
emptylapnum=[];
for lap=1:length(delay_lap)
        cur_empty=PF_empty{lap};
        cur_empty(cur_empty<delay_lap(lap))=[];
        emptylapnum=[emptylapnum length(cur_empty)]; 
         
         
     end

end
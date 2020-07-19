% this script calculate the COM of each PF
clear all; close all;
%COM familiar

% delay lap of familiar PF found by last 15 laps
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose info files to load:','MultiSelect','on');

field_name={'plane1','plane2','plane3'};

PF_mean=[];
PF_count=1;
PF_start_lap=[];

load([temp behavior_filepaths]);


%% COM in f
for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_f_')
      
     if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          f_combine_pfid{1}=PF_info(i).combinePC;
          f_late15_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          f_combine_pfid{2}=PF_info(i).combinePC;
          f_late15_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          f_combine_pfid{3}=PF_info(i).combinePC;
          f_late15_id{3}=i;      
          end
     elseif  contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          f_combine_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          f_combine_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          f_combine_id{3}=i;      
      end
     end
  end
    
end

%% find emptylapjk COM for f
f_PF_empty=cell(1,3);
for plane=1:3
 
 if ~isempty(f_combine_pfid{plane})

 start_bin=PF_info(f_late15_id{plane}).combinePC_start;
 end_bin=PF_info(f_late15_id{plane}).combinePC_end;
 
 end
numbins=50;


%%
count=1;
pf_COM=[];
pf_SP=[];
PF_empty={};

for ii = 1:length(f_combine_pfid{plane})
        
binMean=PF_info(f_combine_id{plane}).combinebinMean(:,:,ii);
%     figure;
%    imagesc(binMean);
%    title([num2str(plane) ' ' num2str(f_combine_pfid{plane}(ii))])
%%%%%%%%%%%%%%%%%%%%%%
PF_denoise=binMean;


start_b=start_bin(ii);
end_b=end_bin(ii);
[meanCOM cur_SP]=meanCOMandSP(PF_denoise,start_b,end_b,numbins);
pf_COM=[pf_COM meanCOM ]; 
pf_SP=[pf_SP cur_SP ];
%find(sum(PF_withoutnoise(:,start_b:end_b),2)==0)
PF_empty{ii}=find(sum(PF_denoise(:,start_b:end_b),2)==0);

f_PF_empty{plane}=PF_empty;
%f_sig_PF{plane}=PF_info(f_combine_id{plane}).combinebinMean;
end

f_COM{plane}=pf_COM;
f_SP{plane}=pf_SP;
PF_info(f_combine_id{plane}).COM=f_COM{plane};
PF_info(f_combine_id{plane}).SP=f_SP{plane};
PF_info(f_combine_id{plane}).PF_empty=f_PF_empty{plane};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_n_')
      
     if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_combine_pfid{1}=PF_info(i).combinePC;
          n_late15_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_combine_pfid{2}=PF_info(i).combinePC;
          n_late15_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          n_combine_pfid{3}=PF_info(i).combinePC;
          n_late15_id{3}=i;      
          end
     elseif  contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_combine_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_combine_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          n_combine_id{3}=i;      
      end
     end
  end
    
end

%% find delay COM for n
n_PF_empty=cell(1,3);
for plane=1:3
 
 if ~isempty(n_combine_pfid{plane})

 start_bin=PF_info(n_late15_id{plane}).combinePC_start;
 end_bin=PF_info(n_late15_id{plane}).combinePC_end;
 
 end
numbins=50;


count=1;
pf_COM=[];
pf_SP=[];
PF_empty={};
n_PF_empty=cell(1,3);
for ii = 1:length(n_combine_pfid{plane})
        
binMean=PF_info(n_combine_id{plane}).combinebinMean(:,:,ii);
%     figure;
%    imagesc(binMean);
%    title([num2str(plane) ' ' num2str(n_combine_pfid{plane}(ii))])
%%%%%%%%%%%%%%%%%%%%%%
PF_denoise=binMean;


start_b=start_bin(ii);
end_b=end_bin(ii);
[meanCOM cur_SP]=meanCOMandSP(PF_denoise,start_b,end_b,numbins);
pf_COM=[pf_COM meanCOM ]; 
pf_SP=[pf_SP cur_SP ];
%find(sum(PF_withoutnoise(:,start_b:end_b),2)==0)
PF_empty{ii}=find(sum(PF_denoise(:,start_b:end_b),2)==0);

n_PF_empty{plane}=PF_empty;
%n_sig_PF{plane}=PF_info(n_combine_id{plane}).combinebinMean;
end
clear PF_empty
clear sig_PF
n_COM{plane}=pf_COM;
n_SP{plane}=pf_SP;
PF_info(n_combine_id{plane}).COM=n_COM{plane};
PF_info(n_combine_id{plane}).SP=n_SP{plane};
PF_info(n_combine_id{plane}).PF_empty=n_PF_empty{plane};

end



 save([PF_info(1).filepaths(1:17) '_COM_and_SP'],'PF_info');%'f_SP','f_COM','n_SP','f_COM','f_combine_pfid','n_combine_pfid','f_PF_empty','n_PF_empty','f_sig_PF','n_sig_PF');

 



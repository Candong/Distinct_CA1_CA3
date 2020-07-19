clear all; close all;


%% 
[transient_path_day1, temp]=uigetfile('*.mat', 'Chose day1 transient files to load:','MultiSelect','on');

[PF_f_path_day1, temp]=uigetfile('*.mat', 'Chose day1 f PF files to load:','MultiSelect','on');
[PF_n_path_day1, temp]=uigetfile('*.mat', 'Chose day1 n PF files to load:','MultiSelect','on');


if ~isa(transient_path_day1,'cell') 
  transient_path_day1={transient_path_day1};  
  PF_f_path_day1={PF_f_path_day1};  
  PF_n_path_day1={PF_n_path_day1}; 
end 

PF_n_path_day1=sort(PF_n_path_day1);
PF_f_path_day1=sort(PF_f_path_day1);
transient_path_day1=sort(transient_path_day1);


%%
sample_rate=10.4;
total_peak=[];
total_width=[];
total_frq=[];

total_PF_id=[];
PF_total_peak=[];
PF_total_width=[];
PF_total_frq=[];

for f=1:size(transient_path_day1,2)
    
    if transient_path_day1{f}~=0
    load([temp transient_path_day1{f}]);
    load([temp PF_f_path_day1{f}]);
    f_PF_id=~isnan(number_of_PFs);
    load([temp PF_n_path_day1{f}]);
    n_PF_id=~isnan(number_of_PFs);
    
    PF_id=f_PF_id|n_PF_id;
    
    Fc3=data.Fc3;
    %PF_Fc3=Fc3(PF_id);
    
    time=size(Fc3,1)/sample_rate;
    cur_width=[];
    cur_peak=[];
    %PF_width=[];
    %PF_peak=[];
    for i=1:size(Fc3,2)
        E=max(bwlabel(Fc3(:,i)));
        trans_label=bwlabel(Fc3(:,i));

        for j=1:E
            peak=max(Fc3(trans_label==j));
            width=length(Fc3(trans_label==j))/sample_rate;
            cur_width=[cur_width width];
            cur_peak=[cur_peak peak];

        end
        if PF_id(i)==1
            PF_total_width=[PF_total_width cur_width];
            PF_total_peak=[PF_total_peak cur_peak];
            PF_total_frq=[PF_total_frq E/time];
             
        end
        
            
            
            

        total_width=[total_width cur_width];
        total_peak=[total_peak cur_peak];
        total_frq=[total_frq E/time];
     
    end
    sum(PF_id)
    

    % figure;hold on;
    % plot(f_v.bylap);
    % plot(n_v.bylap)
    end
end
%behavior_filepaths=sort(behavior);
%% seperate laps by reward

[transient_path_day2, temp]=uigetfile('*.mat', 'Chose Chose day2 transient files files to load:','MultiSelect','on');
[PF_f_path_day2, temp]=uigetfile('*.mat', 'Chose day1 f PF files to load:','MultiSelect','on');
[PF_n_path_day2, temp]=uigetfile('*.mat', 'Chose day1 n PF files to load:','MultiSelect','on');

if ~isa(transient_path_day2,'cell') 
  transient_path_day2={transient_path_day2};  
  PF_f_path_day2={PF_f_path_day2};  
  PF_n_path_day2={PF_n_path_day2};   
end 

if transient_path_day2{1}~=0
PF_n_path_day2=sort(PF_n_path_day2);
PF_f_path_day2=sort(PF_f_path_day2);
transient_path_day2=sort(transient_path_day2);
end
% total_peak=[];
% total_width=[];
% total_frq=[];
for f=1:size(transient_path_day2,2)
    
    if transient_path_day2{f}~=0
    load([temp transient_path_day2{f}]);
    load([temp PF_f_path_day2{f}]);
    f_PF_id=~isnan(number_of_PFs);
    load([temp PF_n_path_day2{f}]);
    n_PF_id=~isnan(number_of_PFs);
    
    PF_id=f_PF_id|n_PF_id;    
    
    
    
    Fc3=data.Fc3;

    time=size(Fc3,1)/sample_rate;
    cur_width=[];
    cur_peak=[];
    for i=1:size(Fc3,2)
        E=max(bwlabel(Fc3(:,i)));
        trans_label=bwlabel(Fc3(:,i));

        for j=1:E
            peak=max(Fc3(trans_label==j));
            width=length(Fc3(trans_label==j))/sample_rate;
            cur_width=[cur_width width];
            cur_peak=[cur_peak peak];
 
        end
        
        if PF_id(i)==1
            PF_total_width=[PF_total_width cur_width];
            PF_total_peak=[PF_total_peak cur_peak];
            PF_total_frq=[PF_total_frq E/time];
             
        end
        
        
        total_width=[total_width cur_width];
        total_peak=[total_peak cur_peak];
        total_frq=[total_frq E/time];
        
        
        
        
        
    end

    % figure;hold on;
    % plot(f_v.bylap);
    % plot(n_v.bylap)
    end
end
%%
save([transient_path_day1{1}(1:10) 'transient_info_withPF'],'total_frq','total_width','total_peak','PF_total_frq','PF_total_width','PF_total_peak');

%%


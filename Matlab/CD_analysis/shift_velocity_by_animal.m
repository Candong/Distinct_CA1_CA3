clear all; %close all;
[behavior, temp]=uigetfile('*.mat', 'Chose velocity files to load:','MultiSelect','on');

load(behavior);
day1=~isempty(fieldnames(f1_v));
day2=~isempty(fieldnames(f2_v));
nday2=~isempty(fieldnames(nday2_v));

window=1;
onsetlap=0;
day1_f_COM=[];
day1_n_COM=[];


if day1
    [day1_file, temp]=uigetfile('*.mat', 'Chose day1 pf files f and n to load:','MultiSelect','on');
    day1_file=sort(day1_file);
    for f=1:size(day1_file,2)
        if contains(day1_file{f},'_f_')
            load([temp day1_file{f}]);
            day1_f_sig_PFs=sig_PFs;
            cur_day1_f_COM=calculate_COM_by_window(day1_f_sig_PFs,window,onsetlap);
            day1_f_COM=[day1_f_COM; cur_day1_f_COM];
        elseif contains(day1_file{f},'_n_')
            load([temp day1_file{f}]);
            day1_n_sig_PFs=sig_PFs;
            cur_day1_n_COM=calculate_COM_by_window(day1_n_sig_PFs,window,onsetlap);
            day1_n_COM=[day1_n_COM; cur_day1_n_COM];
        end
        
        
        
    end
    

end
day1_f_COM_diff=diff(day1_f_COM,1,2);
day1_n_COM_diff=diff(day1_n_COM,1,2);
% figure; plot(f1_v.bylap(1:size(day1_f_COM_diff,2)),day1_f_COM_diff,'o');  
% title('familiar')
% figure; plot(n1_v.bylap(1:size(day1_n_COM_diff,2)),day1_n_COM_diff,'o'); 

day2_f_COM=[];
day2_n_COM=[];

if day2
    [day2_file, temp]=uigetfile('*.mat', 'Chose day2 pf files f and n to load:','MultiSelect','on');
    day2_file=sort(day2_file);
    for f=1:size(day2_file,2)
        if contains(day2_file{f},'_f_')
            load([temp day2_file{f}]);
            day2_f_sig_PFs=sig_PFs;
            cur_day2_f_COM=calculate_COM_by_window(day2_f_sig_PFs,window,onsetlap);
            day2_f_COM=[day2_f_COM; cur_day2_f_COM];
        elseif contains(day2_file{f},'_n_')
            load([temp day2_file{f}]);
            day2_n_sig_PFs=sig_PFs;
            cur_day2_n_COM=calculate_COM_by_window(day2_n_sig_PFs,window,onsetlap);
            day2_n_COM=[day2_n_COM; cur_day2_n_COM];
        end
      
    end    

end
day2_f_COM_diff=diff(day2_f_COM,1,2);
day2_n_COM_diff=diff(day2_n_COM,1,2);
%     figure; plot(f2_v.bylap(1:size(day2_f_COM_diff,2)),day2_f_COM_diff,'o');  
% title('familiar')
% figure; plot(n2_v.bylap(1:size(day2_n_COM_diff,2)),day2_n_COM_diff,'o'); 


nday2_f_COM=[];
nday2_n_COM=[];

if nday2
    [nday2_file, temp]=uigetfile('*.mat', 'Chose nday2 pf files f and n to load:','MultiSelect','on');
    nday2_file=sort(nday2_file);
    for f=1:size(nday2_file,2)
        if contains(nday2_file{f},'_n_')
            load([temp nday2_file{f}]);
            nday2_n_sig_PFs=sig_PFs;
            cur_nday2_n_COM=calculate_COM_by_window(nday2_n_sig_PFs,window,onsetlap);
            nday2_n_COM=[nday2_n_COM; cur_nday2_n_COM];
        end
      
    end    

end
nday2_n_COM_diff=diff(nday2_n_COM,1,2);
%     figure; plot(nday2_v.bylap(1:size(nday2_n_COM_diff,2)),nday2_n_COM_diff,'o'); 

    
save([behavior(1:end-17) '_shift_and_velocity_info'],'f1_v','n1_v','f2_v','n2_v','nday2_v',...
    'nday2_n_COM_diff','day2_f_COM_diff','day2_n_COM_diff','day1_f_COM_diff','day1_n_COM_diff');
    
    
%%

function COM=calculate_COM_by_window(sig_PFs,window,onsetlap)
COM=[];
    for i =1:size(sig_PFs,2)
       for j=1:size(sig_PFs,1)
           if ~isempty(sig_PFs{j,i})
               cur_binmean=sig_PFs{j,i};
               start_bin=1;
               end_bin=size(cur_binmean,1);
               start_lap=find_delaylap(cur_binmean,start_bin,end_bin,6,3,size(cur_binmean,2));
            if ~isclipped(cur_binmean) & start_lap>1
            if onsetlap==1
                cur_binmean=cur_binmean(:,startlap:end);
            end
            
            lapnum=size(cur_binmean,2);  
            %lapnum=min(lapnum,stoplap+stepsize);
            
            cur_COMv=NaN(1,lapnum);

              
            for n=1:(lapnum-window+1)
                step_binmean=cur_binmean(:,n:n+window-1);
                [meanCOM SP]=meanCOMandSP(step_binmean',1,50,50);
                cur_COMv(n)=meanCOM;
            end
            COM=[COM; cur_COMv];
            end
           end
        end
     end
end





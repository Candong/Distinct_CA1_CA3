function [COM all_start_lap pf_id]=calculate_and_linreg_eachCOM_by_window(sig_PFs,window,onsetlap)
% onsetlap =0: calculate each lap 
%onsetlap0
COM=[];
total_SP=[];
all_start_lap=[];
COM_start=[];
COM_alllaps=[];
COM_end=[];
pf_id=[];
    for i =1:size(sig_PFs,2)
       %for j=1:size(sig_PFs,1)
       for j=1:size(sig_PFs,1)
           if ~isempty(sig_PFs{j,i})
               cur_binmean=sig_PFs{j,i};
               start_bin=1;
               end_bin=size(cur_binmean,1);
               %figure;imagesc(cur_binmean');
               start_lap=find_delaylap(cur_binmean,start_bin,end_bin,6,3,size(cur_binmean,2));
            if ~isclipped(cur_binmean)% & start_lap>1
%             if onsetlap==1
%                 cur_binmean=cur_binmean(:,startlap:end);
%             end
            
            lapnum=size(cur_binmean,2);  
            %lapnum=min(lapnum,stoplap+stepsize);
            
            cur_COMv=NaN(1,lapnum);

            [com sp]=meanCOMandSP(cur_binmean',1,50,50);  
            total_SP=[total_SP sp];
            COM_alllaps=[COM_alllaps com];
            
            for n=1:(lapnum-window+1)
                step_binmean=cur_binmean(:,n:n+window-1);
                [meanCOM SP]=meanCOMandSP(step_binmean',1,50,50);
                cur_COMv(n)=meanCOM;
            end
            COM=[COM; cur_COMv];
            all_start_lap=[all_start_lap start_lap];
            pf_id=[pf_id i];
            end
           end
        end
     end
end
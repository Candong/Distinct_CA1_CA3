function [COM_start COM_end COM_alllaps]=calculate_start_end_alllaps_COM(sig_PFs,step,onsetlap)
% onsetlap =0: calculate each lap 
%onsetlap0

COM_start=[];
COM_alllaps=[];
COM_end=[];
    for i =1:size(sig_PFs,2)
       for j=1:size(sig_PFs,1)
           if ~isempty(sig_PFs{j,i})
               cur_binmean=sig_PFs{j,i};
               start_bin=1;
               end_bin=size(cur_binmean,1);
               start_lap=find_delaylap(cur_binmean,start_bin,end_bin,6,3,size(cur_binmean,2));
            if ~isclipped(cur_binmean) %& start_lap>1
            if onsetlap==1
                cur_binmean=cur_binmean(:,start_lap:end);
            end

            [com sp]=meanCOMandSP(cur_binmean',1,50,50);  
            COM_alllaps=[COM_alllaps com];
            %total_SP=[total_SP sp];
            [com_start sp_start]=meanCOMandSP(cur_binmean(:,1:step)',1,50,50);  
            COM_start=[COM_start com_start];
            [com_end sp_end]=meanCOMandSP(cur_binmean(:,end-step+1:end)',1,50,50);  
            COM_end=[COM_end com_end];            
            

            end
           end
        end
     end
end
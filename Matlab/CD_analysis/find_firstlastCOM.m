function [meanCOM_first5 meanCOM_last5 Act_lapnum norm_act_COMshift norm_from_startCOMshift total_lapstart]=find_firstlastCOM(f_sig_PFs,n_sig_PFs,common_PC_id,PF_start_bins,PF_end_bins,numbins,all_bins,ave_bin)
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
                n_cur_binmean=n_sig_PFs{p,i}{1,common_PC_id{p,i}(j)}';
                if ~isclipped(cur_binmean') & ~isclipped(n_cur_binmean')
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
                end
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
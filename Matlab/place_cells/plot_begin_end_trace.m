clear all; close all;
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose f/n files to load:','MultiSelect','on');
load([temp behavior_filepaths])
step=5;
for i=1:200%1:size(sig_PFs,2)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs{j,i};
            binmean_temp=binmean_temp(:,1:30);
            binmean_temp(:,sum(binmean_temp,1)==0)=[];
            if size(binmean_temp,2)>=5
            start_skewness=find_pfskewness(binmean_temp(:,1:step),step);
            end_skewness=find_pfskewness(binmean_temp(:,end-step+1:end),step);
            if size(binmean_temp,2)>=step+1 & ~isclipped(binmean_temp)
                if (start_skewness-end_skewness)>0
                start_mean=mean(binmean_temp(:,1:step),2);
                end_mean=mean(binmean_temp(:,end-step+1:end),2);
                remove_start=find(start_mean~=0);
                remove_end=find(end_mean~=0);
                start_mean=start_mean(remove_start(1):remove_start(end));
                end_mean = end_mean(remove_end(1):remove_end(end));
                
                begin_trace=smooth(mean(binmean_temp(:,1:step),2));
                end_trace=smooth(mean(binmean_temp(:,end-step+1:end),2));
                
                figure; hold on;
                plot(begin_trace/max(begin_trace),'Linewidth',1);
                %line([mean(start_mean) mean(start_mean)],[0 1])
                plot(end_trace/max(end_trace),'Linewidth',1);
                %line([mean(end_mean) mean(end_mean)],[0 1])
                title(['neuron id=' num2str(i) ' ' num2str(start_skewness)  '  '  num2str(end_skewness)])
                legend({'begin','end'})
                
%                 figure; hold on;
%                 plot(begin_trace,'Linewidth',1);
%                 line([0 50],[mean(start_mean) mean(start_mean)],'Color','black')
%                 plot(end_trace,'Linewidth',1);
%                 line([0 50],[mean(end_mean) mean(end_mean)],'Color','red')
%                 title(['neuron id=' num2str(i) ' ' num2str(start_skewness)  '  '  num2str(end_skewness)])
%                 legend({'begin','begin_mean','end','end_mean'})                
                end
            end
            end
%             f=figure;imagesc((binmean_temp'./max(binmean_temp)'));
%             colorbar;
%             %movegui(f,'southwest');
%             title(['neuron id=' num2str(i)])
%             name=fullfile([temp 'pic\'],num2str(i));
%             %saveas(f,name,'epsc')
%             %pause(1)
%             colormap jet
%             colorbar;
    
        end
    
    end
                    
                    
                    
                    
                    
end
clear all; %close all;
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose f/n files to load:','MultiSelect','on');
load([temp behavior_filepaths])
figure;

subplot_x=10;
subplot_y=10;%ceil(sum(number_of_PFs(~isnan(number_of_PFs)))/subplot_x);
PC_count=0;
PC_id=find(~isnan(number_of_PFs));
for i=PC_id(1:100)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs{j,i};
            PC_count=PC_count+1;
%             if j==1
%                 PC_count=PC_count+1;
%             elseif j>1
%                 PC_count=PC_count+1;
%             end
%             f=figure;imagesc(binmean_temp');
%             colorbar;
%             %movegui(f,'southwest');
%             title(['neuron id=' num2str(i)])
%             name=fullfile([temp 'pic\'],num2str(i));
            %saveas(f,name,'epsc')
            %pause(1)
    
        end
subplot(subplot_y,subplot_x,PC_count);
imagesc(binmean_temp');
title(['neuron id=' num2str(i)])    
    end
                    

                    
                    
                    
end
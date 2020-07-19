clear all; close all;
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose f/n files to load:','MultiSelect','on');
load([temp behavior_filepaths])

for i=1:size(sig_PFs,2)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs{j,i};
            f=figure;imagesc((binmean_temp'./max(binmean_temp)'));
            colorbar;
            %movegui(f,'southwest');
            title(['neuron id=' num2str(i)])
            name=fullfile([temp 'pic\'],num2str(i));
            %saveas(f,name,'epsc')
            %pause(1)
            colormap jet
            colorbar;
    
        end
    
    end
                    
                    
                    
                    
                    
end

%%
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose f/n files to load:','MultiSelect','on');
load([temp behavior_filepaths])

for i=1:size(sig_PFs,2)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs{j,i};
            f=figure;imagesc((binmean_temp'./max(binmean_temp)'));
            colorbar;
            %movegui(f,'southwest');
            title(['neuron id=' num2str(i)])
            name=fullfile([temp 'pic\'],num2str(i));
            %saveas(f,name,'epsc')
            %pause(1)
            colormap jet
            colorbar;
    
        end
    
    end
                    
                    
                    
                    
                    
end
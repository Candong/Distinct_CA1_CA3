clear all; close all;

[pf_file, temp]=uigetfile('*.mat', 'Chose PF files to load:','MultiSelect','on');

load(pf_file);

window=2;
onsetlap=1;
%function [slopeall_animal_shift_para(sig_PFs,window,onsetlap)
[COM all_start_lap]=calculate_and_linreg_eachCOM_by_window(sig_PFs,window,onsetlap);

delta_COM=[];
Rsquare=[];
slope=[];
Rsquare2=[];
slope2=[];
var_y=[];
var_y2=[];
total_var_var=[];
%figure; hold on;
    for i=1:size(COM,1)
        cur_deltaCOM= COM(i,:)-COM(i,all_start_lap(i)); 
        delta_COM=[delta_COM; cur_deltaCOM ];
        % figure;hold on;
        % plot(1:size(COM,2), cur_deltaCOM,'o');
        % ylim([-10 10])

        x=1:size(COM,2);
        x(isnan(cur_deltaCOM))=[];
        y=cur_deltaCOM(~isnan(cur_deltaCOM));
        if onsetlap==0 %start regression from first lap with activity
            cur_var_y=var(y);
            var_y=[var_y cur_var_y];
            model=fitlm(x, y);
            % plot(x,model.('Coefficients').('Estimate')(2).*x+ model.('Coefficients').('Estimate')(1));
            % xlim([1 size(COM,2)]);

            slope=[slope model.('Coefficients').('Estimate')(2)];
            cur_r_squared=model.Rsquared.Ordinary;
            Rsquare=[Rsquare cur_r_squared];

        else
            x2=x(x>=all_start_lap(i));
            y2=y(x>=all_start_lap(i));
            cur_var_y2=var(y2);
            var_y2=[var_y2 cur_var_y2];
            model2=fitlm(x2,y2);
            slope2=[slope2 model2.('Coefficients').('Estimate')(2)];
            cur_r_squared2=model2.Rsquared.Ordinary;
            Rsquare2=[Rsquare2 model2.Rsquared.Ordinary];
            var_var=windowed_var(y2,3);
            total_var_var=[total_var_var var_var];

        end




    % x2=x(x>=all_start_lap(i));
    % y2=y(x>=all_start_lap(i));
    % cur_var_y2=var(y2);
    % var_y2=[var_y2 cur_var_y2];
    % model2=fitlm(x2,y2);
    % slope2=[slope2 model2.('Coefficients').('Estimate')(2)];
    % cur_r_squared2=model2.Rsquared.Ordinary;
    % Rsquare2=[Rsquare2 model2.Rsquared.Ordinary];
    % mean_var=windowed_var(y2,3);
    % total_mean_var=[total_mean_var mean_var];
    % 
    if var_var>=95%1.5 & cur_r_squared2<=0.3 % figure;hold on;
    figure;hold on;
    plot(1:size(COM,2), cur_deltaCOM,'o');
    ylim([-10 10])
    plot(x2,model2.('Coefficients').('Estimate')(2).*x2+ model2.('Coefficients').('Estimate')(1));
    xlim([1 size(COM,2)]);
    title([num2str(cur_r_squared2) ' ' num2str(var_var)]);
    end

    end
%end

figure;histogram(total_var_var,'Normalization','probability')

% figure;plot(slope2,Rsquare2,'o')
% %re;plot(slope2,Rsquare2,'o')
% figure;plot3(slope2, Rsquare2, var_y2,'o')
% xlabel('slope');ylabel('rsquare');zlabel('var')
% figure;plot3(slope2, Rsquare2, total_mean_var,'o')
% xlabel('slope');ylabel('rsquare');zlabel('mean_var')

%%
% function mean_var=windowed_var(y,window)
% total_var=[];
%     for n=1:(length(y)-window+1)
%         cur_var=var(y(:,n:n+window-1));
%         total_var=[total_var cur_var];
%     end
% mean_var=mean(total_var);
% 
% 
% end
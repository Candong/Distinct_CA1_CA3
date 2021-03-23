function [slope pvalue all_start_lap COM_start COM_end COM_alllaps onset_deltaCOM all_deltaCOM pf_id Rsquare]=caculate_shift_parameters(sig_PFs,window,onsetlap,step,fit_lapnum)
% if fit_lapnum is empty, use all laps to regress

delta_COM=[];
Rsquare=[];
pvalue=[];
slope=[];
%Rsquare2=[];
slope=[];
var_y=[];
var_y2=[];
onset_deltaCOM=[];
all_deltaCOM=[];
if ~isempty(fit_lapnum)
fit_lapnum_start=fit_lapnum(1);
fit_lapnum_end=fit_lapnum(2);
end


[COM all_start_lap pf_id]=calculate_and_linreg_eachCOM_by_window(sig_PFs,window,onsetlap);
[COM_start COM_end COM_alllaps]=calculate_start_end_alllaps_COM(sig_PFs,step,onsetlap);

% total_mean_var=[];
%figure; hold on;
    for i=1:size(COM,1)
        cur_deltaCOM= COM(i,:)-COM(i,all_start_lap(i)); 
        delta_COM=[delta_COM; cur_deltaCOM ];
        onset_deltaCOM=[onset_deltaCOM cur_deltaCOM(all_start_lap(i)+1)];
        all_deltaCOM=[all_deltaCOM; diff(COM(i,:))];
%         figure;hold on;
%         plot(1:size(COM,2), cur_deltaCOM,'o');
%         ylim([-10 10])

        x=1:size(COM,2);
        x(isnan(cur_deltaCOM))=[];
        y=cur_deltaCOM(~isnan(cur_deltaCOM));
        if onsetlap==0 %start regression from first lap with activity
%             cur_var_y=var(y);
%             var_y=[var_y cur_var_y];
            model=fitlm(x, y);
            plot(x,model.('Coefficients').('Estimate')(2).*x+ model.('Coefficients').('Estimate')(1));
            xlim([1 size(COM,2)]);
            model_summary=anova(model,'summary');
            cur_p=model_summary.('pValue')(2);
            pvalue=[pvalue cur_p];
            slope=[slope model.('Coefficients').('Estimate')(2)];
            cur_r_squared=model.Rsquared.Ordinary;
            Rsquare=[Rsquare cur_r_squared];

        else
            
            x2=x(x>=all_start_lap(i));
            y2=y(x>=all_start_lap(i));
            if ~isempty(fit_lapnum)
                select_id=(x2<=all_start_lap(i)+fit_lapnum_end-1)& (x2>=all_start_lap(i)+fit_lapnum_start-1);
                y2=y2(select_id);
                x2=x2(select_id);
                
            end
%             cur_var_y2=var(y2);
%             var_y2=[var_y2 cur_var_y2];
            if length(x2)>=3
            model=fitlm(x2,y2);
            slope=[slope model.('Coefficients').('Estimate')(2)];
            model_summary=anova(model,'summary');
            cur_p=model_summary.('pValue')(2);
            pvalue=[pvalue cur_p];
            cur_r_squared=model.Rsquared.Ordinary;
            Rsquare=[Rsquare cur_r_squared];
%             cur_r_squared2=model.Rsquared.Ordinary;
%             Rsquare2=[Rsquare2 model.Rsquared.Ordinary];
%             mean_var=windowed_var(y2,3);
%             total_mean_var=[total_mean_var mean_var];
%             if model.('Coefficients').('Estimate')(2)<-0.05 & cur_r_squared>0.5 & length(x2)>20 & cur_p<0.05 & model.('Coefficients').('Estimate')(2)> -0.13
%             figure;hold on;
%             plot(1:size(COM,2), cur_deltaCOM,'o');
%             ylim([-10 10])
%             plot(x2,model.('Coefficients').('Estimate')(2).*x2+ model.('Coefficients').('Estimate')(1));
%             xlim([1 size(COM,2)]);
%             title(['slope' num2str(model.('Coefficients').('Estimate')(2)),'R-s  ' num2str(model.Rsquared.Ordinary)])
%             text(25,-6, ['p-value=' num2str(cur_p)])
%             end
            else
                slope=[slope -100];
            end
        end





    end
end

function [f_width f_skewness f_COM f_COM_slope slope_sig f_skew_slope skew_sig c_r s_r]=calculate_by_lap_skewness_COM(f_binmean,onsetlap,stepsize_com,stepsize_skew,plot_on)
window = 6;
threshold = 3;
          
    if onsetlap~=0
        cur_startlap_f = find_delaylap(f_binmean,1,size(f_binmean,1),window,threshold,size(f_binmean,2));
        f_binmean = f_binmean(:,cur_startlap_f:end);

    end
    
    lapnum = size(f_binmean,2);
    f_skewness = find_pfskewness(f_binmean,stepsize_skew);
    f_COM = find_pfCOM (f_binmean, stepsize_com);
    f_width=find_pfwidth(f_binmean,stepsize_skew);
    
    x=1:length(f_COM);
    y=f_COM;
    x(isnan(f_COM))=[];
    y(isnan(f_COM))=[];
    model=fitlm(x, y);
% plot(x,model.('Coefficients').('Estimate')(2).*x+ model.('Coefficients').('Estimate')(1));
% xlim([1 size(COM,2)]);

    f_COM_slope=model.('Coefficients').('Estimate')(2);
    c_r=model.Rsquared.Ordinary;
    A=coefCI(model,0.05);
    A(1,:)=[];
    if sum(sign(A))==0
        slope_sig=0;
    else
        slope_sig=1;
    end
%     model_summary=anova(model,'summary');
%     slope_sig=model_summary.('pValue')(2);

    
    s_x=1:length(f_skewness);
    s_y=f_skewness;
    s_x(isnan(f_skewness))=[];
    s_y(isnan(f_skewness))=[];
    skew_model=fitlm(s_x, s_y);
% plot(x,model.('Coefficients').('Estimate')(2).*x+ model.('Coefficients').('Estimate')(1));
% xlim([1 size(COM,2)]);

    f_skew_slope=skew_model.('Coefficients').('Estimate')(2);
    B=coefCI(skew_model,0.05);
    B(1,:)=[];
    if sum(sign(B))==0
        skew_sig=0;
    else
        skew_sig=1;
    end
    s_r=skew_model.Rsquared.Ordinary;
%     model_summary=anova(skew_model,'summary');
% %     skew_sig=model_summary.('pValue')(2);
%     if plot_on==1
%     figure;hold on;
%     plot(s_x,s_y,'o')
%     plot(s_x,s_x.*skew_model.('Coefficients').('Estimate')(2)+skew_model.('Coefficients').('Estimate')(1));
%     title([num2str(skew_sig) ' ' num2str(s_r)])
%     end


end
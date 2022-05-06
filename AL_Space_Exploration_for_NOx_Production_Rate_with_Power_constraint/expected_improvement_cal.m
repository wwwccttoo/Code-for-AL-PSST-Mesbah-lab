function EI_value = expected_improvement_cal(x,model,current_y)
    [current_mean,current_variance] = uq_evalModel(model,x);
    current_min = min(current_y);
    EI_value = (current_min-current_mean).*normcdf((current_min-current_mean)./sqrt(current_variance))+sqrt(current_variance).*normpdf((current_min-current_mean)./sqrt(current_variance));
end
function outputs = demof_for_GP(inputs)
    outputs = 2-sin(2.2*cos(inputs)+3)-1./(inputs.*inputs+5)-exp(-(inputs-5).^2/2);
    outputs = outputs/4;
end


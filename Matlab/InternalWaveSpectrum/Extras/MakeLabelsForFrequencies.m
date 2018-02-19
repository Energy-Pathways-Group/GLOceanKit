function labels = MakeLabelsForFrequencies(omega, f0, N0)

labels = cell(length(omega),1);
for i=1:length(omega)
    if round(omega(i)/f0) == 0.0
        labels{i} = '0';
    elseif round(omega(i)/f0) == 1
        labels{i} = 'f_0';
    elseif round(omega(i)/f0) == -1
        labels{i} = '-f_0';
    elseif round(omega(i)/f0) == round(N0/f0)
        labels{i} = 'N_0';
    elseif round(omega(i)/f0) == -round(N0/f0)
        labels{i} = '-N_0';
    else
        labels{i} = sprintf('%df_0',round(omega(i)/f0));
    end
end

end
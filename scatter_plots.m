function scatter_plots(obj, targetPath, paramPath, thresh)
    targetPathParts = strsplit(targetPath, '.');
    paramPathParts  = strsplit(paramPath, '.');
    % Iterate through each object instance in the array
    % Dynamically navigate down to the nested structure for this instance
    targetStruct = obj;
    for partIdx = 1:length(targetPathParts)
        if isprop(targetStruct, targetPathParts{partIdx}) || isfield(targetStruct, targetPathParts{partIdx})
            targetStruct = targetStruct.(targetPathParts{partIdx});
        else
            error(['Path ' targetPath ' does not exist in the object at index ' num2str(idx)]);
        end
    end
    paramStruct = obj;
    for partIdx = 1:length(paramPathParts)
        if isprop(paramStruct, paramPathParts{partIdx}) || isfield(paramStruct, paramPathParts{partIdx})
            paramStruct = paramStruct.(paramPathParts{partIdx});
        else
            error(['Path ' paramPath ' does not exist in the object at index ' num2str(idx)]);
        end
    end
    

    targetFields     = fieldnames(targetStruct);
    paramFields      = fieldnames(paramStruct);

    targetNumFields  = length(targetFields);
    numFigs  = targetNumFields; % Determine how many figures (each with 3x3 subplots) we need

    paramNumFields   = length(paramFields);
    factors  = 1:floor(sqrt(paramNumFields));
    validFactors = factors(mod(paramNumFields, factors)==0);
    corrFactors  = paramNumFields./validFactors;
    [~,idx] = min(abs(validFactors - corrFactors));
    numSubFigX = validFactors(idx);
    numSubFigY = corrFactors(idx);
    

    for i = 1:targetNumFields % Loop through each field
        figure; % Create a new figure for each set of 9 (or fewer) histograms
        targetField = targetFields{i};
        for j = 1:length(targetStruct.(targetField))
            if (~isstruct(targetStruct.(targetField).val))
                data1 = targetStruct.(targetField).val*targetStruct.(targetField).mult; % Extract data for current field
            end
        end
        for j = 1:paramNumFields
            paramField = paramFields{j};
            if (~isstruct(paramStruct.(paramField).val))
                data2 = paramStruct.(paramField).val*paramStruct.(paramField).mult;
                ax = subplot(numSubFigX, numSubFigY, j); % Create subplot
                scatter(data2, data1, 'filled');
                corr = corrcoef(data1, data2);
                if (abs(corr(1,2)) > thresh)
                    set(ax, 'Color', [0.7 0.9 0.7]);  % RGB triplet for very light green
                    hold on;
                    p = polyfit(data1, data2, 1);  % Fit a first-degree polynomial (linear)
                    yfit = polyval(p, data1);
                    plot(yfit, data1, '-r');  % Plot the regression line
                    text(max(get(gca, 'xlim')), max(get(gca, 'ylim')), ...
                        sprintf('Corr: %.2f\nSens: %.2f', corr(1,2), p(1)), ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
                elseif (abs(corr(1,2)) > thresh-0.1)
                    set(ax, 'Color', [1.0 1.0 0.8]);  % RGB triplet for very light green
                    text(max(get(gca, 'xlim')), max(get(gca, 'ylim')), ...
                        sprintf('Corr: %.2f', corr(1,2)), ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
                else
                    text(max(get(gca, 'xlim')), max(get(gca, 'ylim')), ...
                        sprintf('Corr: %.2f', corr(1,2)), ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
                end
                % fprintf("\ntargetField = %s, paramField = %s", targetField, paramField);
                ylabel(sprintf('%s [%s]', strrep(targetField, '_', '\_'), string(targetStruct.(targetField).units)), 'FontSize', 9);
                xlabel(sprintf('%s [%s]', strrep(paramField, '_', '\_'), string(paramStruct.(paramField).units)), 'FontSize', 9);
            end
        end
    end
end

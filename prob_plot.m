function prob_plot(obj, nestedPath)
    % Split the nested path into parts for dynamic field reference
    pathParts = strsplit(nestedPath, '.');
    
    % Iterate through each object instance in the array
    % Dynamically navigate down to the nested structure for this instance
    currentStruct = obj;
    for partIdx = 1:length(pathParts)
        if isprop(currentStruct, pathParts{partIdx}) || isfield(currentStruct, pathParts{partIdx})
            currentStruct = currentStruct.(pathParts{partIdx});
        else
            error(['Path ' nestedPath ' does not exist in the object at index ' num2str(idx)]);
        end
    end

    fields     = fieldnames(currentStruct); % Get all field names
    numFields  = length(fields); % Number of fields in the struc
    numFigures = ceil(numFields / 9); % Determine how many figures (each with 3x3 subplots) we need
    
    for f = 1:numFigures
        figure; % Create a new figure for each set of 9 (or fewer) histograms
        for i = 1:min(9, numFields - (f-1)*9) % Loop through each field
            subplot(3, 3, i); % Create subplot
            currentField = fields{(f-1)*9+i};
            for j = 1:length(currentStruct.(currentField))
                if (~isstruct(currentStruct.(currentField).val))
                    data = currentStruct.(currentField).val*currentStruct.(currentField).mult; % Extract data for current field
%                   spec_min = currentStruct.spec_min.(currentField);
%                   spec_max = currentStruct.spec_max.(currentField);
                end
            end
            % Calculate statistics
            mu = mean(data);
%           spec_min_avg = mean(spec_min);
%           spec_max_avg = mean(spec_max);
            sigma = std(data);
            medianVal = median(data);
            pctOutside = sum(data < mu - 3*sigma | data > mu + 3*sigma) / length(data) * 100;
            
            % Create probplot and add lines for mean ± 3σ
            probplot(data);
            hold on;
            %grid on;
            %xline(mu, 'LineWidth', 2);
%           xline(spec_min_avg, 'r', 'LineWidth', 2);
%           xline(spec_max_avg, 'r', 'LineWidth', 2);
            hold off;
            
            % Annotate the histogram
            %xlabel(sprintf('%s', strrep(currentField,'_', '\_')));
            xlabel(sprintf('%s [%s]', currentStruct.(currentField).name, string(currentStruct.(currentField).units)))
            %text(max(get(gca, 'xlim')), max(get(gca, 'ylim')), ...
            %    sprintf('Mean: %.2f\n3σ: %.2f\nMedian: %.2f\nOut>3σ: %.2f%%', mu, 3*sigma, medianVal, pctOutside), ...
            %    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
            %h = findobj(gca, 'Type', 'histogram');
            %h.FaceColor = string(currentStruct.(currentField).color);
        end
    end
end

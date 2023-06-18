%path to the REIMS data
dataPath = 'F:\Beatson\Beatson Epithelial Pellets REIMS';
%name of the data
dataName = 'V1_SI_data_New.mat';
%load the data
load([dataPath filesep dataName])

%set the mass range to make the comparison from
lowerMassRange = 100:100:700;
higherMassRange = 500:100:1200;

for ii = 1:length(lowerMassRange)
    for jj = 1:length(higherMassRange)
        if lowerMassRange(ii) < higherMassRange(jj)
            
            
            %APC KRAS: * magenta
            %APC: d cyan
            %APC PTEN: + red
            %APC KRAS PTEN: s blue
            %KRAS: o green
            %WT: o yellow
            
            %set the labels for male vs. female
            markerLabels = {'.','.','x', '.',...
                'x','x','.','x',...
                '.','x','x','x',...
                'x','x','.','.',...
                '.','x','x','x',...
                '.','.','x','x',...
                'x','.','.','.','x'};
            
            %set the labels for each genetic variant
            geneticLables = {'APC KRAS PTEN', 'APC PTEN', 'APC',...
                'APC KRAS', 'WT', 'KRAS'};
            %index of what row is what variant (number matches the index of
            %geneticLables)
            geneticIdx = [1 1 1 1 2 2 1 3 3 3 3 3 3 3 3 3 4 3 3 4 4 4 6 5 5 ...
                6 5 6 6];
            
            %set the colours for each of the matching labels
            blue = [1,2,3,4,7];
            cyan = [8,9,10,11,12,13,14,15,16,18,19];
            magenta = [17,20,21,22];
            yellow = [24,25,27];
            green = [23,26,28,29];
            %adjust for removal of the APC/PTEN
            blue(blue>6) = blue(blue>6)-2;
            cyan(cyan>6) = cyan(cyan>6)-2;
            magenta(magenta>6) = magenta(magenta>6)-2;
            yellow(yellow>6) = yellow(yellow>6)-2;
            green(green>6) = green(green>6)-2;
            
            %create a cell for each row of data and add the rgb for the
            %appropriate colour
            [colourLabels{blue}] = deal([0 0 1]);
            [colourLabels{green}] = deal([0 1 0]);
            [colourLabels{magenta}] = deal([1 0 1]);
            [colourLabels{yellow}] = deal([1 1 0]);
            [colourLabels{cyan}] = deal([0 1 1]);
            
            
            
            %load in data
            data = dataRepresentation2.data;
            %load in the spectral channels
            spectralChannels = dataRepresentation2.spectralChannels;
            
            %load the mask
            mask = dataRepresentation2.regionOfInterest.pixelSelection';
            %add an empty pixel on the end of the mask to fix indexing
            %problems
            mask(size(mask,1)+1,:) =0;
            
            
            %% remove APC/Pten data as this was not used in the end as it had only two pellets
            %remove the rows that match APC/Pten
            toRemove = [5 6];
            %get the rows in the mask and make false
            seg = mask;
            seg(:,toRemove) = 0;
            dataToKeep = seg(mask==1);
            data = data(dataToKeep,:);
            mask = seg;
            mask(:,toRemove) = [];
            mask(size(mask,1)+1,:) =0;
            
            %% Reduced mass range
            %get the selected mass range data to look at
            reducedMassRange = (double(spectralChannels<higherMassRange(jj)) + double(spectralChannels>lowerMassRange(ii))) ==2;
            data = data(:, reducedMassRange);
            
            % get the corresponding spectral channels
            
            spectralChannels = spectralChannels(reducedMassRange);
            
            % Scores = tsne(newData, [] ,3);
            tsneReduced = tsne(data,'Algorithm','exact', 'NumDimensions',3, 'Distance', 'correlation',...
                'NumPCAComponents', 0, 'Perplexity', 30);
            %% make the figure
            h = figure;
            hold all
            %remove temp2 variable if it exists
            try
                clear temp2;
            end
            %plot one of each type first for labelling the legend correctly
            for j = 1:length(unique(geneticIdx))
                %get the current genetic index
                idx = find(geneticIdx==j,1);
                %get the t-SNE data from that row
                for i = 1:3
                    image = double(mask);
                    image(image==1) = tsneReduced(:,i);
                    temp = image(:,idx);
                    temp2(i,:) = temp(temp~=0);
                end
                %plot that data
                d = scatter3(temp2(1,:), temp2(2,:), temp2(3,:), 1, colourLabels{idx}, '.');
                clear temp2
            end
            %plot all the data now that the legend will be correct
            for j = 1:size(mask,2)
                %go through each row and get the corresponding data
                for i = 1:3
                    image = double(mask);
                    image(image==1) = tsneReduced(:,i);
                    temp = image(:,j);
                    temp2(i,:) = temp(temp~=0);
                end
                
                %plot the data with the appropriate label
                if markerLabels{j} == 'x'
                    d = scatter3(temp2(1,:), temp2(2,:), temp2(3,:), 50, colourLabels{j}, markerLabels{j});
                else
                    d = scatter3(temp2(1,:), temp2(2,:), temp2(3,:), 200, colourLabels{j}, markerLabels{j});
                end
                %set the label line size
                set(d, 'LineWidth', 2)
                clear temp2
            end
            %add axes and legends
            xlabel({['t-SNE dimension 1']})
            ylabel({['t-SNE dimension 2']})
            zlabel({['t-SNE dimension 3']})
            legend(geneticLables)
            
            
            
            %save the figure
            saveas(h, [cd filesep 'labelled tSNE mz' num2str(lowerMassRange(ii)) ' to' num2str(higherMassRange(jj)) ], 'fig')
            saveas(h, [cd filesep 'labelled tSNE mz' num2str(lowerMassRange(ii)) ' to' num2str(higherMassRange(jj)) ], 'jpg')
            
            %close the figure
            close(h)
        end
    end
end
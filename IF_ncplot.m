function IF_ncplot(plotflag, resvec, scatx, scaty, pathlist_labels, ChannelLabel, calclbl, file, slashtype)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
colors = brewermap(length(file), 'Set2');

%% Per-Treatment Scatter plot
if plotflag(1,2)==1
    treatmentNo = plotflag(2,1);
    fig2 = figure
    scatter(scatx, scaty, 'filled', ...
            'MarkerFaceColor', colors(treatmentNo,:), ...
            'MarkerFaceAlpha',3/9)
    xlim([0 255])
    ylim([0 255])
    xlabel([ChannelLabel{1,1} ' ' ChannelLabel{2,1} ' Intensity (a.u.)'])
    ylabel([ChannelLabel{1,2} ' ' ChannelLabel{2,2} ' Intensity (a.u.)'])
    title(pathlist_labels)
    legend([pathlist_labels '; No. of cells: ' num2str(length(scatx))])
    
    if ~isdir([file(1).maindir slashtype 'Scatter'])
        mkdir([file(1).maindir slashtype 'Scatter']);
    end
    print(fig2,[file(1).maindir slashtype 'Scatter' slashtype 'sca_' pathlist_labels], '-painters', '-dsvg','-r200')
    close gcf
end


%% Consolidated Scatter Plot
if plotflag(1,3)==1
    treatmentNo = plotflag(2,1);
    treatments = nonzeros(plotflag(3,1:end))';
    fig3 = figure
    hold on
    
    for ctr1 = treatments
        xvals = scatx{ctr1}(:,1);
        yvals = scaty{ctr1}(:,1);
        scatter(xvals, yvals, 'filled', ...
                'MarkerFaceColor', colors(ctr1,:), ...
                'MarkerFaceAlpha',3/9)
        legendary{ctr1} = [char(pathlist_labels(ctr1)) '; No. of cells: ' num2str(length(xvals))];
    end
    xlim([0 255])
    ylim([0 255])
    xlabel([ChannelLabel{1} ' intensity (a.u.)'])
    ylabel([ChannelLabel{2} ' intensity (a.u.)'])
    title('Consolidated Scatter')
    %Clean up empty legends:
    cleaner = cellfun(@(x) ~isempty(x), legendary);
    legendary = legendary(cleaner);
    legend(legendary)
    
    if ~isdir([file(1).maindir slashtype 'Con_Scatter'])
        mkdir([file(1).maindir slashtype 'Con_Scatter']);
    end
    print(fig3,[file(1).maindir slashtype 'Con_Scatter' slashtype 'ConsolidatedScat2'], '-painters', '-dsvg','-r200')
    close gcf
end




%% Box-Plot
%Changing the sort labels to group data by treatment


if plotflag(1,1)==1
    for ctr1 = 1:length(resvec)
        resvec{1,ctr1}(:,2) = ctr1;
    end
    tempvec = cell2mat(resvec');
    
    
    fig1 = figure
    notBoxPlot(tempvec(:,1), tempvec(:,2))
    set(gca, 'XTickLabel', pathlist_labels)
    title(['Median ' char(ChannelLabel) ' ' calclbl ' Intensity Per Cell'])
    ylabel('Median Intensity')
    if ~isdir([file(1).maindir slashtype 'BoxPlot'])
        mkdir([file(1).maindir slashtype 'BoxPlot']);
    end
    print(fig1,[file(1).maindir slashtype 'BoxPlot' slashtype 'box_' calclbl '_' char(ChannelLabel)], '-painters', '-dsvg','-r200')
    close gcf
end
end





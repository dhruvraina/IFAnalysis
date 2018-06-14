function IF_ncplot(plotflag, resvec, scatx, scaty, pathlist_labels, ChannelLabel, calclbl, file, file2, slashtype, lims)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
colors = brewermap(length(file), 'Set2');

switch plotflag.type
    %% Per-Treatment Scatter plot
    case('SingleScatter')
        treatmentNo = plotflag.singTreat;
        fig2 = figure
        scatter(scatx, scaty, 'filled', ...
            'MarkerFaceColor', colors(treatmentNo,:), ...
            'MarkerFaceAlpha',3/7)
        xlim(lims.x)
        ylim(lims.y)
        xlabel([ChannelLabel{1,1} ' ' ChannelLabel{2,1} ' Intensity (a.u.)'])
        ylabel([ChannelLabel{1,2} ' ' ChannelLabel{2,2} ' Intensity (a.u.)'])
        title(pathlist_labels)
        h = legend([pathlist_labels '; No. of cells: ' num2str(length(scatx))])
        h.FontSize = 14
        
        
        if plotflag.corrprint ==1
            hold on
            [clinx cliny] = corrline(scatx, scaty);
            [pearsn pval] = corr(scatx, scaty);
            rsq = pearsn^2;
            plot(clinx, cliny)
            plotstrg = {['r^2 = ' num2str(rsq)], ['p-value = ' num2str(pval)]};
            dim = [.2 .5 .3 .3];
            annotation('textbox',dim,'String',plotstrg,'FitBoxToText','on');
        end
        
        set(gca,'FontSize', 14)
        
        if ~isdir([file(1).maindir slashtype 'Scatter'])
            mkdir([file(1).maindir slashtype 'Scatter']);
        end
        print(fig2,[file(1).maindir slashtype 'Scatter' slashtype file2.savename 'sca_' pathlist_labels], '-painters', '-dsvg','-r200')
        close gcf
        
        
        
        %% Consolidated Scatter Plot
    case('ConScatter')
        treatments = nonzeros(plotflag.conTreat)';
        fig3 = figure
        hold on
        
        for ctr1 = treatments
            xvals = scatx{ctr1}(:,1);
            yvals = scaty{ctr1}(:,1);
            scatter(xvals, yvals, 'filled', ...
                'MarkerFaceColor', colors(ctr1,:), ...
                'MarkerFaceAlpha',3/7)
            legendary{ctr1} = [char(pathlist_labels(ctr1)) '; No. of cells: ' num2str(length(xvals))];
        end
        xlim(lims.x)
        ylim(lims.y)
        xlabel([ChannelLabel{1,1} ' intensity (a.u.)'])
        ylabel([ChannelLabel{1,2} ' intensity (a.u.)'])
        title('Consolidated Scatter')
        %Clean up empty legends:
        cleaner = cellfun(@(x) ~isempty(x), legendary);
        legendary = legendary(cleaner);
        h = legend(legendary)
        h.FontSize = 14
        
        set(gca,'FontSize', 14)
        if ~isdir([file(1).maindir slashtype 'Con_Scatter'])
            mkdir([file(1).maindir slashtype 'Con_Scatter']);
        end
        print(fig3,[file(1).maindir slashtype 'Con_Scatter' slashtype file2.savename 'ConScat2'], '-painters', '-dsvg','-r200')
        close gcf
        
        
        
        
        
        %% Box-Plot
        %Changing the sort labels to group data by treatment
        
        
    case ('boxplot')
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
        set(gca,'FontSize', 14)
        print(fig1,[file(1).maindir slashtype 'BoxPlot' slashtype 'box_' calclbl '_' char(ChannelLabel)], '-painters', '-dsvg','-r200')
        close gcf
end
end





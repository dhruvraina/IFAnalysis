%IF_ncplot
%Description: Plotter function for IF_gui. Called by IF_gui_backend.m.
%Dependencies: brewermap, notBoxPlot, UnivarScatter, cleanNames
%Author: draina
%Last Edit: 180822 

function IF_ncplot(plotflag, resvec, scatx, scaty, pathlist_labels, ChannelLabel, calclbl, file, lims)

colors = brewermap(length(file.treatmentLabels), 'Set2');
switch plotflag.type
    
   
% --------------  A. Single Treatment Scatter Plots  --------:
    case('SingleScatter')
        treatmentNo = plotflag.singTreat;
        fig1        = figure
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
        
        %Plot line of best fit + Pearson correlation:
        if plotflag.corrprint ==1
            hold on
            [clinx cliny] = corrline(scatx, scaty);
            [pearsn pval] = corr(scatx, scaty);
            rsq           = pearsn^2;
            plot(clinx, cliny)
            plotstrg      = {['r^2 = ' num2str(rsq)], ['p-value = ' num2str(pval)]};
            dim           = [.2 .5 .3 .3];
            annotation('textbox',dim,'String',plotstrg,'FitBoxToText','on');
        end
        
        set(gca,'FontSize', 14)
        fig1.PaperUnits    = 'inches';
        fig1.PaperPosition = [0 0 8 8];
        
        if ~isdir([file.outpath file.slashtype 'Scatter'])
            mkdir([file.outpath file.slashtype 'Scatter']);
        end
        
        switch plotflag.imageFormat
            case 'svg'
                print(fig1,[file.outpath file.slashtype 'Scatter' file.slashtype 'sca_' pathlist_labels], '-painters', '-dsvg','-r200')
            case'png'
                print(fig1,[file.outpath file.slashtype 'Scatter' file.slashtype 'sca_' pathlist_labels], '-painters', '-dpng','-r200')
        end
        close gcf
        
        
        
% --------------  B. Multiple Treatments In One Scatter  --------:
    case('ConScatter')
        treatments = nonzeros(plotflag.conTreat)';
        fig1       = figure
        hold on
        
        %Plot each treatment sequentially:
        for ctr1  = treatments
            xvals = scatx{ctr1}(:,1);
            yvals = scaty{ctr1}(:,1);
            scatter(xvals, yvals, 'filled', ...
                'MarkerFaceColor', colors(ctr1,:), ...
                'MarkerFaceAlpha',3/7)
            legendary{ctr1} = [pathlist_labels{ctr1} '; No. of cells: ' num2str(length(xvals))];
        end
        
        xlim(lims.x)
        ylim(lims.y)
        xlabel([ChannelLabel{1,1} ' intensity (a.u.)'])
        ylabel([ChannelLabel{1,2} ' intensity (a.u.)'])
        title('Consolidated Scatter')
        
        %Clean up empty legends + set props:
        cleaner    = cellfun(@(x) ~isempty(x), legendary);
        legendary  = legendary(cleaner);
        h          = legend(legendary)
        h.FontSize = 14
        
        set(gca,'FontSize', 14)
        fig1.PaperUnits    = 'inches';
        fig1.PaperPosition = [0 0 8 8];
        if ~isdir([file.outpath file.slashtype 'Con_Scatter'])
            mkdir([file.outpath file.slashtype 'Con_Scatter']);
        end
        
        switch plotflag.imageFormat
            case 'svg'
                print(fig1,[file.outpath file.slashtype 'Con_Scatter' file.slashtype 'ConScat'], '-painters', '-dsvg','-r200')
            case'png'
                print(fig1,[file.outpath file.slashtype 'Con_Scatter' file.slashtype 'ConScat'], '-painters', '-dpng','-r200')
        end

        close gcf

        
        
        
% --------------  C. Boxplot - notBoxPlot  --------:
    case ('boxplot1')
        %Unpack resvec
        for ctr1 = 1:length(resvec)
            resvec{1,ctr1}(:,2) = ctr1;
        end
        tempvec = cell2mat(resvec');
        
        %Plotting:
        fig2 = figure
        notBoxPlot(tempvec(:,1), tempvec(:,2))
        set(gca, 'XTickLabel', pathlist_labels)
        title(['Median ' char(ChannelLabel) ' ' calclbl ' Intensity Per Cell'])
        ylabel('Median Intensity')
        
        fig2.PaperUnits    = 'inches';
        fig2.PaperPosition = [0 0 10 6];
        set(gca,'FontSize', 14)
        if ~isdir([file.outpath file.slashtype 'BoxPlot1'])
            mkdir([file.outpath file.slashtype 'BoxPlot1']);
        end
        
        switch plotflag.imageFormat
            case 'svg'
                print(fig2,[file.outpath file.slashtype 'BoxPlot1' file.slashtype 'box1_' calclbl '_' ChannelLabel], '-painters', '-dsvg','-r200')
            case'png'
                print(fig2,[file.outpath file.slashtype 'BoxPlot1' file.slashtype 'box1_' calclbl '_' ChannelLabel], '-painters', '-dpng','-r200')
        end
        close gcf

        
        
        
% --------------  D. Boxplot - UnivarScatter  --------:
    case('boxplot2')
        %unpack resvec
        for ctr1 = 1:length(resvec)
            resvec{1,ctr1}(:,2) = ctr1;
        end
        tempvec = cell2mat(resvec');
        
        %Make an empty matrix with nans
        maxLength = max(cell2mat(cellfun(@(x) length(x), resvec, 'UniformOutput', 0)));
        nTreatments = size(resvec,2);
        dataArray = nan(maxLength, nTreatments);
        
        %Fill in values from tempvec
        for cc = 1:nTreatments
            tvec2 = tempvec(tempvec(:,2)==cc,1);
            dataArray(1:length(tvec2),cc) = tvec2;
            clear tvec2
        end
        
        %Plotting:
        fig2 = figure, UnivarScatter(dataArray, 'RangeCut', 30, 'Label', pathlist_labels, 'Whiskers', 'lines', 'PointSize', 15);
        ylabel('Median Intensity')
        title(['Median ' char(ChannelLabel) ' ' calclbl ' Intensity Per Cell'])
        
        fig2.PaperUnits = 'inches';
        fig2.PaperPosition = [0 0 10 6];
        set(gca,'FontSize', 14)
        if ~isdir([file.outpath file.slashtype 'BoxPlot2'])
            mkdir([file.outpath file.slashtype 'BoxPlot2']);
        end
        
        switch plotflag.imageFormat
            case 'svg'
                print(fig2,[file.outpath file.slashtype 'BoxPlot2' file.slashtype 'univ_'  calclbl '_' ChannelLabel], '-painters', '-dsvg','-r200')
            case 'png'
                print(fig2,[file.outpath file.slashtype 'BoxPlot2' file.slashtype 'univ_'  calclbl '_' ChannelLabel], '-painters', '-dpng','-r200')
        end
        close gcf
end
end









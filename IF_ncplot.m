%IF_ncplot
%Description: Plotter function for IF_gui. Called by IF_gui_backend.m.
%Dependencies: brewermap, notBoxPlot, UnivarScatter, cleanNames
%Author: draina
%Last Edit: 180822

function IF_ncplot(plotflag, resvec, scatx, scaty, scatz, pathlist_labels, ChannelLabel, calclbl, file, lims, calcs)
%Using Brewermap
colors = brewermap(length(file.treatmentLabels), 'Spectral');  %'Set2'
addpath([file.codeparent file.slashtype 'ExchangeTools']);
addpath([file.outpath file.slashtype]);
%Using UniformPercepCols from matplotlib
%col1   = viridis();
col1   = magma();
col2   = floor(length(col1)/length(file.treatmentLabels));
colors = flipud(col1(1:col2:end,:));

switch plotflag.type
    
    
    % --------------  A. Single Treatment Scatter Plots  --------:
    case('SingleScatter')
        treatmentNo = plotflag.singTreat;
        fig1        = figure
        
        if plotflag.margDist %Inbuilt fx scatterhist lacks MarkerCol and FaceAlpha, so use subplot
            
            %Build the marginal distributions
            nbins = 20;
            xbins = (lims.x(2)-lims.x(1))/nbins;
            ybins = (lims.y(2)-lims.y(1))/nbins;
 
            %Plot marginal distributions
            sbplt1 = subplot(2,2,2)
            histogram(scaty, lims.y(1):ybins:lims.y(2), 'Normalization', 'probability', 'FaceColor', colors(treatmentNo,:));
            view([90 -90])
            sbplt2 = subplot(2,2,3)
            histogram(scatx, lims.x(1):xbins:lims.x(2), 'Normalization', 'probability', 'FaceColor', colors(treatmentNo,:));
            view([0,-90])
            
            %Plot scatter last so its the gca for xlim ylim etc. later
            subplot(2,2,1)
            scatter(scatx, scaty, 'filled', ...
                'MarkerFaceColor', colors(treatmentNo,:), ...
                'MarkerFaceAlpha',3/7)
            
            %Setting the aspect ratio and repositioning the subplots
            sbplt1.Position = [sbplt1.Position(1), sbplt1.Position(2)    sbplt1.Position(3)/4, sbplt1.Position(4)  ]
            sbplt2.Position = [sbplt2.Position(1), sbplt2.Position(2)*3,   sbplt2.Position(3), sbplt2.Position(4)/4]
        else
            scatter(scatx, scaty, 'filled', ...
                'MarkerFaceColor', colors(treatmentNo,:), ...
                'MarkerFaceAlpha',3/7)
        
        if plotflag.logscale
            set(gca, 'xscale', 'log')
            set(gca, 'yscale', 'log')
            hold on
            logOffst = 1; %plotting 0 on log axis doesn't work, thresholds don't show.
        end
            
        %Setting the cutoff lines
        if calcs.Ratios
            hold on
            plot([1 1].*calcs.XAxisThresh, lims.y+logOffst, '--k')
            plot(lims.x+logOffst, [1 1].*calcs.YAxisThresh, '--k')
        end
        
        end
        
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
        
        %Load figure props for different projects:
        ff.fontSize      = 14;
        ff.PaperPosition = [0 0 10 10];
        ff.PaperUnits    = 'inches';
        ff.pointSize     = 200;
        addProps = 0;
        if addProps
            ff = figprops; %Potentially save various properties as a function in the outpath
        end
        %get markersize handles and adjust
        handleMarker = findall(gca, 'marker', 'o');
        set(handleMarker, 'sizedata', ff.pointSize);
        
        set(gca,'FontSize', ff.fontSize)
        set(gca, 'TickDir', 'out');
        set(gcf, 'Units', ff.PaperUnits, 'Position', ff.PaperPosition);
       
        
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
        
        %Create the subplots out of loop to keep from overwriting.
        %each subplot also needs its own explicit hold!
        if plotflag.margDist && ~plotflag.tdplot
            sbplt3 = subplot(2,2,1)
            hold on
            sbplt1 = subplot(2,2,2)
            view([90 -90])
            hold on
            sbplt2 = subplot(2,2,3)
            view([0,-90])
            hold on
            
            %Setting the aspect ratio and repositioning the subplots
            sbplt1.Position = [sbplt1.Position(1), sbplt1.Position(2)    sbplt1.Position(3)/4, sbplt1.Position(4)  ]
            sbplt2.Position = [sbplt2.Position(1), sbplt2.Position(2)*3,   sbplt2.Position(3), sbplt2.Position(4)/4]
            
%             %hard coding to save as fig
%             plotflag.imageFormat = 'fig';
        end
        
        cnt1 = 1;
        
        %Plot each treatment sequentially:
        for ctr1  = treatments
              xvals           = scatx{ctr1}(:,1);
              yvals           = scaty{ctr1}(:,1);
              zvals           = scatz{ctr1}(:,1);
              legendary{cnt1} = [pathlist_labels{ctr1} '; No. of cells: ' num2str(length(xvals))];
              cnt1            = cnt1+1;
              
            if plotflag.margDist && ~plotflag.tdplot
                %Build the marginal distributions
                nbins = 20;
                xbins = (lims.x(2)-lims.x(1))/nbins;
                ybins = (lims.y(2)-lims.y(1))/nbins;

                %Plot marginal distributions
                histogram(sbplt1, yvals, lims.y(1):ybins:lims.y(2),'Normalization', 'probability', 'FaceColor', colors(ctr1,:)); %changed from countdensity
                histogram(sbplt2, xvals, lims.x(1):xbins:lims.x(2),'Normalization', 'probability', 'FaceColor', colors(ctr1,:));
                
                %Plot scatter
                scatter(sbplt3, xvals, yvals, 'filled', ...
                    'MarkerFaceColor', colors(ctr1,:), ...
                    'MarkerFaceAlpha',3/7)
                
                %Set scatter as gca so axis props will be correctly
                %assigned later
                figure(fig1)
                subplot(sbplt3)


            elseif plotflag.tdplot
                %3d Plot
                scatter3(xvals, yvals, zvals, 'filled', ...
                    'MarkerFaceColor', colors(ctr1,:), ...
                    'MarkerFaceAlpha', 5/7)
                view(18,33)
                zlim(lims.z)
                zlabel([ChannelLabel{1,3} ' intensity (a.u.)'])
                grid on
                
            else
                scatter(xvals, yvals, 100, 'filled', ...
                    'MarkerFaceColor', colors(ctr1,:), ...
                    'MarkerFaceAlpha',3/7)

            end
        end
        
        %Draw cutoff lines:
        if calcs.Ratios
            plot([1 1].*calcs.XAxisThresh, lims.y, '--k')
            plot(lims.x, [1 1].*calcs.YAxisThresh, '--k')
        end
        
        xlim(lims.x)
        ylim(lims.y)
        xlabel([ChannelLabel{1,1} ' intensity (a.u.)'])
        ylabel([ChannelLabel{1,2} ' intensity (a.u.)'])
        title('Consolidated Scatter')
        
        
        %Clean up empty legends + set props:
        cleaner    = cellfun(@(x) ~isempty(x), legendary);
        legendary  = legendary(cleaner);
        h          = legend(legendary);
        h.FontSize = 14;
        
        set(gca, 'TickDir', 'out');
        set(gca,'FontSize', 14)
        fig1.PaperUnits    = 'inches';
        fig1.PaperPosition = [0 0 10 10];
        if ~isdir([file.outpath file.slashtype 'Con_Scatter'])
            mkdir([file.outpath file.slashtype 'Con_Scatter']);
        end
        
        switch plotflag.imageFormat
            case 'svg'
                print(fig1,[file.outpath file.slashtype 'Con_Scatter' file.slashtype 'ConScat'], '-painters', '-dsvg','-r200')
            case 'png'
                print(fig1,[file.outpath file.slashtype 'Con_Scatter' file.slashtype 'ConScat'], '-painters', '-dpng','-r200')
            case 'fig'
                savefig(fig1, [file.outpath file.slashtype 'Con_Scatter' file.slashtype 'ConScat'])
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
        title(['Mean ' char(ChannelLabel) ' ' calclbl ' Intensity Per Cell'])
        ylabel('Mean Intensity')
        
        ylim([lims.boxmin lims.boxmax])
        fig2.PaperUnits    = 'inches';
        fig2.PaperPosition = [0 0 10 6];
        
        set(gca, 'TickDir', 'out');
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
        maxLength   = max(cell2mat(cellfun(@(x) length(x), resvec, 'UniformOutput', 0)));
        nTreatments = size(resvec,2);
        dataArray   = nan(maxLength, nTreatments);
        
        %Fill in values from tempvec
        for cc = 1:nTreatments
            tvec2 = tempvec(tempvec(:,2)==cc,1);
            
            if all(tvec2)==0
                tvec2(1) = 0.00001;
                msgbox('Warning: setting one value from 0 to 0.0001 so the plots work properly')
            elseif max(isinf(tvec2))
                tvec2(isinf(tvec2)) = 0.00001;
                msgbox('Warning: setting one value from Inf to 0.0001 so plots work properly')
            end
            
            dataArray(1:length(tvec2),cc) = tvec2;
            clear tvec2
        end
        
        %Plotting:
        figWide = length(pathlist_labels)*2.3;
        figTall = 7;
        
        fig2 = figure, UnivarScatter(dataArray, 'RangeCut', 30, 'Label', pathlist_labels, 'Whiskers', 'lines', 'PointSize', 12);
        ylabel('Mean Intensity')
        title(['Mean ' char(ChannelLabel) ' ' calclbl ' Intensity Per Cell'])
        
        ylim([lims.boxmin lims.boxmax])
        fig2.PaperUnits    = 'inches';
        fig2.PaperPosition =  [0 0 figWide figTall];
      %  annotation('textbox', [0.3 0.3 .9 .7], 'String', num2str(lims.boxmean), 'FitBoxToText','on')
        
      %Mean Annotations:
        for cc2 = 1:size(dataArray, 2)
           labels_stacked = num2str(floor(lims.boxmean(cc2).*10)./10);
           hText          = text(cc2, lims.maxLabel(cc2), labels_stacked);
           set(hText,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center','FontSize',14);
        end

        set(gca, 'TickDir', 'out');
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
        
         fileID = fopen([file.outpath file.slashtype file.experimentName '_Reformat.txt'],'w');
         fprintf(fileID, '%1$s\t %2$s\t %3$s\t %4$s\t %5$s\t %6$s\t %7$s\r\n', 'Treatment', 'Mean', 'St.Dev', 'Channel', 'Measurement','Number of Cells', 'Experiment');
         
         for tnum = 1:size(pathlist_labels,2)
             fprintf(fileID,'%1$s\t', pathlist_labels{tnum});
             fprintf(fileID,'%4.2f\t %4.2f\t', lims.boxmean(tnum), lims.boxStd(tnum));
             fprintf(fileID, '%1$s\t', ChannelLabel);
             fprintf(fileID, '%1$s\t', calclbl);
             fprintf(fileID, '%1$s\t', num2str(length(resvec{tnum})));
             fprintf(fileID, '%1$s\r\n', file.experimentName);
         end

        % --------------  E. Violin Plots  -----------------------: 
    case('violin')
        %unpack resvec
        tempvec = cellfun(@(x) x(:,1), resvec, 'UniformOutput', 0);
        
        %Make an empty matrix with nans
        maxLength   = max(cell2mat(cellfun(@(x) length(x), resvec, 'UniformOutput', 0)));
        nTreatments = size(resvec,2);
        
        %Colors:
        for nn = 1:nTreatments
            colCell{nn} = colors(nn,:);
        end
        
        figN = figure
        violinPlot(tempvec,'addSpread',true,'showMM',6, 'color', colCell); 
        violinPlot(tempvec,'xyOri','flipped','histOri','right','showMM',6),

        %Set markerSize in plotSpread.m
        %Set Quartile Markers size in violinPlots.m
        
        figWide = nTreatments*2.3;
        figTall = 7;
        
        ylabel('Mean Intensity')
        title(['Mean ' char(ChannelLabel) ' ' calclbl ' Intensity Per Cell'])
        
        ylim([lims.boxmin lims.boxmax])
        figN.PaperUnits    = 'inches';
        figN.PaperPosition =  [0 0 figWide figTall];
        
      %Mean Annotations:
        for cc2 = 1:nTreatments
           labels_stacked = num2str(floor(lims.boxmean(cc2).*10)./10);
           hText          = text(cc2, lims.maxLabel(cc2), labels_stacked);
           set(hText,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center','FontSize',14);
        end

        set(gca, 'TickDir', 'out');
        set(gca,'FontSize', 14)
        set(gca, 'XTickLabel', pathlist_labels);

        if ~isdir([file.outpath file.slashtype 'BoxPlot2'])
            mkdir([file.outpath file.slashtype 'BoxPlot2']);
        end
        
        switch plotflag.imageFormat
            case 'svg'
                print(figN,[file.outpath file.slashtype 'BoxPlot2' file.slashtype 'viol_'  calclbl '_' ChannelLabel], '-painters', '-dsvg','-r200')
            case 'png'
                print(figN,[file.outpath file.slashtype 'BoxPlot2' file.slashtype 'viol_'  calclbl '_' ChannelLabel], '-painters', '-dpng','-r200')
        end
        close gcf
        
               
        
        % --------------  F. Stacked Bar Plot for Ratios  --------:        
    case('stacked')
        %Redfining ChannelLabel as a struct for the stacked bar plots since they
        %require a different order of the data.
        resvec = resvec';
        fig2 = figure
        barWide = 10;
        barTall = 10;
        horzBar = 0;
        vertBar = 1;
        
        %Horizontal Bar:
        if horzBar
            h = barh(resvec, 0.9, 'stacked')
            set(gca, 'YTick', 1:length(pathlist_labels), 'YTickLabels', pathlist_labels)
            barTall = length(pathlist_labels)*2;
            xName   = 'Fraction';
            yName   = '';
            
            % Print the text labels inside the bar graph
            for i=1:size(resvec,1)
                for j=1:size(resvec,2)
                    if resvec(i,j)>0.1                                           %Don't print values less than 10% due to label space constraints
                        labels_stacked = num2str((floor(resvec(i,j)*100))/100);  %Round digits to nearest lower full percent
                        hText          = text(sum(resvec(i,1:j),2), i, labels_stacked);
                        set(hText,'HorizontalAlignment', 'right','FontSize',14, 'Color','w');
                    end
                end
            end
        end
        
        %Vertical Bar:
        if vertBar
            h = bar(resvec, 0.9, 'stacked')
            set(gca, 'XTick', 1:length(pathlist_labels), 'XTickLabels', pathlist_labels)
            barWide = length(pathlist_labels)*2;
            xName   = '';
            yName   = 'Fraction';
            
          % Print the text labels inside the bar graph
            for i=1:size(resvec,1)
                for j=1:size(resvec,2)
                    if resvec(i,j)>0.06                                          %Don't print values less than x% due to label space constraints
                        labels_stacked = num2str((floor(resvec(i,j)*100))/100);  %Round digits to nearest lower full percent
                        hText          = text(i, sum(resvec(i,1:j),2), labels_stacked);
                        set(hText,'VerticalAlignment', 'top', 'HorizontalAlignment', 'center','FontSize',14, 'Color','w');
                    end
                end
            end
            
        end
        %Setting custom colours:
        colmap = [217 111 171;  %Pink     GATA
                  210 188  44;  %Yellow   DoublePositive
                    0 168 110;  %Green    NANOG
                  140 198 236]; %Blue     DoubleNegative
        colmap = colmap./256;   %convert from rgb space to [0 1]
        
        h(1).FaceColor = colmap(1,:)
        h(2).FaceColor = colmap(2,:)
        h(3).FaceColor = colmap(3,:)
        h(4).FaceColor = colmap(4,:)
        
        %legend(ChannelLabel.xChan, ChannelLabel.dPos, ChannelLabel.yChan, ChannelLabel.dNeg)
        xlabel(xName)
        ylabel(yName)
        title(ChannelLabel.title)
        
        fig2.PaperUnits = 'inches';
        fig2.PaperPosition = [0 0 barWide barTall];
        set(gca,'FontSize', 14)
        if ~isdir([file.outpath file.slashtype 'StackedBar'])
            mkdir([file.outpath file.slashtype 'StackedBar']);
        end
        
        switch plotflag.imageFormat
            case 'svg'
                print(fig2,[file.outpath file.slashtype 'StackedBar' file.slashtype 'RatioStacked'], '-painters', '-dsvg','-r200')
            case 'png'
                print(fig2,[file.outpath file.slashtype 'StackedBar' file.slashtype 'RatioStacked'], '-painters', '-dpng','-r200')
        end
        
        fileID = fopen([file.outpath file.slashtype file.experimentName '_Reformat.txt'],'w');
        fprintf(fileID, '%1$s\t %2$s\t %3$s\t %4$s\t %5$s\r\n', 'Treatment', ChannelLabel.xChan, ChannelLabel.dPos, ChannelLabel.yChan, ChannelLabel.dNeg);
        
        for tnum = 1:size(pathlist_labels,2)
            fprintf(fileID,'%1$s\t', pathlist_labels{tnum});
            fprintf(fileID,'%5.4f\t %5.4f\t %5.4f\t %5.4f\r\n', resvec(tnum, 1:4));
        end
        fprintf(fileID, '%1$s\t', ChannelLabel.title);
        fprintf(fileID,'%5.4f\t %5.4f\t %5.4f\t %5.4f\r\n', zeros(1, 4));
            

    
        close gcf
        
        
end
end









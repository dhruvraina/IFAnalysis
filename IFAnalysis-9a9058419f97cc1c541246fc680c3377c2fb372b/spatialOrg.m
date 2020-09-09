%TBD: extend to 4 state variables for DN and DP.
%States 0- unclassified, 1-Nanog, 2-Gata6

%% Init
addpath(genpath('/Users/draina/Documents/Code/MATLAB/ExchangeTools/'));

testChanA   = resvec_calc4; %List the two channels that need testing
testChanB   = resvec_calc2; %List the two channels that need testing
labChanA    = [inputs.ch4lab 'intensity (log a.u.)'];
labChanB    = [inputs.ch2lab 'intensity (log a.u.)'];
nTreatments = size(resvec_another, 2);
fitComp     = 2;
lims.x      = [5 12];
lims.y      = [5 12];
%colormaps:
cols        = linspecer(nTreatments);    %GMfit
cols2       = brewermap(130, '*YlGnBu'); %GMfit gradient
cols2       = cols2(1:110,:); %Remove *very* yellow points
colmap      = [217 111 171;  %Pink     GATA
    210 188  44;  %Yellow   DoublePositive
    0 168 110;  %Green    NANOG
    140 198 236]; %Blue     DoubleNegative
colmap      = colmap./256;   %convert from rgb space to [0 1]
colGrad     = colorGradient(colmap(1,:), colmap(3,:), 101); %Gradient from G+ to N+


%% Assign state variables:
%Currently doing a simple gmm fit. Ask Angel how he factored in DP and DNs
for tt = outputs.limconscatID
    %Init loop vars:
    testVecA = [];
    testVecB = [];
    testVecT = [];
    
    %Log transform
    testVecA = log(testChanA{1, tt}(:,1));
    testVecB = log(testChanB{1, tt}(:,1));
    
    %fit gaussian mixture model
    testVecT = horzcat(testVecA, testVecB);
    gmfit    = fitgmdist(testVecT, fitComp, 'RegularizationValue', 0.15); %0.3
    post     = posterior(gmfit, testVecT);
    
    %Build complex figure
    figure
    p = panel();
    p.pack('h', {2/3 1/3}) %Divide main panel width
    p.fontsize = 12;
    
    %addressing:
    p1 = p(1);
    p2 = p(2);
    
    %pack again to include colorbar:
    p2.pack('h', {0.9 0.1});
    p2(1).pack(2,1);
    
    
    %Plot PDFs
    p1.select()
    scatter(testVecA, testVecB, 15, 'filled', 'MarkerFaceColor', [cols(3,:)])
    xlim(lims.x)
    ylim(lims.y)
    hold on
    gmPDF = @(x,y)reshape(pdf(gmfit,[x(:) y(:)]),size(x));                                      %Thanks @ A.S
    fc = fcontour(gmPDF,[[lims.x], [lims.y]],'LevelList',[0:0.01:0.8], 'LineColor', 'flat');    %Thanks @ A.S.
    xlabel(labChanA)
    ylabel(labChanB)
    title([file.treatmentLabels{tt} ' pdf of gm fit. P(0.9) is identity thresh'])
    caxis([0 1])
    colormap(cols2)
    
    %Component with highest yaxis mu is N+, i.e. stateVar = 1
    [~, nanIdx]= max(gmfit.mu);
    [~, gatIdx]= min(gmfit.mu);
    nanIdx     = nanIdx(2);
    gatIdx     = gatIdx(2);
    
    %Plot N+ prediction
    p2(1,1,1).select()
    scatter(testVecA, testVecB, 10, post(:,nanIdx), 'filled')
    title('Nanog+ prediction (state 1)')
    xlim(lims.x)
    ylim(lims.y)
    xlabel(labChanA)
    ylabel(labChanB)
    caxis([0 1])
    colormap(cols2)
    
    %Plot G+ prediction
    p2(1,2,1).select()
    scatter(testVecA, testVecB, 10, post(:,gatIdx), 'filled')
    title('Gata+ prediction (state 2)')
    xlim(lims.x)
    ylim(lims.y)
    xlabel(labChanA)
    ylabel(labChanB)
    caxis([0 1])
    colormap(cols2)
    
    %Colorbar
    p2(2).select()
    caxis([0 1])
    c1 = colorbar;
    c1.Label.String = 'Probability';
    c1.Location     = 'westoutside';
    set(gca, 'visible', 'off');
    
    %Layout adjustments
    p.de.margin  = 15;
    p2.de.margin = 25;
    
    %Set state variables
    stateVar                     = zeros(size(post,1),1);   %Unclassified
    stateVar(post(:,nanIdx)>0.9) = 1;                       %Nanog
    stateVar(post(:,gatIdx)>0.9) = 2;                       %Gata
    
    %Storage variable
    predict{tt}(:,1)  = resvec_another{6,tt}(:,1); %X values
    predict{tt}(:,2)  = resvec_another{6,tt}(:,2); %Y Values
    predict{tt}(:,3)  = resvec_another{6,tt}(:,3); %Image Field number
    predict{tt}(:,4)  = stateVar;                  %State variable
end

%% Distance plot
%Note iterate over fields is the corrent place to find

%Inputs:
radiusCutoff = 30; %Local radius uM or px based on image calib
tempFieldNanog = [];

%Iterate over treatments
for tt = outputs.limconscatID
    fieldNums            = unique(predict{tt}(:,3));
    fieldGlobalCellN{tt} = length(nonzeros(predict{tt}(:,4)));
    fieldGlobalNanog{tt} = sum(predict{tt}(:,4)==1)/fieldGlobalCellN{tt};
    compCounts           = [];
    
    %Iterate over fields @predict{tt}(:,3)
    for ff = fieldNums'
        tempVec = predict{tt}(predict{tt}(:,3)==ff, :);
        tDist   = eucDist(tempVec(:,1),tempVec(:,2));                      %Matrix of euclidean distances
        resvec_distMat{tt,ff}     = tDist;                                 %Storage
        tDist(tDist>radiusCutoff) = 0;                                     %set nonLocal to 0 distance
        
        %Iterate over cells in field
        for cc = 1:length(tDist)
            
            %Store values for cells in local radius
            nearbyLogical    = tDist(cc,:)>0;                              %Ignore cells with '0' distance
            nearbyStates     = tempVec(nearbyLogical, 4);                  %index the state variable @(:,4)
            nHoodLoc         = length(nonzeros(nearbyStates));             %nz ignores state 0 cells from calcs
            nearbyNanog      = sum(nearbyStates==1)/nHoodLoc;              %Fraction of Nanog neighbors (1-Gata strict)
            
            %Storage
            selfState        = tempVec(cc, 4);
            compCounts       = [compCounts; horzcat(nearbyNanog, selfState, nHoodLoc)];
            
            clear nearbyStates globalStates
        end
        clear tDistGlob tDist tempVec
    end
    
    %Storage [nearbyNanogFraction; selfState; n-Neighbours]
    resvec_compCounts{tt} = compCounts;
    resvec_compCountsG{tt}= compCountsG;
end

%% Neighborhood plots - Calc composition for N+ and G+ cells
%Note: this is always NANOG proportion for N+ and G+ cells. Variables with
% 'Nanog' or 'Gata' in their names refers to self identity, not nhood proportion

tt = 0;

%Iterate over active treatments
for tt2 = nonzeros(outputs.limconscatID)'
    tt = tt+1;
    
    %Index self N+ and G+ cells
    idxNanLog = resvec_compCounts{tt2}(:,2)==1;
    idxGatLog = resvec_compCounts{tt2}(:,2)==2;
    
    %Sort by fraction of Nanog neighbors for N+ and G+ cells
    selfNanog = sort(resvec_compCounts{tt2}(idxNanLog, 1), 'descend');
    selfGata  = sort(resvec_compCounts{tt2}(idxGatLog, 1), 'descend');
    
    %Delete NaNs - i.e. cells with no close neighbors
    selfNanog(isnan(selfNanog)) = [];
    selfGata(isnan(selfGata))   = [];
       
    %Storage for plotting
    selfNanogRes{tt} = selfNanog;
    selfGataRes{tt}  = selfGata;
    figLabs{tt}      = file.treatmentLabels{tt2};
end

%% Generate synthetic data
replacementType = 'Replacement';
for rep = 1:1000 %Number of iterations
    for tt = outputs.limconscatID
        
        %Synthetic data is same size as number of cells in treatment tt
        tempN        =  floor(fieldGlobalCellN{tt}*fieldGlobalNanog{tt});
        tempG        =  fieldGlobalCellN{tt} - tempN;
        
        if tempN>0 %Control conditions fuck with indexing
            
            %state var vector has same frac of N+ cells as fieldGlobalCellN
            allCellCloud =  ones(fieldGlobalCellN{tt}, 1);
            allCellCloud(1:tempN)   = 1; %N+ is state1
            allCellCloud(tempN:end) = 2; %G+ is state2
            
            %Randomize state variable vector
            allCellCloud = allCellCloud(randperm(length(allCellCloud)));
            
            %With or without replacement
            switch replacementType
                case 'Replacement'
                    %Pick n cells from allCellCloud (with replacement)
                    synthNHood = [];
                    for cc = 1:length(resvec_compCounts{tt})
                        tempnHood   = resvec_compCounts{tt}(cc,3);                   %nCells in nHood
                        randIdx     = randperm(length(allCellCloud), tempnHood);     %pick nCells from synth dat
                        synthSubset = allCellCloud(randIdx);
                        synthFrac   = sum(synthSubset==1)/length(synthSubset);       %Fraction of N+ cells
                        synthNHood  = [synthNHood; horzcat(synthFrac,resvec_compCounts{tt}(cc,2))];
                        clear tempnHood randIdx
                    end
                    
                    %[synthFrac; selfIdentity]
                    synthNanogNhood{tt} = synthNHood;
                    
                case 'NoReplacement'
                    %Pick n cells from allCellCloud (without replacement)
                    synthNHood = [];
                    for cc = 1:length(resvec_compCounts{tt})
                        tempnHood   = resvec_compCounts{tt}(cc,3);                   %nCells in nHood
                        randIdx     = randperm(length(allCellCloud), tempnHood);     %pick nCells from synth dat
                        synthSubset = allCellCloud(randIdx);
                        synthFrac   = sum(synthSubset==1)/length(synthSubset);       %Fraction of N+ cells
                        synthNHood  = [synthNHood; horzcat(synthFrac,resvec_compCounts{tt}(cc,2))];
                        
                        %delete chosen values from allCellCloud
                        allCellCloud(randIdx) = [];
                        clear tempnHood randIdx
                    end
                    
                    synthNanogNhood{tt} = synthNHood;
                    
            end
        end
    end
    
    for tt = outputs.limconscatID
        
        %Control conditions are empty by previous decress
        if ~isempty(synthNanogNhood{tt})
            
            %Index self N+ and G+ cells
            nanIdx = find(synthNanogNhood{tt}(:,2)==1);
            gatIdx = find(synthNanogNhood{tt}(:,2)==2);
            
            %N+ and G+ neighborhood in different variables
            synthSelfNanog = sort(synthNanogNhood{tt}(nanIdx,1), 'ascend');
            synthSelfGata  = sort(synthNanogNhood{tt}(gatIdx,1), 'ascend');
            
            %Delete NaNs
            synthSelfNanog(isnan(synthSelfNanog)) = [];
            synthSelfGata(isnan(synthSelfGata))   = [];
            
            %Storage
            synthResSelfNanog{tt}(:,rep)= synthSelfNanog;
            synthResSelfGata{tt}(:,rep) = synthSelfGata;

        end
    end
end

%Average all rnd repeats
synthResSelfNanog = cellfun(@(x) mean(x,2),synthResSelfNanog, 'UniformOutput', 0);
synthResSelfGata  = cellfun(@(x) mean(x,2), synthResSelfGata, 'UniformOutput', 0);


%% Neighborhood plots - Make figures
figure
pp  = panel();
pp.pack(3, length(outputs.limconscatID));
pp1 = pp(1);
pp2 = pp(2);

for nRow = 1:3
    
    switch nRow
        
        %nHood of N+ cells
        case 1
            for nCol = 1:length(outputs.limconscatID)
                
                pp(nRow, nCol).select();
                nItems  = length(selfNanogRes{nCol});
                nTr     = outputs.limconscatID(nCol);
                
                %Sequential plotting for multicolor bars
                hold on
                for ss = 1:nItems
                    colIdx = floor((selfNanogRes{nCol}(ss)*100)+1);     %101 idx for colGrad avoids indexing error here
                    barh(ss, selfNanogRes{nCol}(ss), ...
                        'BaseValue', fieldGlobalNanog{nTr}, ...
                        'BarWidth', 1, 'FaceColor', colGrad(colIdx,:),...
                        'EdgeColor', 'none')
                end
                
                plot(synthResSelfNanog{nTr},length(synthResSelfNanog{nTr}):-1:1, 'Color', [0 0 132]./255, 'LineWidth', 2)
                set(gca, 'XTick', [0;0.5;1], 'XTickLabel', {[]})
                %set(gca, 'XTick', 0:0.5:1, 'XTickLabel', {' Only \newline GATA6', 'Mixed', ' Only \newline NANOG' })
                ylabel('Nanog+ cell index')
                xlim([-0.1 1.1])
                title('Mean neighborhood: N+ cell')
                set(gca, 'FontSize', 12, 'TickDir', 'out')
                if nCol>1; ylabel([]); end
                
            end
            
            %nHood of G+ cells
        case 2
            for nCol = 1:length(outputs.limconscatID)
                
                pp(nRow, nCol).select();
                nItems  = length(selfGataRes{nCol});
                nTr     = outputs.limconscatID(nCol);

                %Sequential plotting for multicolor bars
                hold on
                for ss = 1:nItems
                    colIdx = floor((selfGataRes{nCol}(ss)*100)+1);     %101 idx for colGrad avoids indexing error here
                    barh(ss, selfGataRes{nCol}(ss), ...
                        'BaseValue', fieldGlobalNanog{nTr},...
                        'BarWidth', 1, 'FaceColor', colGrad(colIdx,:), ...
                        'EdgeColor', 'none')
                end
                
                plot(synthResSelfGata{nTr},length(synthResSelfGata{nTr}):-1:1, 'Color', [0 0 132]./255, 'LineWidth', 2)
                set(gca, 'XTick', [0;0.5;1], 'XTickLabel', {'0:100\newlineN+:G+', '50:50\newlineN+:G+', '100:0\newlineN+:G+' })
                ylabel('Gata+ cell index')
                xlim([-0.1 1.1])
                title('Mean neighborhood: G+ cell' )
                set(gca, 'FontSize', 12, 'TickDir', 'out')
                if nCol>1; ylabel([]); end
                
            end
            
        case 3
            for nCol = 1:length(outputs.limconscatID)
                pp(nRow, nCol).select();
                
                %Build stacked bar vector
                barvec(:,1) = [mean(selfNanogRes{nCol}); 1-mean(selfNanogRes{nCol})];
                barvec(:,2) = [mean(selfGataRes{nCol}); 1-mean(selfGataRes{nCol})];
                %barerr      = horzcat(std(selfNanogRes{nCol}), std(selfGataRes{nCol}));
                
                h = bar(barvec', 'Stacked', 'BarWidth', 0.8)
                hold on
                h(1).FaceColor = colmap(3,:)
                h(2).FaceColor = colmap(1,:)
                %errbar([1 2], barvec(1,:),barerr, barerr, 'k+-','linewidth', 1);
                set(gca, 'XTickLabel', {'Nanog \newlinecells', 'Gata \newlinecells'})
                set(gca, 'YTick', 0:0.5:1)
                ylabel('Avg. neighborhood composition')
                title(figLabs{nCol})
                set(gca, 'FontSize', 12, 'TickDir', 'out')
                if nCol>1; ylabel([]); end
                if nCol<2; legend({'Nanog neighbors', 'Gata neighbors'}); end
                
            end
            
            
    end
    
end

%% Plot number of neighbors:
for tt = outputs.limconscatID
    tempNanIdx   = resvec_compCounts{tt}(:,2)==1;
    tempGatIdx   = resvec_compCounts{tt}(:,2)==2;
    barVecNan{tt}= resvec_compCounts{tt}(tempNanIdx,3);
    barVecGat{tt}= resvec_compCounts{tt}(tempGatIdx,3);
    
end
figure
g = panel();
g.pack(2, length(outputs.limconscatID));

for nRow = 1:2
    switch nRow
        
        %N+ cell neighborhood
        case 1
            tt = 0;
            for nCol = 1:length(outputs.limconscatID)
                g(nRow, nCol).select();
                hist(barVecNan{outputs.limconscatID(nCol)})
                xlim([0 15])
                h = findobj(gca,'Type','patch');
                h.FaceColor = [0 0.5 0.5];
                h.EdgeColor = 'w';
                title(figLabs{nCol})
                if nCol==1; ylabel('Counts'); end
            end
            
            %G+ cell neighborhood
        case 2
            for nCol = 1:length(outputs.limconscatID)
                g(nRow, nCol).select();
                hist(barVecGat{outputs.limconscatID(nCol)})
                xlim([0 15])
                h = findobj(gca,'Type','patch');
                h.FaceColor = [180 55 87]./255;
                h.EdgeColor = 'w';
                title(figLabs{nCol})
                if nCol==1; ylabel('Counts'); end
                xlabel('Number of neighbors')
            end
    end
end




%% for me
% figure
% gg = panel();
% gg.pack(1,5);
%
% for jj = 1:5
% gg(1,jj).select();
% imagesc(resvec_distMat{jj+5, 1})
% title(figLabs{jj})
% end

%Legacy:
%barh(selfNanogRes{nCol}, 'BaseValue', 0.5, 'BarWidth', 1, 'FaceColor', colmap(3,:), 'EdgeColor', 'none')
%barh(selfGataRes{nCol}, 'BaseValue', 0.5, 'BarWidth', 1, 'FaceColor', colmap(1,:), 'EdgeColor', 'none')

%% local fx
function distMat = eucDist(x, y)
distMat = zeros(length(x));

for cc1 = 1:length(x)
    for cc2 = 1:length(x)  %Triangle = cc1, square = length(x)
        distMat(cc1, cc2) = sqrt( ((x(cc1) - x(cc2))^2) + ((y(cc1)-y(cc2))^2) );
    end
end


end
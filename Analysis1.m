%Dhruv Raina
%Nucelar Reader for ES Cell Images, Tester code v1 for image analysis.
%eventual intention is to run this in python
%dependencies: bioformats reader, Test12 (nuc extraction code), and
%dataplotter(data plotting and unpackaging code)

base_dir = '/Users/draina/Desktop/AnalysisImages/';
output_dir = '/Users/draina/Desktop/AnalysisOutputs/';

NoChannels = 3; %Can only be 2 or 3
StackSize = 1;
PrintIm = 0; %Yes or No
ResizeFlag = 1;
batchmode = 35;
%Background correction for Ch2 or Ch3
ch2_sub = 2;
ch3_sub = 16;


main_dirlist = dir(base_dir);
for ctr1 = 4:length(main_dirlist)
    temp_read = bfopen([base_dir main_dirlist(ctr1).name]);  %"Bioformats Reader"
    
    switch StackSize
        case (1)
            %how do you iterate over all images in a stack
            stack1 = temp_read{1,1};
            %stack2 = temp_read{2,1}; %etc.
    end
    
    switch NoChannels
        case(3)
            c1im = double(stack1{1,1}); %Nuc
            c2im = double(stack1{2,1}); %Sensor
            c3im = double(stack1{3,1}); %Antibody
    end
    
%     %Background Correction: Assuming flatfield correction
%     c2im = c2im-ch2_sub;
%     c3im = c3im-ch3_sub;
%     
%     %Bringing negatives to zero:
%     c2im(c2im<0) = 0;
%     c3im(c3im<0) = 0;
%         
    %Find Nucleus
    [nuclbl num] = Test12(c1im, batchmode);
    
    %Define Cytoplasm
    cell_lbl = imdilate(nuclbl, strel('disk', 40));
    nuclbldil = imdilate(nuclbl, strel('disk', 10));
    cell_lbl(nuclbldil>0) = 0;
    
    %Extract Data from raw images
    nucprops_ch1 = regionprops(nuclbl, c1im, 'Area', 'MeanIntensity');
    nucprops_ch2 = regionprops(nuclbl, c2im, 'Area', 'MeanIntensity');
    cytprops_ch2 = regionprops(cell_lbl, c2im, 'Area', 'MeanIntensity');
    nucprops_ch3 = regionprops(nuclbl, c3im, 'Area', 'MeanIntensity');
    cytprops_ch3 = regionprops(cell_lbl, c3im, 'Area', 'MeanIntensity');
        
    %Packaging Data:
    %File with Labels
    lbl_datavec{1,1} = 'Ch1_NucMInt';
    lbl_datavec{2,1} = 'Ch2_NucMInt';
    lbl_datavec{3,1} = 'Ch2_CytMInt';
    lbl_datavec{4,1} = 'Ch2_NucMInt';
    lbl_datavec{5,1} = 'Ch2_CytMInt';
    lbl_datavec{6,1} = 'ImgNum';
   % lbl_datavec{1, ctr1-2} = (main_dirlist(ctr1).name(1:end-4));
    lbl_list{ctr1-3} = main_dirlist(ctr1).name(1:end-4);
    
    %File with data: NO CALCULATIONS HERE>> Everything's been moved to the
    %file datanal.m
    datavec{2, ctr1-2} = [nucprops_ch1.MeanIntensity];
    datavec{3, ctr1-2} = [nucprops_ch2.MeanIntensity];
    datavec{4, ctr1-2} = [cytprops_ch2.MeanIntensity];
    datavec{5, ctr1-2} = [nucprops_ch3.MeanIntensity];
    datavec{6, ctr1-2} = [cytprops_ch3.MeanIntensity];
    datavec{7, ctr1-2} = repmat(ctr1-3,1, num);
    
    %Printing Mask Images:
    if PrintIm == 1
        figure, imagesc(nuclbl)
        set(gcf, 'visible', 'off');
        hgexport(gcf, [output_dir main_dirlist(ctr1).name(1:end-4) '_NucMask.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        close gcf;
        
        figure, imagesc(cell_lbl)
        set(gcf, 'visible', 'off');
        hgexport(gcf, [output_dir main_dirlist(ctr1).name(1:end-4) '_CytMask.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        close gcf;
        
        figure, imagesc(c1im)
        set(gcf, 'visible', 'off');
        hgexport(gcf, [output_dir main_dirlist(ctr1).name(1:end-4) '_OrigImg.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        close gcf;
        
        %Overlay Image
        imolay = nuclbl;
        imolay(nuclbl>0)=1;
        imolay = imolay.*c1im;
        figure, imagesc(imolay)
        set(gcf, 'visible', 'off');
        hgexport(gcf, [output_dir main_dirlist(ctr1).name(1:end-4) '_NucOverlay.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        close gcf;
    end
end
keyboard

%to extract the data, just do cell2mat(datavec) to get an easy-to-work-with
%matrix
dataplotter(datavec, lbl_list)
datanal(datavec, lbl_list, lbl_datavec)







%load needed files
dir_cell = {'/Users/draina/Desktop/2016_E2_7_tiff'}; %Add more entries as needed

%User Inputs:

%Define Excel Data:
% flagdft  = 1;  %0 - No Plot; 1 - Plot; 2 - Keyboard in fft
% flagpeaky= 1;  %As yet undefined
%% Main:
%Parse directory, find the right files:

for c1 = 1:length(dir_cell)
    base_dir = dir_cell{c1};                                               %no slash needed (write the slashchanger fn to handle errors here)
    dir_list = dir(base_dir);
    
    for c2 = 4:length(dir_list)   %stupid .dsstore file
        file.path = base_dir;
        file.name = dir_list(c2).name;
        keyboard
        im = tiffread2([file.path '/' file.name]);
        nuc_orig = mat2gray(im(1).data);
        cyt_orig = mat2gray(im(4).data);
        erk_orig = mat2gray(im(3).data);
        gfp_orig = mat2gray(im(2).data);
        
        [nucmask, nucnum, cytmask, cnum] = ecadseg(nuc_orig, cyt_orig);
        
        
        figure, imagesc(nucmask)
        figure, imagesc(cytmask)
        figure, imagesc(cyt_orig)
        %calculate overlaps of nuc and cyt
        
        
        %measurements:
        for c3 = 1:cnum
            
            
            
            
        end
    end
end
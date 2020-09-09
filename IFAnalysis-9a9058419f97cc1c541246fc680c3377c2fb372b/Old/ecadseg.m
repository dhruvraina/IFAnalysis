%E-Cad + Nuc segmenter
%im = tiffread2('/Users/draina/Desktop/2016_E2_7_NoPDO3_Pos001_TEST.tif');

function varargout = ecadseg(cyt, nuc)

% cyt = mat2gray(im(4).data);
% nuc = mat2gray(im(1).data);
% ch2 = mat2gray(im(2).data);
cyt_open = imopen(cyt, strel('disk', 3));
cyt_blur = mat2gray(imgaussfilt(cyt, 0.5));

% tt_flat = reshape(cyt_blur, size(cyt_blur,1)*size(cyt_blur,2), 1);
% [idx, c, sumcent] = kmeans(tt_flat,6, 'emptyaction', 'drop');
% if max(isnan(c))
%     [idx, ~, sumcent] = kmeans(tt_flat,3, 'emptyaction', 'drop');
%     z=1;
%     msgbox('Please check image',  'Whatevs' )
% end
% reshape_im = reshape(idx, size(cyt_blur));
% figure, imagesc(reshape_im)

%b2 = imopen(reshape_im, strel('disk', 1));

%For Nuc:
bw2 = edge(mat2gray(cyt), 'canny', [0.1 0.30]);
bw3 = imdilate(bw2, strel('disk', 2));
bw4 = bwmorph(bw3, 'thin', Inf);
figure, imagesc(bw4)
bw4 = ~bw4;
bw5 = cyt.*bw4;
[nmask nnum] = internuc(nuc, bw4);
lblbw=nmask;
lblbw(nmask>0)=1;

%For Cyt:
% hy = fspecial('sobel');
% hx = hy';
% Iy = imfilter(double(cyt), hy, 'replicate');
% Ix = imfilter(double(cyt), hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);
%
% D = bwdist(lblbw);
% DL = watershed(D);
% bgm = DL == 0;
%
% ww = imimposemin(gradmag, bgm | lblbw);
% w = watershed(ww);

figure, imagesc(bw4)
hline = imline(gca)

binline = hline.createMask();
bw4 = ~bw4;
new_bw4 = binline+bw4;
new_bw4(new_bw4>1)=1;
new_bw4 = ~new_bw4;
new_bw5 = bwareaopen(new_bw4, 50, 4);


[c_lbl c_num] = bwlabel(new_bw5,4); %have to bwlabel because regionprops looks for 8 connectivity not 4
cprop = regionprops(c_lbl, c_lbl, 'Area', 'MeanIntensity');

[filt_hv filt_hidx] = max([cprop.Area]);
c_lbl(c_lbl==filt_hidx) = 0;

cprop(filt_hidx).Area = 0; %Set background to zero

filt_l = ceil(mean([cprop(:).Area])/3);

imclear = ismember(c_lbl, find([cprop.Area] >= filt_l));


for ct1 = 1:c_num
    tmask = zeros(size(c_lbl));
    tmask(c_lbl==ct1) = 1;
    tmask = tmask.*imclear;
    tmask2 = tmask.*nmask;
    assoc(ct1).cytmask = ct1;
  assoc(ct1).nucmask = max(tmask2(:));
    if max(tmask2(:))==0
        c_lbl(c_lbl==ct1)=0;
    end
end


varargout{1} = nmask;
varargout{2} = nnum;
varargout{3} = c_lbl;
varargout{4} = c_num;



% hy = fspecial('sobel');
% hx = hy';
% Iy = imfilter(double(cyt_blur), hy, 'replicate');
% Ix = imfilter(double(cyt_blur), hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);

%Step 2: Watershed:



%eliminate all cells smaller than a certain relative value, then imposeMAX, use nuclei as localMIN and watershed.

%Pop up image, click on useless cells to exclude them
function varargout = getnucsize(nuc)

%User question:
dlg_prompt = {'Use saved strel size?', 'Save new strel?'};
dlg_title = 'Nucleus Size?';
dlg_num_lines = 1;
dlg_defans = {'0', '0'};
dlg_userAns = inputdlg(dlg_prompt,dlg_title,dlg_num_lines,dlg_defans);
dlg_userAns = cellfun(@(x) str2double(x), dlg_userAns);

%LoadFile for last used strel: is the user an idiot?:
if dlg_userAns(1)==1
    iniFile = fullfile(pwd, 'nucseg_ini.mat');
    if exist(iniFile, 'file')==2 %If it exists as a file
        segStrel = load(iniFile);
        if isfield(segStrel,'segStrel')
        else
            dlg_userAns(1) = 0;
        end
    else
        dlg_userAns(1) = 0;
    end
end

%LoadFile for last used strel: actually getting the size of nucleus:
if dlg_userAns(1)==0
figure, imagesc(nuc)
[xdat, ydat] = getline;
segStrel = ceil(sqrt((xdat(1)-xdat(2))^2+(ydat(1)-ydat(2))^2));
close gcf
end

if dlg_userAns(2)==1
    save('nucseg_ini.mat', 'segStrel');
end    


varargout{1} = segStrel

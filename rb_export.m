function files = rb_export(fig,basename)
[folder,~,ext]=fileparts(basename);
if ~isempty(ext),basename=basename(1:end-numel(ext));end
if isempty(folder),folder=pwd;end
if ~exist(folder,'dir'),mkdir(folder);end
set(fig,'Visible','on','Color','white','InvertHardcopy','off', ...
    'ToolBar','none','MenuBar','none','Units','centimeters');
sz=get(fig,'Position');sz=sz(3:4);paper=sz/2.54;
set(fig,'PaperUnits','inches','PaperPositionMode','manual', ...
    'PaperSize',paper,'PaperPosition',[0 0 paper]);
drawnow;
pdfTemporary=[tempname(folder) '.pdf'];
pngTemporary=[tempname(folder) '.png'];
cleanup=onCleanup(@() cleanup_files(pdfTemporary,pngTemporary));
files={};
fprintf('  PNG export at 600 dpi: %s.png\n',basename);
exportClock=tic;
try
    if exist('OCTAVE_VERSION','builtin')
        print(fig,pngTemporary,'-dpng','-r600');
    else
        exportgraphics(fig,pngTemporary,'Resolution',600,'BackgroundColor','white');
    end
    assert(valid_png(pngTemporary),'The PNG export is incomplete.');
    if isappdata(fig,'rb_check_coefficient_maps'),check_maps(pngTemporary);end
    destination=[basename '.png'];
    [ok,msg]=movefile(pngTemporary,destination,'f');assert(ok,'PNG save failed: %s',msg);
    files{end+1,1}=destination;
    fprintf('  PNG saved in %.3f s. Numerical results and PNG are available before PDF export.\n',toc(exportClock));
catch exception
    warning('rb:PNGExportFailed','PNG export failed: %s. The native figure remains open.',exception.message);
end
fprintf('  PDF export: %s.pdf\n',basename);
exportClock=tic;
try
    if exist('OCTAVE_VERSION','builtin')
        fprintf('  Octave print fallback; MATLAB image-PDF export is not exercised.\n');
        print(fig,pdfTemporary,'-dpdf','-r600');
    else
        fprintf('  MATLAB image-PDF export at 600 dpi, including text and curves.\n');
        exportgraphics(fig,pdfTemporary,'ContentType','image','Resolution',600,'BackgroundColor','white');
    end
    assert(valid_header(pdfTemporary,'%PDF-'),'The PDF export is incomplete.');
    destination=[basename '.pdf'];
    [ok,msg]=movefile(pdfTemporary,destination,'f');assert(ok,'PDF save failed: %s',msg);
    files{end+1,1}=destination;
    fprintf('  PDF saved in %.3f s.\n',toc(exportClock));
catch exception
    warning('rb:PDFExportFailed','PDF export failed: %s. Saved PNG and numerical data are retained; the native figure remains open.',exception.message);
end
set(fig,'Visible','on');drawnow;
end

function yes=valid_png(file)
yes=false;information=dir(file);
if isempty(information)||information.bytes<128,return;end
metadata=imfinfo(file);
yes=metadata.Width>0&&metadata.Height>0;
end

function check_maps(file)
pixels=imread(file);
assert(ndims(pixels)==3&&size(pixels,3)>=3,'The exported PNG must contain RGB pixels.');
rows=1:floor(.56*size(pixels,1));
halves={1:floor(size(pixels,2)/2),floor(size(pixels,2)/2)+1:size(pixels,2)};
for k=1:2
    sample=double(pixels(rows,halves{k},1:3));
    colored=max(sample,[],3)-min(sample,[],3)>32;
    assert(mean(colored(:))>.04,'Coefficient map %d is missing from the exported PNG.',k);
end
end

function yes=valid_header(file,expected)
yes=false;information=dir(file);
if isempty(information)||information.bytes<128,return;end
fid=fopen(file,'r');if fid<0,return;end
header=fread(fid,numel(expected),'*char').';fclose(fid);
yes=strcmp(header,expected);
end

function cleanup_files(varargin)
for k=1:numel(varargin)
    if exist(varargin{k},'file'),delete(varargin{k});end
end
end

function files=rb_plot(result,folder)
source='current in-memory reconstructions';
if nargin<2||isempty(folder),folder=fullfile(pwd,'figures');end
[results,studies,selected]=prepare_results(result);
R=results{selected};runs=R.runs;
cfg=struct;if isfield(R,'cfg'),cfg=R.cfg;end
noise=studies(selected).noiseRuns;mesh=studies(selected).meshRuns;
reconstructionPower=studies(selected).p;
if ~exist(folder,'dir'),mkdir(folder);end
seed=option(cfg,'plot_seed',20261012);
levels=sort(unique([runs(noise).delta]));
if isfield(R,'curveDeltas')&&~isempty(R.curveDeltas)
    requested=R.curveDeltas(:).';
else
    requested=levels(unique(round(linspace(1,numel(levels),min(3,numel(levels))))));
    requested=fliplr(requested);
end
representatives=zeros(size(requested));
for k=1:numel(requested)
    ix=noise(close_to([runs(noise).delta],requested(k)));
    assert(~isempty(ix),'A requested reconstruction noise level is not present.');
    representatives(k)=select_seed(runs,ix,seed);
end
ix=noise(close_to([runs(noise).delta],min(levels)));
reconstruction=select_seed(runs,ix,seed);
assert(all(isfinite([runs(noise).relative_q_error]))&&all(isfinite([runs(mesh).relative_q_error])), ...
    'Nonfinite study errors must be resolved before publication plotting.');
if R.d==1
    fig=new_figure([16.8 6.5]);
    ax1=new_axes(fig,[.066 .245 .270 .590]);
    ax2=new_axes(fig,[.442 .245 .212 .590]);
    ax3=new_axes(fig,[.778 .245 .210 .590]);
    draw_curves(ax1,R,runs,representatives,cfg);
    title(ax1,'(a) Coefficient recovery','FontWeight','normal','FontSize',10);
    [noiseStats,meshStats]=draw_studies(ax2,ax3,studies);
    title(ax2,math_display('(b) Errors against $\delta$'),'Interpreter',math_interpreter(),'FontWeight','normal','FontSize',10);
    title(ax3,math_display('(c) Errors against $h$'),'Interpreter',math_interpreter(),'FontWeight','normal','FontSize',10);
else
    fig=new_figure([16.8 9.6]);
    setappdata(fig,'rb_check_coefficient_maps',true);
    ax1=new_axes(fig,[.150 .56 .200 .350]);
    ax2=new_axes(fig,[.635 .56 .200 .350]);
    ax3=new_axes(fig,[.10 .130 .31 .270]);
    ax4=new_axes(fig,[.59 .130 .31 .270]);
    limits=map_limits(R,runs(reconstruction));
    draw_map(ax1,R,[],cfg,limits);
    title(ax1,{math_display('(a) True coefficient'),math_display(sprintf('$p=%g$',reconstructionPower))}, ...
        'Interpreter',math_interpreter(),'FontWeight','normal','FontSize',10);
    draw_map(ax2,R,runs(reconstruction),cfg,limits);
    title(ax2,{math_display('(b) Recovered'),math_display(sprintf('$p=%g, \\delta=%s$',reconstructionPower,scientific(runs(reconstruction).delta)))}, ...
        'Interpreter',math_interpreter(),'FontWeight','normal','FontSize',10);
    draw_colourbar(fig,[.885 .570 .012 .325],limits);
    [noiseStats,meshStats]=draw_studies(ax3,ax4,studies);
    title(ax3,math_display('(c) Errors against $\delta$'),'Interpreter',math_interpreter(),'FontWeight','normal','FontSize',10);
    title(ax4,math_display('(d) Errors against $h$'),'Interpreter',math_interpreter(),'FontWeight','normal','FontSize',10);
end
files=rb_export(fig,fullfile(folder,sprintf('ex%d_panels',R.d)));
manifest=struct('source',source,'powers',[studies.p],'reconstructionPower',reconstructionPower, ...
    'reconstructionResult',selected,'reconstructionRun',reconstruction,'studiesByPower',studies, ...
    'curveRuns',representatives,'noiseRuns',noise,'meshRuns',mesh, ...
    'noiseStatistics',noiseStats,'meshStatistics',meshStats, ...
    'selection','Prespecified noise level and seed; no error-based selection.', ...
    'fieldRepresentation','Exact radial P0 values; 1D stairs, 2D RGB pixel sampling with nearest-neighbour display when supported; no reconstruction is assigned outside D.', ...
    'mapResolution',option(cfg,'map_pixels',1000), ...
    'mapMask','Grey in the unfitted annulus D<r<=1; white outside the unit disk.', ...
    'displaySampling','Recovered pixels take the containing cell value without smoothing; display sampling is independent of physical error quadrature.', ...
    'exportRepresentation',export_representation(), ...
    'errorAggregation','Separate arithmetic mean of the five relative errors for each p; minima and maxima are retained in the saved statistics but not drawn.', ...
    'files',{files});
save(fullfile(folder,sprintf('ex%d_plot_manifest.mat',R.d)),'manifest','-v7');
end

function [results,studies,selected]=prepare_results(input)
if iscell(input)
    results=input(:).';
elseif isstruct(input)&&isscalar(input)
    results={input};
else
    error('rb_plot requires a result structure or a cell array of result structures.');
end
assert(~isempty(results),'At least one result structure is required.');
powers=zeros(1,numel(results));studyCells=cell(size(results));
for k=1:numel(results)
    R=results{k};
    assert(isstruct(R)&&isscalar(R)&&all(isfield(R,{'d','runs','truthRadius','truthQ'})), ...
        'Each plot input requires d, runs, truthRadius and truthQ.');
    assert(ismember(R.d,[1 2]),'Only dimensions 1 and 2 are supported.');
    cfg=struct;if isfield(R,'cfg'),cfg=R.cfg;end
    if numel(results)>1,assert(isfield(cfg,'p'),'Every power-comparison result needs cfg.p.');end
    powers(k)=option(cfg,'p',2);
    assert(ismember(powers(k),[1 2 3]),'The plotted powers must be 1, 2 or 3.');
    runs=R.runs;
    assert(~isempty(runs)&&all(isfield(runs,{'study','h','delta','seed','q','edges','relative_q_error'})), ...
        'Run records must contain the required plotting fields.');
    noise=find(strcmp({runs.study},'noise')&[runs.delta]>0);mesh=find(strcmp({runs.study},'mesh'));
    assert(~isempty(noise)&&~isempty(mesh),'Both studies are required for every power.');
    noiseh=option(cfg,'noise_h',min([runs(noise).h]));
    meshdelta=option(cfg,'mesh_delta',min([runs(mesh).delta]));
    noise=noise(close_to([runs(noise).h],noiseh));mesh=mesh(close_to([runs(mesh).delta],meshdelta));
    assert(~isempty(noise)&&~isempty(mesh),'Configured study parameters were not found.');
    ns=aggregate(runs,noise,'delta');ms=aggregate(runs,mesh,'h');
    seeds=sort(unique([runs(noise).seed]));
    if numel(results)>1
        assert(numel(seeds)==5&&all(ns.count==5)&&all(ms.count==5), ...
            'Each power comparison requires all five recorded directions at every level.');
    end
    check_directions(runs,noise,'delta',ns.x,seeds);
    check_directions(runs,mesh,'h',ms.x,seeds);
    studyCells{k}=struct('p',powers(k),'dimension',R.d,'noiseRuns',noise,'meshRuns',mesh, ...
        'noiseStatistics',ns,'meshStatistics',ms,'noise_h',noiseh,'mesh_delta',meshdelta,'seeds',seeds);
end
[powers,order]=sort(powers);results=results(order);studies=[studyCells{order}];
if numel(results)>1
    assert(isequal(powers,[1 2 3]),'The power comparison requires exactly p=1, p=2 and p=3.');
    for k=2:numel(studies)
        assert(studies(k).dimension==studies(1).dimension,'All compared results must have the same physical dimension.');
        assert(isequal(studies(k).noiseStatistics.x,studies(1).noiseStatistics.x)&& ...
            isequal(studies(k).meshStatistics.x,studies(1).meshStatistics.x), ...
            'All compared powers must use the same noise and mesh grids.');
        assert(studies(k).noise_h==studies(1).noise_h&&studies(k).mesh_delta==studies(1).mesh_delta&& ...
            isequal(studies(k).seeds,studies(1).seeds),'All compared powers must use the same study slices and directions.');
    end
    selected=find(powers==2);
else
    selected=1;
end
end

function check_directions(runs,idx,field,levels,seeds)
values=[runs(idx).(field)];
for k=1:numel(levels)
    selected=idx(close_to(values,levels(k)));
    assert(isequal(sort([runs(selected).seed]),seeds), ...
        'Study levels must contain the same recorded directions exactly once.');
end
end

function fig=new_figure(sz)
fig=figure('Visible','on','Color','white','Units','centimeters', ...
    'Position',[1 1 sz],'MenuBar','none','ToolBar','none');
end

function ax=new_axes(fig,position)
ax=axes('Parent',fig,'Units','normalized','Position',position);hold(ax,'on');
set(ax,'FontName','Times New Roman','FontSize',9.4,'LineWidth',.65,'Box','on', ...
    'TickDir','out','TickLength',[.007 .007],'Layer','top');
grid(ax,'on');
if isprop(ax,'GridColor'),set(ax,'GridColor',[.87 .87 .87]);end
if isprop(ax,'GridAlpha'),set(ax,'GridAlpha',1);end
if isprop(ax,'MinorGridAlpha'),set(ax,'MinorGridAlpha',0);end
if isprop(ax,'TickLabelInterpreter'),set(ax,'TickLabelInterpreter',math_interpreter());end
end

function draw_curves(ax,R,runs,idx,cfg)
[r,order]=sort(R.truthRadius(:));q=R.truthQ(order);
assert(r(1)==0,'The saved true coefficient must include r=0.');
handles=cell(1,numel(idx)+1);
handles{1}=plot(ax,[-flipud(r(2:end));r],[flipud(q(2:end));q],'-', ...
    'Color',[.12 .12 .12],'LineWidth',1.4);
limits=[q(:);0];
for k=1:numel(idx)
    run=runs(idx(k));e=run.edges(:);v=run.q(:);
    assert(numel(e)==numel(v)+1,'Radial cell edges and coefficient values do not align.');
    edges=[-flipud(e(2:end));e];values=[flipud(v);v];
    handles{k+1}=stairs(ax,edges,[values;values(end)],'Color',palette(k), ...
        'LineWidth',1.05,'LineStyle',style(k));
    limits=[limits;v];
end
radius=option(cfg,'curve_radius',max(runs(idx(1)).edges));
xlim(ax,[-radius radius]);
set(ax,'XTick',[-radius 0 radius]);
span=max(limits)-min(limits);if span==0,span=1;end
ylim(ax,[min(limits)-.035*span max(limits)+.75*span]);
positioned_labels(ax,'$x$','$q$',-.21,-.22);
labels=[{'$q$ (true)'},arrayfun(@(j) ['$\delta=' scientific(runs(j).delta) '$'],idx,'UniformOutput',false)];
inside_legend(ax,[handles{:}],labels,'north',2,9.2,sprintf('$p=%g$',option(cfg,'p',2)));
end

function [ns,ms]=draw_studies(noiseax,meshax,studies)
ns=[studies.noiseStatistics];ms=[studies.meshStatistics];
noiseHandles=draw_error(noiseax,ns);meshHandles=draw_error(meshax,ms);
labels=arrayfun(@(p)sprintf('$p=%g$',p),[studies.p],'UniformOutput',false);
inside_legend(noiseax,noiseHandles,labels,'northwest',1,9.2);
inside_legend(meshax,meshHandles,labels,'northwest',1,9.2);
figurepos=get(get(noiseax,'Parent'),'Position');
if figurepos(4)>8
    positioned_labels(noiseax,'$\delta$','$e_q$',-.40,-.20);
else
    positioned_labels(noiseax,'$\delta$','$e_q$',-.21,-.30);
end
values=ns(1).x;pick=unique(round(linspace(1,numel(values),min(4,numel(values)))));
set(noiseax,'XTick',values(pick),'XTickLabel',arrayfun(@(v) math_display(['$' scientific(v) '$']),values(pick),'UniformOutput',false));
if isprop(noiseax,'XTickLabelRotation'),set(noiseax,'XTickLabelRotation',0);end
if figurepos(4)>8
    positioned_labels(meshax,'$h$','$e_q$',-.40,-.20);
else
    positioned_labels(meshax,'$h$','$e_q$',-.30,-.29);
end
values=ms(1).x;set(meshax,'XTick',values,'XTickLabel',repmat({''},size(values)));
xl=get(meshax,'XLim');
for k=1:numel(values)
    xpos=(log(values(k))-log(xl(1)))/(log(xl(2))-log(xl(1)));
    text(meshax,xpos,-.075,math_display(sprintf('$%.4g$',values(k))), ...
        'Interpreter',math_interpreter(),'Units','normalized','Rotation',30,'FontName','Times New Roman','FontSize',9, ...
        'HorizontalAlignment','right','VerticalAlignment','top','Clipping','off');
end
end

function S=aggregate(runs,idx,field)
x=[runs(idx).(field)];error=[runs(idx).relative_q_error];
S=struct('x',sort(unique(x)),'mean',[],'minimum',[],'maximum',[],'count',[]);
for k=1:numel(S.x)
    values=error(close_to(x,S.x(k)));
    assert(all(isfinite(values)),'Nonfinite errors cannot be omitted from a study mean.');
    S.mean(k)=mean(values);S.minimum(k)=min(values);S.maximum(k)=max(values);S.count(k)=numel(values);
end
end

function handles=draw_error(ax,S)
set(ax,'XScale','log','YScale','log');handles=cell(1,numel(S));
styles={'-','--','-.'};markers={'o','s','^'};
for k=1:numel(S)
    assert(all(S(k).x>0)&&all(isfinite(S(k).mean))&&all(S(k).mean>0), ...
        'Error studies need finite positive abscissae and errors.');
    handles{k}=plot(ax,S(k).x,S(k).mean,'Color',palette(k),'LineWidth',1.35, ...
        'LineStyle',styles{1+mod(k-1,numel(styles))},'Marker',markers{1+mod(k-1,numel(markers))}, ...
        'MarkerSize',4.2,'MarkerFaceColor','white');
end
handles=[handles{:}];x=[S.x];means=[S.mean];
xlim(ax,[min(x)/1.35 max(x)*1.35]);
low=min(means);high=max(means);
factor=max(1.25,(high/low)^.12);ylim(ax,[low/factor high*factor]);
decades=10.^(ceil(log10(low/factor)):floor(log10(high*factor)));
if isempty(decades),decades=unique([low high]);end
set(ax,'YTick',decades,'YTickLabel',arrayfun(@(v) ...
    math_display(['$' scientific(v) '$']),decades,'UniformOutput',false));
end

function draw_map(ax,R,run,cfg,limits)
radius=option(cfg,'plot_radius',1);D=option(cfg,'D',.8);
np=option(cfg,'map_pixels',1000);
assert(isscalar(np)&&isfinite(np)&&np>=64&&np==round(np),'map_pixels must be an integer at least 64.');
x=((1:np)-.5)*(2*radius/np)-radius;
[X,Y]=meshgrid(x,x);r=hypot(X,Y);
inside=r<=D;disk=r<=1;
q=zeros(size(r));
if isempty(run)
    if isfield(cfg,'q_radius')&&isfield(cfg,'q_amplitude')
        q(inside)=cfg.q_amplitude*sj_bump(r(inside)/cfg.q_radius);
    else
        [tr,order]=sort(R.truthRadius(:));tq=R.truthQ(order);
        q(inside)=interp1(tr,tq,r(inside),'linear');
    end
else
    e=run.edges(:);v=run.q(:);
    assert(numel(e)==numel(v)+1&&e(1)==0&&all(diff(e)>0), ...
        'Coefficient cell edges and values must be aligned and strictly ordered.');
    assert(abs(e(end)-D)<=64*eps(max(1,D)),'The coefficient cells must cover D.');
    [~,index]=histc(r(inside),e);
    index=min(index,numel(v));
    assert(all(index>=1),'A display pixel is outside the coefficient cells.');
    q(inside)=v(index);
end
assert(all(isfinite(q(inside))),'Coefficient map values must be finite.');
map=colourmap();nc=size(map,1);
index=1+floor((nc-1)*(q-limits(1))/(limits(2)-limits(1)));
index=max(1,min(nc,index));
rgb=reshape(map(index(:),:),[np np 3]);
for k=1:3
    plane=rgb(:,:,k);plane(disk&~inside)=.92;plane(~disk)=1;rgb(:,:,k)=plane;
end
im=image(ax,[x(1) x(end)],[x(1) x(end)],rgb);
if isprop(im,'Interpolation'),set(im,'Interpolation','nearest');end
setappdata(im,'rb_radial_cell_sampling',true);
theta=linspace(0,2*pi,721);ct=cos(theta);st=sin(theta);
plot(ax,ct,st,'Color',[.4 .4 .4],'LineWidth',.5);
plot(ax,D*ct,D*st,'-','Color',[1 1 1],'LineWidth',1.15);
plot(ax,D*ct,D*st,'--','Color',[.25 .25 .25],'LineWidth',.70);
axis(ax,'equal');xlim(ax,[-radius radius]);ylim(ax,[-radius radius]);
set(ax,'YDir','normal','XTick',[-radius 0 radius],'YTick',[-radius 0 radius]);
grid(ax,'off');caxis(ax,limits);colormap(ax,colourmap());
positioned_labels(ax,'$x_1$','$x_2$',-.21,-.31);
end

function draw_colourbar(fig,position,limits)
ax=axes('Parent',fig,'Units','normalized','Position',position);
imagesc(ax,[0 1],limits,linspace(limits(1),limits(2),256)');
set(ax,'YDir','normal','XTick',[],'YAxisLocation','right','FontName','Times New Roman', ...
    'FontSize',9.2,'TickDir','out','TickLength',[.012 .012],'LineWidth',.5,'Box','on');
if isprop(ax,'TickLabelInterpreter'),set(ax,'TickLabelInterpreter',math_interpreter());end
ticks=linspace(limits(1),limits(2),5);
set(ax,'YTick',ticks,'YTickLabel',arrayfun(@(v)sprintf('%.2g',v),ticks,'UniformOutput',false));
ylim(ax,limits);caxis(ax,limits);colormap(ax,colourmap());
text(ax,.5,1.08,math_display('$q$'),'Interpreter',math_interpreter(), ...
    'Units','normalized','FontName','Times New Roman','FontSize',10.2, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','off');
end

function positioned_labels(ax,xtext,ytext,xoffset,yoffset)
xlabel(ax,'');ylabel(ax,'');
text(ax,.5,xoffset,math_display(xtext),'Units','normalized','Interpreter',math_interpreter(), ...
    'FontName','Times New Roman','FontSize',10.2,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','Clipping','off');
text(ax,yoffset,.5,math_display(ytext),'Units','normalized','Interpreter',math_interpreter(), ...
    'FontName','Times New Roman','FontSize',10.2,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','Clipping','off','Rotation',90);
end

function inside_legend(ax,handles,labels,location,columns,fontsize,heading)
if nargin<7,heading='';end
labels=cellfun(@math_display,labels,'UniformOutput',false);
fig=ancestor(ax,'figure');sz=get(fig,'Position');points=sz(3:4)*72/2.54;
pos=get(ax,'Position');rows=ceil(numel(labels)/columns);
headingInCell=~isempty(heading)&&numel(labels)<rows*columns;
headingRows=~isempty(heading)&&~headingInCell;
height=((rows+headingRows)*(fontsize+3)+4)/points(2);
if columns>1,width=.94*pos(3);else,width=55/points(1);end
left=pos(1)+.025*pos(3);
if strcmp(location,'north'),left=pos(1)+(pos(3)-width)/2;end
bottom=pos(2)+pos(4)-height-.025*pos(4);
lg=axes('Parent',fig,'Units','normalized','Position',[left bottom width height], ...
    'XLim',[0 1],'YLim',[0 1],'Box','on','Color','white','LineWidth',.5, ...
    'XTick',[],'YTick',[],'XColor',[.35 .35 .35],'YColor',[.35 .35 .35]);
hold(lg,'on');
if headingInCell
    text(lg,(columns-.5)/columns,.5/rows,math_display(heading),'Interpreter',math_interpreter(), ...
        'FontName','Times New Roman','FontSize',fontsize,'HorizontalAlignment','center','VerticalAlignment','middle');
end
if headingRows
    text(lg,.5,1-.5/(rows+headingRows),math_display(heading),'Interpreter',math_interpreter(), ...
        'FontName','Times New Roman','FontSize',fontsize,'HorizontalAlignment','center','VerticalAlignment','middle');
end
for k=1:numel(labels)
    row=floor((k-1)/columns);column=mod(k-1,columns);
    x=column/columns;y=1-(row+headingRows+.5)/(rows+headingRows);
    token=[3 16]/(width*points(1));labelx=20/(width*points(1));
    if strcmp(get(handles(k),'Type'),'patch')
        patch(lg,x+token([1 2 2 1]),y+[-.15 -.15 .15 .15]/(rows+headingRows), ...
            get(handles(k),'FaceColor'),'EdgeColor','none');
    else
        line(lg,x+token,[y y],'Color',get(handles(k),'Color'), ...
            'LineStyle',get(handles(k),'LineStyle'),'LineWidth',1.2);
        if isprop(handles(k),'Marker')&&~strcmp(get(handles(k),'Marker'),'none')
            line(lg,x+mean(token),y,'Color',get(handles(k),'Color'),'LineStyle','none', ...
                'Marker',get(handles(k),'Marker'),'MarkerSize',4.2,'MarkerFaceColor','white');
        end
    end
    text(lg,x+labelx,y,labels{k},'Interpreter',math_interpreter(), ...
        'FontName','Times New Roman','FontSize',fontsize, ...
        'HorizontalAlignment','left','VerticalAlignment','middle');
end
end

function value=math_interpreter()
value='latex';if exist('OCTAVE_VERSION','builtin'),value='tex';end
end

function value=math_display(value)
if exist('OCTAVE_VERSION','builtin')
    value=strrep(value,'$','');
end
end

function limits=map_limits(R,run)
values=[0;R.truthQ(:);run.q(:)];assert(all(isfinite(values)),'Nonfinite coefficient values cannot be plotted.');
limits=[min(values) max(values)];
if limits(1)==limits(2),limits=limits+[-1 1]*max(1,abs(limits(1)))*.01;end
end

function map=colourmap()
anchors=[.208 .166 .529;.118 .348 .725;0 .506 .744;0 .619 .690; ...
    .230 .696 .557;.590 .750 .351;.902 .765 .178;.976 .983 .054];
map=interp1(linspace(0,1,size(anchors,1)),anchors,linspace(0,1,256));
end

function c=palette(k)
colors=[0 114 178;213 94 0;0 158 115;204 121 167;230 159 0]/255;
c=colors(1+mod(k-1,size(colors,1)),:);
end

function s=style(k)
styles={'--','-.',':','-'};s=styles{1+mod(k-1,numel(styles))};
end

function value=option(cfg,key,default)
if isfield(cfg,key)&&~isempty(cfg.(key)),value=cfg.(key);else,value=default;end
end

function mask=close_to(values,target)
mask=abs(values-target)<=64*eps(max(abs(target),realmin));
end

function index=select_seed(runs,indices,seed)
match=indices([runs(indices).seed]==seed);
if isempty(match),[~,k]=min([runs(indices).seed]);index=indices(k);else,index=match(1);end
end

function s=scientific(x)
e=floor(log10(x));m=x/10^e;
if abs(m-1)<1e-10,s=sprintf('10^{%d}',e);else,s=sprintf('%.2g\\times10^{%d}',m,e);end
end

function description=export_representation()
if exist('OCTAVE_VERSION','builtin')
    description='PNG first at 600 dpi; Octave print PDF fallback, not MATLAB image-PDF export.';
else
    description='PNG first; PDF contains the full figure as a 600 dpi raster, including text and curves.';
end
end

function [x_mc,y_mc,x_mc0,y_mc0,S]=prepare_coastline(S)
if ~isempty(S.LDBcoastline)&~isfield(S,'x_mc')
    xy_mc=load(S.LDBcoastline);
    if isempty(S.XYoffset)
        S.XYoffset = [floor(min(xy_mc(:,1))/1000)*1000 , floor(min(xy_mc(:,2))/1000)*1000];
    end
    x_mc=xy_mc(:,1)' - S.XYoffset(1);   %x_mc=xy_mc(end:-1:1,1)';   % SHIFT COASLTINE
    y_mc=xy_mc(:,2)' - S.XYoffset(2);   %y_mc=xy_mc(end:-1:1,2)';   % SHIFT COASLTINE
    x_mc0=x_mc;
    y_mc0=y_mc;
    
elseif ~isfield(S,'x_mc')
    figure(99);clf;
    plot([S.xlimits(1) S.xlimits(2) S.xlimits(2) S.xlimits(1) S.xlimits(1)], ...
        [S.ylimits(1) S.ylimits(1) S.ylimits(2) S.ylimits(2) S.ylimits(1)],'k:');
    axis equal;
    xl=xlim;yl=ylim;
    xlabel('Easting [m]');
    ylabel('Northing [m]');
    htxt=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add coastline (LMB); Next segment (RMB); Exit (q)');set(htxt,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.1 0.6]);
    [x_mc,y_mc]=select_multi_polygon('r');
    set(htxt,'Visible','off');
    x_mc0=x_mc;
    y_mc0=y_mc;
    
else
    x_mc=S.x_mc;
    y_mc=S.y_mc;
    x_mc0=S.x_mc0;
    y_mc0=S.y_mc0;
    
end
if isempty(S.xlimits)
    S.xlimits = [min(x_mc),max(x_mc)];
else
    S.xlimits = S.xlimits - S.XYoffset(1);
end
if isempty(S.ylimits)
    S.ylimits = [min(y_mc),max(y_mc)];
else
    S.ylimits = S.ylimits - S.XYoffset(2);
end
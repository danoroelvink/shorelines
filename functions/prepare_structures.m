function [x_hard,y_hard]=prepare_structures(S)
if S.struct
    if isfield(S,'x_hard')
        x_hard=S.x_hard;
        y_hard=S.y_hard;
    elseif ~isempty(S.LDBstructures)
        xy_hard=load(S.LDBstructures);
        x_hard=xy_hard(:,1)'-S.XYoffset(1);
        y_hard=xy_hard(:,2)'-S.XYoffset(2);
        S.x_hard=x_hard;
        S.y_hard=y_hard;
    else
        figure(99);
        axis equal;
        xl=xlim;yl=ylim;
        htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add structure (LMB); Next structure (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
        [x_hard,y_hard]=select_multi_polygon('k');
        pset=set(htxt2,'Visible','off');
        S.x_hard=x_hard;
        S.y_hard=y_hard;
    end
else
    x_hard=[];
    y_hard=[];
end
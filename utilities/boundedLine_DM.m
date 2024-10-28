function [hl,hp] = boundedLine_DM(x,y,color)
%Simple adaptation of boundedline function from:
%https://github.com/kakearney/boundedline-pkg
%Ensure that x is a row vector, and y has size [3 x length(x)]
    
if numel(x) == 1    %special case
    hp = plot([x x],y([1 3]),'-','Linewidth',1.5,'Color',[color 0.2]);
    hl = plot(x,y(2,:),'.','Color',color,'MarkerSize',10);
else                %default behavior
    x2 = [x fliplr(x)];
    y2 = [y(1,:) fliplr(y(3,:))];
    hp = patch(x2,y2,color,'facealpha', 0.2, 'edgecolor', 'none');
    hl = plot(x,y(2,:),'-','Color',color,'Linewidth',1.5);
end

end %[EoF]
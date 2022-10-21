% Function turns the figure background white, increases font sizes, and
% allows mathematical formulae be inputted as figure titles. 
% (Makes figures appropriate for a Latex report)
function plot_latex(plt, xtitle, ytitle,ztitle, tl ,str)



%set the plot legend 
if isempty(str) == false 
    lgd = legend(str, 'Interpreter','latex','Orientation', 'Horizontal'); 
    lgd.FontSize = 15;
end 


% set x label, y label and title of the graph, as interpreted using latex 
xlabel(xtitle,'fontsize',15,'fontweight','bold', 'Interpreter','latex');
ylabel(ytitle,'fontsize',15,'fontweight','bold','Interpreter','latex');
zlabel(ztitle,'fontsize',15,'fontweight','bold','Interpreter','latex');
title(tl,'fontsize',15,'fontweight','bold','Interpreter','latex');

% sets the fontsize of the values on the axes.
set(gca,'fontsize',15)

% the default matlab plot background is grey, this turns it to white 
set(gcf,'color','w');


% increase the width of lines plotted on the graph
% set(plt , 'linewidth',2);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

box on
grid minor

%fullscreen 
%figure('units','normalized','outerposition',[0 0 6 1])
fig=gcf;
fig.Position(3:4)=[550,400];
end 

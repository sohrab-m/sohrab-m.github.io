function [y1_group,y2_group]=gtfill(x,y1,y2,n,c1,c2)
%This function fills the areas between two curves that intersect at multple
%points. The area where curve 1 > curve 2 is filled in one color, and the
%area where curve 2 > curve 1 is filled in another color.
%
%Handles for each grouped area set are returned
%
%x  = x values, length(x) must equal length(y1) && length(y2)
%y1 = y values of first curve
%y2 = y values of second curve
%n  = number of points to interpolate, n should be >> length(x) 
%
%Will Hobbs, January 2012
%Inspired by jbfill.m by John A. Bockstege, November 2006
%
%EXAMPLE:
% x=-pi:1:pi;
% y1=sin(x);
% y2=cos(x);
% n=200;
% figure
% hold on
% [h1,h2]=gtfill(x,y1,y2,n);
% set(get(h1,'Children'),'FaceColor',[1 .5 .2])
% plot(x,y1,'g','LineWidth',4)
% plot(x,y2,'m','LineWidth',4)
% legend('y1 > y2','y2 > y1','y1','y2')

if nargin<5;c1='r';c2='b';end %default colors are red, blue

%interpolate for fill: http://blogs.mathworks.com/loren/2008/08/25/piecewise-linear-interpolation/
x_interp=linspace(x(1),x(end),n);
y1_interp=interp1(x,y1,x_interp);
y2_interp=interp1(x,y2,x_interp);

x_intercepts=find(diff([0,sign(y1_interp-y2_interp)])~=0); %estimate intersections of y1_interp and y2_interp

if max(x_intercepts)<n %if the last intersection is not the last point
    x_intercepts(end+1)=n; %add end boundary
end

y1_greater_start=find(diff([0,sign(y1_interp(1:end-1)-y2_interp(1:end-1))])>0); %points where y1 crosses above y2, ignoring last point
y1_greater_end=x_intercepts(find(ismember(x_intercepts,y1_greater_start))+1); %points where y1 crosses below y2
y1_greater_bounds=[y1_greater_start',y1_greater_end'];

y2_greater_start=find(diff([0,sign(y2_interp(1:end-1)-y1_interp(1:end-1))])>0); %points where y2 crosses above y1, ignoring last point
y2_greater_end=x_intercepts(find(ismember(x_intercepts,y2_greater_start))+1); %points where y2 crosses below y1
y2_greater_bounds=[y2_greater_start',y2_greater_end'];

hold on;

for i=1:size(y1_greater_bounds,1)
    y1_fill_handle(i)=fill([x_interp(y1_greater_bounds(i,1):y1_greater_bounds(i,2)),fliplr(x_interp(y1_greater_bounds(i,1):y1_greater_bounds(i,2)))],...
        [y1_interp(y1_greater_bounds(i,1):y1_greater_bounds(i,2)),fliplr(y2_interp(y1_greater_bounds(i,1):y1_greater_bounds(i,2)))],c1);
end

for j=1:size(y2_greater_bounds,1)
    y2_fill_handle(j)=fill([x_interp(y2_greater_bounds(j,1):y2_greater_bounds(j,2)),fliplr(x_interp(y2_greater_bounds(j,1):y2_greater_bounds(j,2)))],...
        [y2_interp(y2_greater_bounds(j,1):y2_greater_bounds(j,2)),fliplr(y1_interp(y2_greater_bounds(j,1):y2_greater_bounds(j,2)))],c2);
end

%http://www.mathworks.com/help/techdoc/creating_plots/braliom.html
if exist('y2_fill_handle') && exist('y1_fill_handle')
    y1_group=hggroup; %group handles
    set(y1_fill_handle(:),'Parent',y1_group)
    y2_group=hggroup; %group handles
    set(y2_fill_handle(:),'Parent',y2_group)
    % Include these hggroups in the legend:
    set(get(get(y1_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    set(get(get(y2_group,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
end
hold off
end
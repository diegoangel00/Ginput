%% Template for worksheet set answers

%% Answers to in-class worksheet ##
Name = 'Angel_Garcia'; % << your name here...
WS = 10;              % << number of problem set or work sheet...

%% THE LINES BELOW ARE TO HELP ME GRADE
%  PLEASE COPY AND PASTE THIS (WITH ALL THE "%" SIGNS) 
%  INTO YOUR SCRIPT BELOW YOUR NAME!
%  save('wsid','Name','WS');
%  clear; close all; clc;
%  load('wsid.mat'); delete('wsid.mat')
%  unqdir = sprintf('%s_ws%02.0f',Name,WS);
%  mkdir(unqdir)
%  cd(unqdir)
%% LAST LINE YOU SHOULD COPY... BELOW HERE IS ALL YOU :)


%% Citations and collaborations
% On this worksheet set I collaborated with the following students:

% And I used the following websites:
%   http://www.matlabtips.com
%   Matlab Help Documentation
%% Question A
clearvars;
clear variables;
disp('==================================================')
disp(' >>>> Question A')
% part i
disp(' >> Part i')

interactive
%% Question B
disp('==================================================')
disp(' >>>> Question B')
% part i
disp(' >> Part i')
x = [0:pi/90:2*pi]';
y = x';
Z = sin(x.*y)./(1+x.*y);
[X,Y] = meshgrid(x,y);
fg = figure(2); clf
subplot(2,1,1)
mesh(X,Y,Z)
view(gca,[45 65])
zlabel('N','fontsize',19,'fontweight','bold');
ylabel('y','fontsize',19,'fontweight','bold');
xlabel('x','fontsize',19,'fontweight','bold');
fg.Position = [50 50 450 720];

% part ii
disp(' >> Part ii')

subplot(2,1,2)
contourf(X,Y,Z)
colorbar
view(gca,[0 90])
zlabel('N','fontsize',19,'fontweight','bold');
ylabel('y','fontsize',19,'fontweight','bold');
xlabel('x','fontsize',19,'fontweight','bold');
ylim([0 6])
%% Question C
disp('==================================================')
disp(' >>>> Question C')
data = load('WS10_data.mat');
xyear = data.population.year;
pops = data.population.pop;
co2 = mean(table2array(data.co2(:,2:4)),2);
c = 1;
for i = 1:height(xyear)
   if xyear(i) > 0
       x(c) = xyear(i);
       pop(c) = pops(i);
       c = c + 1;
   end
end
x = x';
pop = pop';
%plots to both axis
figure(3); clf
[ax, h1, h2] = plotyy(x,pop,x,co2);
set(ax(2),'ylim',[265 415],'xlim',[0 2010]);
set(ax(1),'ylim',[0 10*10^9],'xlim',[0 2010]);
ylabel('y','fontsize',19,'fontweight','bold');
ylabel(ax(1),'World Population'); % left y-axis 
ylabel(ax(2),'Mean global atm CO2(ppm)'); % right y-axis
title('Comparison of Population and CO2 levels, y0-2015');

%%
C = 1:1000;
F = 3:1003';









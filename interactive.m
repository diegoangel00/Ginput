clear all 
close all
clc

javaaddpath('IRIS-WS-2.0.18.jar'); %  add irisFetch Java
% Fetch data form IRIS
%   args: Net, Sta, Loc, Cha, Starttime, Endtime [,quality][,includePZ][,verbosity]
mytrace=irisFetch.Traces('IU','ANMO','00','BHZ','2010-02-27 06:30:00','2010-02-27 07:30:00');
sampletimes=linspace(mytrace.startTime,mytrace.endTime,mytrace.sampleCount);

%% Make figure
figure(1)
% clear figure
clf
% set figure position on screen
set(gcf,'pos',[1 1 1350 370])
hold on % allow >1 line

%% plot seismic data
htr = plot(sampletimes,mytrace.data,'linewidth',1.3);

% make plot pretty
ylim(max(abs(mytrace.data))*1.2*[-1 1])
datetick;
set(gca,'fontsize',14,'linewidth',1.5,'box','on')
title('USE CURSOR TO SELECT THE FIRST SEISMIC WAVE ARRIVAL','fontsize',22,'fontweight','bold')

%% pick P wave arrival time
okPtime = 0;
while okPtime ~= 1
    % let user select P arrival time
    figure(1)
    x = ginput(1);
    tP = x(1);
    % display time
    fprintf('Selected wave arrival time is %s\n',datestr(tP))

    % check in with user - option to re-pick if not ok
    fprintf('Compare the printed time above ^ to the x-axis on your plot.\n')
    yn = input('Does the selected time look right (yes/no)?    ','s');

    if strcmp(yn,'n') || strcmp(yn,'no')
        fprintf('Ok, select new time on plot\n');
    elseif strcmp(yn,'y') || strcmp(yn,'yes')
        fprintf('Great, moving on...\n');
        okPtime = 1;
    end  
end

%% dialog box to customise plot
fprintf('Let''s customize some plot aspects\n');

% variable to decide whether to stay in the while loop
okPlot = 0; 
% dialog box defaults
lw = '3'; cl = 'b'; yl = 'Range'; ysize = '20'; ph = 'P';
% stay in while loop until user sets okPlot to 1 to move on
while okPlot ~= 1
    % set some defaults - not important
    set(0, 'DefaultTextFontWeight', 'bold','DefaultTextFontSize', 16)
    % set up dialog box prompt - cell array
    prompt = {'Enter linewidth','Enter color for trace (r/g/b)','Enter ylabel','Enter yaxis Font Size','Enter seismic phase to highlight (P/S)'};
    % make dialog box using prompt for all the fields and the defaults as
    % set above
    pltopt = inputdlg(prompt,'Plot options',1,{lw,cl,yl,ysize,ph});
    % extract each individual value from the elements of the cell returned
    % by inputdlg
    lw = pltopt{1};
    cl = pltopt{2};
    yl = pltopt{3};
    ysize = pltopt{4}
    ph = pltopt{5};
    
    % hard-coded predicted phase arrival time
    t_ar_P = 7.341962820819386e+05;
    t_ar_S = 7.341962890164197e+05;
    % decide which arrival to plot, based on user input "ph" from above
    if strcmp(ph,'P')
            t_ar = t_ar_P;
    elseif strcmp(ph,'S')
            t_ar = t_ar_S;
    else
        error('You did not select either P or S phase! How silly...\n')
    end
    % plot selected phase on figure
    hp(1) = plot(t_ar*[1 1],max(abs(mytrace.data))*0.2*[-1 1],':m','linewidth',3);
    hp(2) = text(t_ar,max(abs(mytrace.data))*0.23,ph,'verticalalignment','bottom');
    
    % modify the plot properties as per the other user inputs
    set(htr,'linewidth',str2num(lw),'color',cl)
    ylabel(yl,'fontweight','bold','Fontsize',str2num(ysize))
    
    figure(1);

    % check with user if the plot looks ok          
    yn = input('Do you want to change any of the things you just modified? (yes/no)?    ','s');
    if strcmp(yn,'y') || strcmp(yn,'yes')
        fprintf('Ok, select plot values\n');
        delete(hp);
    elseif strcmp(yn,'n') || strcmp(yn,'no')
        fprintf('Great, moving on...\n');
        fprintf('Look back at the figure window and follow the title prompt\n');
        okPlot = 1;
    end  
end

%% select maximum amplitude


% give instruction in figure title
title('USE CURSOR TO SELECT MAXIMUM AMPLITUDE ON TRACE','fontsize',22,'fontweight','bold')

% first ginput input
xma = ginput(1);
maxamp = xma(2);

% add text with value of max amp, just below the plotted line
text(sampletimes(end/3),maxamp*0.95,sprintf('Max amp = %.2e counts',maxamp),...
    'fontsize',17,'verticalalignment','top','color','m','interpreter','latex')

title('USE CURSOR TO SELECT MINIMUM ON TRACE','fontsize',22,'fontweight','bold')

% Second ginput prompt
xm = ginput(1);
minamp = xm(2);

plot([sampletimes(1),sampletimes(end)],minamp*[1 1],'--m','linewidth',1.5)

text(sampletimes(end/3),minamp*0.95,sprintf('Min amp = %.2e counts',minamp),...
    'fontsize',17,'verticalalignment','bottom','color','m','interpreter','latex')

%% finish up
title('NOW YOU ARE A SEISMOLOGIST','fontsize',22,'fontweight','bold')
drawnow

str = input('Would you like to save the data this seismic data occured? \n','s');
if strcmpi(str,'yes') || strcmpi(str,'y')
    save('mytrace.mat','mytrace')
    disp('You have saved "mytrace" Have a great day!')
else
    disp('Have a great day!')
end
function Dimless_unsteadyT  
    % set parameters:
    T0=20;
    T1=100;
    alpha=5e-5;
    b=2;
    
    etaR=linspace(0,1,50);
    etaL=linspace(-1,0,50);
    eta=unique(cat(2,etaL,etaR)); % concatenates arrays along dim=2 (rows) 
                                  % and is passed in unique to remove 
                                  % duplicate zero element in array
    tau=linspace(0,1,400);
    xN=length(etaR);
    
    theta=pdepe(0,@pdefun,@pdeic,@pdebc,eta,tau); % theta(t,x)
    theta(1,:)=0; % enforces IC at tau=0
    theta2=-theta(:,xN:end)+1; % flips the y-coord values about x-axis to 
                               % resemble plots in Bird et al., 
                               % i.e. theta2=T-T0/T1-T0
    T=theta2*(T1-T0)+T0;
    x=etaR*b;
    t=tau*b^2/alpha;
    
    iterationCount=0;
    tauPlot=[0.01 0.04 0.1 0.2 0.4 0.6 1.0];
    
    % the array ind stores the index locations of tauPlot values
    ind=zeros(1,length(tauPlot));
    for i=1:length(tauPlot)
        [~,ind(i)]=closest(tau,tauPlot(i)); % function that finds the index 
                                            % whose value closest matching 
                                            % in tauPlot
    end
    
    % limits the axis values for the plot
    xlim([0 1]);
    ylim([0 1]);
    
    % calls the function for the interactive slider plot
    createGUI;
    
    % defines the model to solve required by pdepe function
    function [c,f,s] = pdefun(eta,tau,theta,DtauDeta)
        c = 1;
        f = DtauDeta;
        s = 0;
    end

    % defines the ICs for the model defined in pdefun
    function theta0 = pdeic(eta)
        theta0 = 1;      
    end

    % defines the BCs for the model defined in pdefun
    function [pl,ql,pr,qr] = pdebc(etal,thetal,etar,thetar,tau)
        pl = thetal;
        ql = 0;
        
        pr = thetar;
        qr = 0;
    end
    
    function createGUI
        % This function takes the results from pdepe function and contains
        % all the code relevant to creating the plots, the interactive
        % sliders, and all related aspects to the UI.

        % creates new (empty figure) window 
        fig1 = figure('ResizeFcn',@figPosition,'Visible','off'); 
        sVal = 1;   % slider index, used to initialize value of indices 
                    % which represent the index value of the dependant 
                    % variable array, eg. c_g(i,j)=f(t,x), (i or j)=sindex 
                    % depending on the slider wishing to display.
        sVal2 = 1;

        theta_tSlide = @(sVal)theta2(round(sVal),:);
        T_tSlide = @(sVal2)T(round(sVal2),:);
        
        ax1=axes('Parent',fig1,'Units','pixel');
        box on
        hold on
        ploth1 = plot(etaR,theta_tSlide(sVal));
        hold off
        xlabel('\eta =x/b');
        ylabel('\theta =(T-T_0)/(T_1-T_0)');
        ylim([0,1]);
        xlim([0,1]);  
        
        ax2=axes('Parent',fig1,'Units','pixel');
        box on
        hold on
        ploth2 = plot(x,T_tSlide(sVal2));
        hold off
        xlabel('x (m)');
        ylabel('T (C)');
        ylim([0,T1]);
        xlim([0,b]);
        
        ax3=axes('Parent',fig1,'Units','pixel','Visible','off');
        
        slider1 = uicontrol(fig1,'Style','slider',...
            'value',sVal,...
            'min',1,'max',length(tau));
        txt1 = uicontrol(fig1,'Style','text',...
            'String',sprintf('tau = %0.1f',tau(round(sVal))));
        hLstn1 = addlistener(slider1,'ContinuousValueChange',@updateplot1);
        
        slider2 = uicontrol(fig1,'Style','slider',...
            'value',sVal2,...
            'min',1,'max',length(t));
        txt2 = uicontrol(fig1,'Style','text',...
            'String',sprintf('t = %0.1f s',t(round(sVal2))));
        hLstn2 = addlistener(slider2,'ContinuousValueChange',@updateplot2);
        
        txt3 = uicontrol(fig1,'Style','text',...
            'String','T0 (C):');
        edit3 = uicontrol(fig1,'Style','edit','Callback',@editCallback3);
        
        function figPosition(varargin)
            figWidth=fig1.Position(3);
            figHeight=fig1.Position(4);
            VertSeg=3;
            
            ax1.OuterPosition(1)=0;
            ax1.OuterPosition(3)=figWidth/VertSeg;
            
            ax2.OuterPosition(1)=ax1.OuterPosition(3)+ax1.OuterPosition(1);
            ax2.OuterPosition(3)=figWidth/VertSeg;
            
            ax3.OuterPosition(1)=ax2.OuterPosition(3)+ax2.OuterPosition(1);
            ax3.OuterPosition(3)=figWidth/VertSeg;
            
            ax1.Position(3)=ax2.OuterPosition(1)-ax1.Position(1);
            ax2.Position(3)=ax3.OuterPosition(1)-ax2.Position(1);

            slider1.Position(1)=ax1.Position(1);
            slider1.Position(3)=ax1.Position(3);
            txt1.Position(1)=ax1.Position(1);
            txt1.Position(2)=slider1.Position(2)-txt1.Position(4);
            txt1.Position(3)=ax1.Position(3);
            
            slider2.Position(1)=ax2.Position(1);
            slider2.Position(3)=ax2.Position(3);
            txt2.Position(1)=ax2.Position(1);
            txt2.Position(2)=slider1.Position(2)-txt1.Position(4);
            txt2.Position(3)=ax2.Position(3);
          
            ShiftUp1=slider1.Position(2)+slider1.Position(4);
            ax1.OuterPosition(2)=ShiftUp1;
            ax1.OuterPosition(4)=figHeight-ShiftUp1;
            
            ShiftUp2=slider2.Position(2)+slider2.Position(4);
            ax2.OuterPosition(2)=ShiftUp2;
            ax2.OuterPosition(4)=figHeight-ShiftUp2;
            
            edit3.Position(1)=ax2.OuterPosition(1)+ax2.OuterPosition(3);
        end
        
        function updateplot1(~,~)
            sVal = get(slider1,'value');
            set(ploth1,'YData',theta_tSlide(sVal));
            set(txt1,'String',sprintf('tau = %0.1f',...
                tau(round(sVal))));
        end
        function updateplot2(~,~)
            sVal2 = get(slider2,'value');
            set(ploth2,'YData',T_tSlide(sVal2));
            set(txt2,'String',sprintf('t = %0.1f s',...
                t(round(sVal2))));
        end
        function editCallback3(~,~)
            sVal3 = get(edit3,'value');
            T0 = sVal3;
            sVal2 = get(slider2,'value');
            T=theta2*(T1-T0)+T0;
            x=etaR*b;
            t=tau*b^2/alpha;
            
            T_tSlide = @(sVal2)T(round(sVal2),:);
            set(ploth2,'YData',T_tSlide(sVal2));
            disp('YES');
        end
        
        set(fig1,'Visible','on');
    end

end

function [MinVal,iMinVal]=closest(Vec,Val)
% This function simply evaluates the smallest difference between values Val
% and each value within a vector Vec to determine the closest index
% position. This is used to determine the index of Vec where the value is
% closest to Val, i.e. the difference is smallest ~0.
    diff=abs(Vec-Val);
    [MinVal,iMinVal]=min(diff);
end

function [htext] = label(h,textString,varargin)
%LABEL places a label next to your data.  
% 
% This function provides an option between the legend and text or annotation commands
% for labeling data that you plot.  Edward Tufte
% says that data shouldn't stray far from its label, because
% the viewer of a graph should not need to repeatedly move his or her eyes
% back and forth between plotted data and the legend to connect the dots of
% which data are which.  In this spirit, label can be used to place a
% label directly on a plot close to the data it describes.  
%
%% Syntax 
% 
%  label(h,'string')
%  label(...,'location',LocationString)
%  label(...,'TextProperty',PropertyValue)
%  label(...,'slope')
%  h = label(...)
%
%% Description 
% 
% label(h,'string') places 'string' near the leftmost data described by
% handle h. 
%
% label(...,'location',LocationString) specifies location of the string.
% LocationString can be any of the following:
% 
% * 'left' or 'west' (default) 
% * 'right' or 'east' 
% * 'top' or 'north' 
% * 'bottom' or 'south' 
% * 'center' or 'middle' 
% 
% label(...,'TextProperty',PropertyValue) specifies text properties as
% name-value pairs. 
%
% label(...,'slope') attempts to angle text following the local slope of
% the data. 
%
% htext = label(...) returns the handle htext of the newly-created text
% object. 
% 
%% Author Info
% Written by Chad A. Greene of the University of Texas at Austin and its 
% Institute for Geophysics, July 2014. 
% Fixed for R2014b in January 2015. 
% 
% See also annotation, text, and legend. 
%% Initial input error checks: 

assert(ishandle(h)==1,'Unrecognized object handle.')
assert(isempty(get(0,'children'))==0,'No current axes are open.') 
assert(isnumeric(textString)==0,'Label given by textString must be a string.') 
assert(nargin>=2,'Must input an object handle and corresponding label string.') 

%% Set defaults: 

location = 'left'; 
followSlope = false; 

%% Parse inputs

tmp = strncmpi(varargin,'loc',3); 
if any(tmp)
    location = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
end

tmp = strcmpi(varargin,'slope'); 
if any(tmp) 
    followSlope = true; 
    varargin = varargin(~tmp); 
end


%% 

color = get(h,'color'); 
xdata = get(h,'XData'); 
assert(isvector(xdata)==1,'Plotted data must be vector or scalar.') 
ydata = get(h,'YData'); 

gcax = get(gca,'xlim'); 
gcay = get(gca,'ylim'); 

if followSlope
    pbp = kearneyplotboxpos(gca); % A modified version of Kelly Kearney's plotboxpos function is included as a subfunction below.  

    % slope is scaled because of axes and plot box may not be equal and square:
    gcaf = pbp(4)/pbp(3); 
    apparentTheta = atand(gcaf*gradient(ydata,xdata).*(gcax(2)-gcax(1))/(gcay(2)-gcay(1)));

end

% Find indices of data within figure window: 
ind = find(xdata>=gcax(1)&xdata<=gcax(2)&ydata>=gcay(1)&ydata<=gcay(2)); 

switch lower(location)
    case {'left','west','leftmost','westmost'}
        horizontalAlignment = 'left'; 
        verticalAlignment = 'bottom'; 
        x = min(xdata(ind));
        y = ydata(xdata==x);
        textString = [' ',textString]; 
        xi = xdata==x; 
        
    case {'right','east','rightmost','eastmost'}
        horizontalAlignment = 'right'; 
        verticalAlignment = 'bottom'; 
        x = max(xdata(ind)); 
        y = ydata(xdata==x);
        textString = [textString,' ']; 
        xi = xdata==x(1); 
        
    case {'top','north','topmost','northmost'}
        horizontalAlignment = 'left'; 
        verticalAlignment = 'top'; 
        y = max(ydata(ind));
        x = xdata(ydata==y);
        xi = xdata==x(1); 
        
    case {'bottom','south','southmost','bottommost'} 
        horizontalAlignment = 'left'; 
        verticalAlignment = 'bottom'; 
        y = min(ydata(ind));
        x = xdata(ydata==y);
        xi = xdata==x(1); 
        
    case {'center','central','middle','centered','middlemost','centermost'}
        horizontalAlignment = 'center'; 
        verticalAlignment = 'bottom'; 
        xi = round(mean(ind)); 
        if ~ismember(xi,ind)
            xi = find(ind<xi,1,'last'); 
        end
        x = xdata(xi); 
        y = ydata(xi);
        
        
    otherwise
        error('Unrecognized location string.') 
end
 
% Set rotation preferences: 
if followSlope
    theta = apparentTheta(xi); 
else
    theta = 0; 
end


% Create the label: 
htext = text(x(1),y(1),textString,'color',color,'horizontalalignment',horizontalAlignment,...
    'verticalalignment',verticalAlignment,'rotation',theta); 

% Add any user-defined preferences: 
if length(varargin)>1 
    set(htext,varargin{:});
end


% Clean up: 
if nargout == 0
    clear htext
end

end

function pos = kearneyplotboxpos(h)
    %PLOTBOXPOS Returns the position of the plotted axis region. THIS IS A
    %SLIGHTLY MODIFIED VERSION OF KELLY KEARNEY'S PLOTBOXPOS FUNCTION. 
    %
    % pos = plotboxpos(h)
    %
    % This function returns the position of the plotted region of an axis,
    % which may differ from the actual axis position, depending on the axis
    % limits, data aspect ratio, and plot box aspect ratio.  The position is
    % returned in the same units as the those used to define the axis itself.
    % This function can only be used for a 2D plot.  
    %
    % Input variables:
    %
    %   h:      axis handle of a 2D axis (if ommitted, current axis is used).
    %
    % Output variables:
    %
    %   pos:    four-element position vector, in same units as h

    % Copyright 2010 Kelly Kearney

    % Check input

    if nargin < 1
        h = gca;
    end

    if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
        error('Input must be an axis handle');
    end

    % Get position of axis in pixels

    currunit = get(h, 'units');
    set(h, 'units', 'pixels');
    axisPos = get(h, 'Position');
    set(h, 'Units', currunit);

    % Calculate box position based axis limits and aspect ratios

    darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
    pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

    if ~darismanual && ~pbarismanual

        pos = axisPos;

    else

        dx = diff(get(h, 'XLim'));
        dy = diff(get(h, 'YLim'));
        dar = get(h, 'DataAspectRatio');
        pbar = get(h, 'PlotBoxAspectRatio');

        limDarRatio = (dx/dar(1))/(dy/dar(2));
        pbarRatio = pbar(1)/pbar(2);
        axisRatio = axisPos(3)/axisPos(4);

        if darismanual
            if limDarRatio > axisRatio
                pos(1) = axisPos(1);
                pos(3) = axisPos(3);
                pos(4) = axisPos(3)/limDarRatio;
                pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
            else
                pos(2) = axisPos(2);
                pos(4) = axisPos(4);
                pos(3) = axisPos(4) * limDarRatio;
                pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
            end
        elseif pbarismanual
            if pbarRatio > axisRatio
                pos(1) = axisPos(1);
                pos(3) = axisPos(3);
                pos(4) = axisPos(3)/pbarRatio;
                pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
            else
                pos(2) = axisPos(2);
                pos(4) = axisPos(4);
                pos(3) = axisPos(4) * pbarRatio;
                pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
            end
        end
    end

    % Convert plot box position to the units used by the axis

    temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', get(h, 'parent'));
    % set(temp, 'Units', currunit); % <-This line commented-out by Chad Greene, specifically for label function.  
    pos = get(temp, 'position');
    delete(temp);
end
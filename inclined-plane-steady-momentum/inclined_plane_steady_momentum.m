function handles = inclined_plane_steady_momentum
    % parameter values
    delta = 2.5/1000; % m; film thickness
    beta  = 0.5;      % --; angle of inclination
    rho   = 0.8e+3;   % kg/m^3; density
    nu    = 2e-4;     % m2/s; kinematic viscosity
    mu    = nu*rho;   % Pa.s; dynamic viscosity
    g     = 9.81;     % m/s^2; gravitational constant
    x     = linspace(0,delta,60); % vector of thickness of film [0,delta]

    % calcs solution array for velocity profile using plotValFcn
    plotVal = plotValFcn('mu',mu); 
    plotMax = max(plotVal'); % find maximum value to set limit of y-axis

    createGUI; % calls createGUI to built interactive plot and sliders
    %----------------------------------------------------------------------
    function fun = plotValFcn(varStr,value)
        eval(sprintf('%s = %d;',varStr,value));
        fun = (rho*g*delta^2*cos(beta)/(2*mu))*(1-(x/delta).^2);
    end
    %----------------------------------------------------------------------
    function createGUI
        hFig1 = figure('ResizeFcn',@figPosition,'Visible','off');
        sVal1 = mu;
        sVal2 = beta;
        sVal3 = rho;
        sVal4 = delta;
        
        ax1    = axes('Parent',hFig1,'Units','pixel');
        hPlot1 = plot(ax1,x,plotVal);
        title('INCLINED MOMENTUM BALANCE');
        xlabel('Distance, x (m)','FontSize',16);
        ylabel('Velocity, v_z(x) (m/s)','FontSize',16);
        ylim([0,plotMax]);
        xlim([0,delta]);
        plotXRange = hPlot1.XData(end) - hPlot1.XData(1);
        plotYRange = hPlot1.YData(end) - hPlot1.YData(1);
        plotXMid   = hPlot1.XData(1) + plotXRange/2;
        ptxt1 = text(plotXMid,0.9*max(hPlot1.YData),'v_z(x) = \rhog\delta^2 cos(\beta)/(2\mu))[1 - (x/\delta)^2]','FontSize',14);
        ptxt1.HorizontalAlignment = 'center';
        
        stxt1 = 'VISCOSITY = %#.3f kg/m.s : ';
        stxt2 = 'INCLINATION ANGLE = %#.2g : ';
        stxt3 = 'DENSITY = %#5.4g kg/m^3 : ';
        stxt4 = 'FILM THICKNESS = %#.4f m : ';
        
        %------------------------------------------------------------------
        slider1 = uicontrol('Parent',hFig1,'Style','slider',...
                            'value',sVal1,...
                            'min',0.5*sVal1,'max',1.5*sVal1);
        txt1    = uicontrol('Style','text','String',...
                            sprintf(stxt1,sVal1));   
        hLstn1  = addlistener(slider1,'ContinuousValueChange',@updateplot1);
        %------------------------------------------------------------------
        slider2 = uicontrol('Parent',hFig1,'Style','slider',...
                            'value',sVal2,...
                            'min',0.5*sVal2,'max',1.5*sVal2);
        txt2    = uicontrol('Style','text','String',...
                            sprintf(stxt2,sVal2));   
        hLstn2  = addlistener(slider2,'ContinuousValueChange',@updateplot2);
        %------------------------------------------------------------------
        slider3 = uicontrol('Parent',hFig1,'Style','slider',...
                            'value',sVal3,...
                            'min',0.5*sVal3,'max',1.5*sVal3);
        txt3    = uicontrol('Style','text','String',...
                            sprintf(stxt3,sVal3));   
        hLstn3  = addlistener(slider3,'ContinuousValueChange',@updateplot3);
        %------------------------------------------------------------------
        slider4 = uicontrol('Parent',hFig1,'Style','slider',...
                            'value',sVal4,...
                            'min',0.5*sVal4,'max',1.5*sVal4);
        txt4    = uicontrol('Style','text','String',...
                            sprintf(stxt4,sVal4));   
        hLstn4  = addlistener(slider4,'ContinuousValueChange',@updateplot4);
        %------------------------------------------------------------------
        function figPosition(varargin)
            figWidth  = hFig1.Position(3);
            figHeight = hFig1.Position(4);
            ax1pos    = getpixelposition(ax1);
            ax1left   = ax1pos(1);
            ax1btm    = ax1pos(2);
            ax1width  = ax1pos(3); 
            ax1height = ax1pos(4);
            
            border    = 20;
            
            slider1.Position(1) = ax1left + ax1width/2;
            slider1.Position(2) = border;
            slider1.Position(3) = ax1width/2;

            txt1.Position(4)    = slider1.Position(4);
            txt1.Position(3)    = ax1width/2;
            txt1.Position(1)    = slider1.Position(1) - txt1.Position(3);
            txt1.Position(2)    = slider1.Position(2);
            txt1.HorizontalAlignment = 'right';
            
            slider2.Position    = slider1.Position;
            slider2.Position(2) = slider1.Position(2) + slider1.Position(4);
            txt2.Position       = txt1.Position;
            txt2.Position(2)    = slider2.Position(2);
            txt2.HorizontalAlignment = 'right';
            
            slider3.Position    = slider2.Position;
            slider3.Position(2) = slider2.Position(2) + slider2.Position(4);
            txt3.Position       = txt2.Position;
            txt3.Position(2)    = slider3.Position(2);
            txt3.HorizontalAlignment = 'right';
            
            slider4.Position    = slider3.Position;
            slider4.Position(2) = slider3.Position(2) + slider3.Position(4);
            txt4.Position       = txt3.Position;
            txt4.Position(2)    = slider4.Position(2);
            txt4.HorizontalAlignment = 'right';
            
            ShiftUp   = slider4.Position(2) + slider4.Position(4);

            ax1pos(2) = ShiftUp;
            ax1pos(4) = figHeight - ShiftUp;
            ax1pos(3) = figWidth;
            ax1pos(1) = 0;
  
            handles   = struct('fig',hFig1,'slider1',slider1,'txt1',txt1,'plot',hPlot1);
            set(ax1,'OuterPosition',ax1pos);
            
        end
        %------------------------------------------------------------------
        function updateplot1(~,~)
            sVal1 = get(slider1,'value');
            set(hPlot1,'YData',plotValFcn('mu',sVal1));
            set(txt1,'String',sprintf(stxt1,sVal1));
        end
        %------------------------------------------------------------------
        function updateplot2(~,~)
            sVal2 = get(slider2,'value');
            set(hPlot1,'YData',plotValFcn('beta',sVal2));
            set(txt2,'String',sprintf(stxt2,sVal2));
        end
        %------------------------------------------------------------------
        function updateplot3(~,~)
            sVal3 = get(slider3,'value');
            set(hPlot1,'YData',plotValFcn('rho',sVal3));
            set(txt3,'String',sprintf(stxt3,sVal3));
        end
        %------------------------------------------------------------------
        function updateplot4(~,~)
            sVal4 = get(slider4,'value');
            set(hPlot1,'YData',plotValFcn('delta',sVal4));
            set(txt4,'String',sprintf(stxt4,sVal4));
        end
        %------------------------------------------------------------------
        set(hFig1,'Visible','on');
    end
    %----------------------------------------------------------------------
end

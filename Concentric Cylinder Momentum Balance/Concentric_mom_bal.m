function Concentric_mom_bal       
    % parameter values
    omega_o = 2;
    omega_i = 4;
    Ro      = 0.8;
    Ri      = 0.5; 
    r       = linspace(0.1,5*Ro,100);
    
    % calcs solution array for temp profile using plotValFcn
    plotVal = plotValFcn('Ro',Ro); 
    YMax    = max(plotVal.Y); % find maximum value to set limit of y-axis
    XMax    = max(plotVal.X); % find maximum value to set limit of y-axis

    createGUI; % calls createGUI to built interactive plot and sliders
    %----------------------------------------------------------------------
    function fun = plotValFcn(varStr,value)
        eval(sprintf('%s = %d;',varStr,value));
        fun = Velocity_profile(omega_o,omega_i,Ro,Ri,r);
    end
    %----------------------------------------------------------------------
    function createGUI
        hFig1 = figure('ResizeFcn',@figPosition,'Visible','off');
        
        sVal1 = omega_o;
        sVal2 = omega_i;
        sVal3 = Ro;
        sVal4 = Ri;
        
        sVal1Name = 'omega_o';
        sVal2Name = 'omega_i';
        sVal3Name = 'Ro';
        sVal4Name = 'Ri';
        
        stxt1 = 'OUTTER ANGULAR VELOCITY = %#.3g rad/s : ';
        stxt2 = 'INNER ANGULAR VELOCITY = %#.3g rad/s : ';
        stxt3 = 'OUTTER CYLINDER RADIUS = %#5.3g m : ';
        stxt4 = 'INNER CYLINDER RADIUS = %#.3g m : ';
        
        sRangeMin = 0.1;
        sRangeMax = 2;
        
        ax1    = axes('Parent',hFig1,'Units','pixel');
        hPlot1 = plot(ax1,plotVal.X,plotVal.Y);
        title('CONCENTRIC CYLINDER ANGULAR VELOCITY PROFILE');
        xlabel('r (m) ');
        ylabel('v_{\theta}(r) (m/s) ');
        ylim([0,YMax]);
        xlim([0,XMax]);
        plotXRange = hPlot1.XData(end) - hPlot1.XData(1);
        plotYRange = hPlot1.YData(end) - hPlot1.YData(1);
        plotXMid   = hPlot1.XData(1) + plotXRange/2;
        ptxt1 = text(plotXMid,0.9*max(hPlot1.YData),'$v_{\theta}(r) = \frac{kR}{1-k^2}[(\Omega_o-\Omega_ik^2)(\frac{r}{kR})+(\Omega_i-\Omega_o)(\frac{kR}{r})$','interpreter','latex');
        ptxt1.FontSize = 15;
        ptxt1.HorizontalAlignment = 'center';
        
        %------------------------------------------------------------------
        slider1 = uicontrol('Parent',hFig1,'Style','slider',...
                            'value',sVal1,...
                            'min',sRangeMin*sVal1,'max',sRangeMax*sVal1);
        txt1    = uicontrol('Style','text','String',...
                            sprintf(stxt1,sVal1));   
        hLstn1  = addlistener(slider1,'ContinuousValueChange',@updateplot1);
        %------------------------------------------------------------------
        slider2 = uicontrol('Parent',hFig1,'Style','slider',...
                            'value',sVal2,...
                            'min',sRangeMin*sVal2,'max',sRangeMax*sVal2);
        txt2    = uicontrol('Style','text','String',...
                            sprintf(stxt2,sVal2));   
        hLstn2  = addlistener(slider2,'ContinuousValueChange',@updateplot2);
        %------------------------------------------------------------------
        slider3 = uicontrol('Parent',hFig1,'Style','slider',...
                            'value',sVal3,...
                            'min',sRangeMin*sVal3,'max',sRangeMax*sVal3);
        txt3    = uicontrol('Style','text','String',...
                            sprintf(stxt3,sVal3));   
        hLstn3  = addlistener(slider3,'ContinuousValueChange',@updateplot3);
        %------------------------------------------------------------------
        slider4 = uicontrol('Parent',hFig1,'Style','slider',...
                            'value',sVal4,...
                            'min',sRangeMin*sVal4,'max',sRangeMax*sVal4);
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
            
            shiftUp   = slider4.Position(2) + slider4.Position(4);

            ax1pos(2) = shiftUp;
            ax1pos(4) = figHeight - shiftUp;
            ax1pos(3) = figWidth;
            ax1pos(1) = 0;

            set(ax1,'OuterPosition',ax1pos);
            
        end
        %------------------------------------------------------------------
        function updateplot1(~,~)
            sVal1 = get(slider1,'value');
            plotVal = plotValFcn(sVal1Name,sVal1);
            set(hPlot1,'YData',plotVal.Y);
            set(txt1,'String',sprintf(stxt1,sVal1));
            ylim([0,YMax]);
            xlim([0,XMax]);
        end
        %------------------------------------------------------------------
        function updateplot2(~,~)
            sVal2 = get(slider2,'value');
            plotVal = plotValFcn(sVal2Name,sVal2);
            set(hPlot1,'YData',plotVal.Y);
            set(txt2,'String',sprintf(stxt2,sVal2));
            ylim([0,YMax]);
            xlim([0,XMax]);
        end
        %------------------------------------------------------------------
        function updateplot3(~,~)
            sVal3 = get(slider3,'value');
            plotVal = plotValFcn(sVal3Name,sVal3);
            set(hPlot1,'YData',plotVal.Y);
            set(txt3,'String',sprintf(stxt3,sVal3));
            ylim([0,YMax]);
            xlim([0,XMax]);
        end
        %------------------------------------------------------------------
        function updateplot4(~,~)
            sVal4 = get(slider4,'value');
            plotVal = plotValFcn(sVal4Name,sVal4);
            set(hPlot1,'YData',plotVal.Y);
            set(txt4,'String',sprintf(stxt4,sVal4));
            ylim([0,YMax]);
            xlim([0,XMax]);
        end
        %------------------------------------------------------------------
        set(hFig1,'Visible','on');
    end
    %----------------------------------------------------------------------
end

function out = Velocity_profile(omega_o,omega_i,Ro,Ri,r)
    wi      = omega_i;
    wo      = omega_o;
    k       = Ri/Ro;
    v_theta = zeros(length(r),1);
    
    for i = 1:length(r)
        v_theta(i) = k*Ro/(1-k^2)*((wo-wi*k^2)*r(i)/(k*Ro)+(wi-wo)*(k*Ro/r(i)));
    end
    
    out.Y = v_theta;
    out.X = r;
end
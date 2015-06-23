function Inclined_mom_bal
    % parameter values
    delta = 1;
    beta = 0.5;
    mu = 1.4e-5;
    rho = 1000;
    g = 9.81;

    % discrete vector along thickness of film [0,delta]
    x = linspace(0,delta,60);

    % calcs solution array for velocity profile using v_zSol
    v_z = v_zSol('mu',mu); 
    v_z0 = max(v_z'); % find maximum value to set limit of y-axis

    createGUI; % calls createGUI to built interactive plot and slider

    function fun = v_zSol(varStr,value)
        eval(sprintf('%s = %d;',varStr,value));
        fun = (rho*g*delta^2*cos(beta)/(2*mu))*(1-(x/delta).^2);
    end

    function createGUI
       
        fig1 = figure('ResizeFcn',@figPosition,'Visible','off');
        sVal1 = mu;
        sVal2 = beta;
        sVal3 = rho;
        sVal4 = delta;
        
        ax1 = axes('Parent',fig1,'Units','pixel');
        ploth1 = plot(ax1,x,v_z);
        title('INCLINED MOMENTUM BALANCE');
        xlabel('Distance x (m)');
        ylabel('Velocity v_z(x) (m/s)');
        ylim([0,v_z0]);
        xlim([0,delta]);
        
        slider1 = uicontrol('Parent',fig1,'Style','slider',...
            'value',sVal1,...
            'min',0.5*sVal1,'max',1.5*sVal1);
        txt1 = uicontrol('Style','text','String',...
            sprintf('viscosity = %0.1d kg/m.s',sVal1));   
        hLstn1 = addlistener(slider1,'ContinuousValueChange',@updateplot1);
        
        slider2 = uicontrol('Parent',fig1,'Style','slider',...
            'value',sVal2,...
            'min',0.5*sVal2,'max',1.5*sVal2);
        txt2 = uicontrol('Style','text','String',...
            sprintf('beta = %0.1d',sVal2));   
        hLstn2 = addlistener(slider2,'ContinuousValueChange',@updateplot2);

        slider3 = uicontrol('Parent',fig1,'Style','slider',...
            'value',sVal3,...
            'min',0.5*sVal3,'max',1.5*sVal3);
        txt3 = uicontrol('Style','text','String',...
            sprintf('density = %0.1d kg/m^3',sVal3));   
        hLstn3 = addlistener(slider3,'ContinuousValueChange',@updateplot3);
        
        slider4 = uicontrol('Parent',fig1,'Style','slider',...
            'value',sVal4,...
            'min',0.5*sVal4,'max',1.5*sVal4);
        txt4 = uicontrol('Style','text','String',...
            sprintf('film thickness = %0.1d m',sVal4));   
        hLstn4 = addlistener(slider4,'ContinuousValueChange',@updateplot4);
        
      
        function figPosition(varargin)
            ax1pos=getpixelposition(ax1);
             
            slider1.Position(1)=ax1pos(1);
            slider1.Position(3)=ax1pos(3);

            txt1.Position(1)=ax1pos(1);
            txt1.Position(2)=slider1.Position(2)-txt1.Position(4);
            txt1.Position(3)=ax1pos(3);
            
            txt2.Position=txt1.Position;
            txt2.Position(2)=slider1.Position(2)+slider1.Position(4);
            slider2.Position=slider1.Position;
            slider2.Position(2)=txt2.Position(2)+txt2.Position(4);
            
            txt3.Position=txt2.Position;
            txt3.Position(2)=slider2.Position(2)+slider2.Position(4);
            slider3.Position=slider2.Position;
            slider3.Position(2)=txt3.Position(2)+txt3.Position(4);
            
            txt4.Position=txt3.Position;
            txt4.Position(2)=slider3.Position(2)+slider3.Position(4);
            slider4.Position=slider3.Position;
            slider4.Position(2)=txt4.Position(2)+txt4.Position(4);
            
            figWidth=fig1.Position(3);
            figHeight=fig1.Position(4);
            
            ShiftUp=slider4.Position(2)+slider4.Position(4);

            ax1pos(2)=ShiftUp;
            ax1pos(4)=figHeight-ShiftUp;
            ax1pos(3)=figWidth;
            ax1pos(1)=0;
  
            set(ax1,'OuterPosition',ax1pos);
        end
        
        function updateplot1(~,~)
            sVal1 = get(slider1,'value');
            set(ploth1,'YData',v_zSol('mu',sVal1));
            set(txt1,'String',sprintf('viscosity = %0.1d kg/m.s',sVal1));
        end
        
        function updateplot2(~,~)
            sVal2 = get(slider2,'value');
            set(ploth1,'YData',v_zSol('beta',sVal2));
            set(txt2,'String',sprintf('beta = %0.1d',sVal2));
        end
        
        function updateplot3(~,~)
            sVal3 = get(slider3,'value');
            set(ploth1,'YData',v_zSol('rho',sVal3));
            set(txt3,'String',sprintf('density = %0.1d kg/m^3',sVal3));
        end
        
        function updateplot4(~,~)
            sVal4 = get(slider4,'value');
            set(ploth1,'YData',v_zSol('delta',sVal4));
            set(txt4,'String',sprintf('film thickness = %0.1d m',sVal4));
        end
        
        set(fig1,'Visible','on');
    end
end

# Inclined Momentum Balance

## Inclined_mom_bal.m
This code plots the solution to the inclined plane momentum balance of a falling film with an interactive slider to allow the user to change some of the key parameters.

### The assumptions:
  1. constant density and viscosity
  2. steady-state
  3. laminar flow (simple shear flow)
  4. fully developed flow
  5. newton's law of viscosity is applicable

### Derivation:

    d(tau)/dx=rho*g*cos(beta) reduces to: d^2(v_z)/dx^2=-(rho*g*cos(beta)/mu)
  
  where: 
  
    rho=density of fluid
    g=gravity
    beta=angle of inclination w.r.t the vertical axis
    mu=dynamic viscosity of fluid
    v_z=viscosity of fluid in the z-direction
    z=direction of flow (parallel to plane)
    x=direction orthogonal to plane
  
  After integrating twice:
  
    v_z=-(rho*g*cos(beta)/mu)*(x^2/2)+c1*x+c2
  
  Boundary conditions:
  
    d(v_z(x=0))/dx=0 
    v_z(x=delta)=0
  
  Final solution:
  
    v_z=(rho*g*delta^2*cos(beta)/(2*mu))*(1-(x/delta)^2)
  
### The MATLAB code:
```
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
        
        ax1 = subplot(2,1,1,'Parent',fig1,...
            'Units','normalized');
        ploth1 = plot(ax1,x,v_z);
        xlabel('Distance x (m)');
        ylabel('Velocity v_z(x) (m/s)');
        ylim([0,v_z0]);
        xlim([0,delta]);

        subplot(2,1,2,'visible','off');
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
            sprintf('viscosity = %0.1d kg/m.s',sVal2));   
        hLstn2 = addlistener(slider2,'ContinuousValueChange',@updateplot2);

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
        
        set(fig1,'Visible','on');
    end
end
```
  

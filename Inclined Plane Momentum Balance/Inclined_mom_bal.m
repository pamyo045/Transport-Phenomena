function Inclined_mom_bal
%   Solution to the inclined plane momentum balance
%   of a falling film of thickness delta
%   with the following assumptions:
%   1) constant density, viscosity
%   2) steady-state
%   3) laminar flow (simple shear flow)
%   4) fully developed flow
%   5) newton's law of viscosity is applicable
%   d(tau)/dx=rho*g*cos(beta) reduces to:
%   d^2(v_z)/dx^2=-(rho*g*cos(beta)/mu)
%   where: rho=density of fluid, g=gravity,
%   beta=angle of inclination w.r.t the vertical axis
%   mu=dynamic viscosity of fluid
%   v_z=viscosity of fluid in the z-direction
%   z=direction of flow (parallel to plane)
%   x=direction orthogonal to plane
%   after integrating twice:
%   v_z=-(rho*g*cos(beta)/mu)*(x^2/2)+c1*x+c2
%   boundary conditions of d(v_z(x=0))/dx=0 and v_z(x=delta)=0:
%   v_z=(rho*g*delta^2*cos(beta)/(2*mu))*(1-(x/delta)^2)

%parameter values
delta = 1;
beta = 0.5;
mu = 1.4e-5;
rho = 1000;
g = 9.81;

%discrete vector along thickness of film [0,delta]
x = linspace(0,delta,60);

v_z = v_zSol(mu); %calcs solution array for velocity profile using v_zSol

v_z0 = max(v_z'); %find maximum value to set limit of y-axis

createGUI; %calls createGUI to built interactive plot and slider

    function createGUI
        fig1 = figure;
        sVal = mu;
        ax1 = axes('Parent',fig1,'position',[0.13 0.39  0.77 0.54]);
        ploth1 = plot(ax1,x,v_z);
        title('v_z vs t');
        xlabel('Distance x (m)');
        ylabel('Velocity v_z(x) (m/s)');
        ylim([0,v_z0]);
        xlim([0,delta]);

        slider1 = uicontrol('Parent',fig1,'Style','slider',...
            'Position',[81 54 419 23],...
            'value',sVal,...
            'min',0.5*mu,'max',1.5*mu);
        
        txt1 = uicontrol('Style','text',...
            'Position',[81 (54-25) 419 23],...
            'String',sprintf('Viscosity = %0.1d kg/m.s',mu));

        hLstn1 = addlistener(slider1,'ContinuousValueChange',@updateplot1);

        function updateplot1(~,~)
            sVal = get(slider1,'value');
            set(ploth1,'YData',v_zSol(sVal));
            set(txt1,'String',sprintf('Viscosity = %0.1d kg/m.s',sVal));
        end
    end

    function fun = v_zSol(mu)
        fun = (rho*g*delta^2*cos(beta)/(2*mu))*(1-(x/delta).^2);
    end

end

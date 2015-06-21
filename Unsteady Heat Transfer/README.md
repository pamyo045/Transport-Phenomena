# Unsteady Heat Transfer

## Dimless_unsteadyT.m
Solves the 1D unsteady-state temperature conduction heat equation in dimentionless form.

### Derivation:
  1-D Unsteady heat conduction equation for a uniform rectangular plane:
  
    d(T)/d(t)=alpha*d^2(T)/(d(x))^2 

  Where:

    T     =temperature as a function of x and t (C)
    t     =time (s)
    x     =position (m)
    alpha =thermal diffusivity (m^2/s)
    T0    =initial temperature of the slab
    T1    =temperature imposed at the slab surfaces for time > 0
    b     =plane slab thickness divided by 2 (i.e. slab thickness = 2b, x=0 is the center)

  Dimensionless variable transformation:
  
    theta=(T1-T)/(T1-T0)  *dimensionless temperature*
    eta=x/b               *dimensionless position*
    tau=alpha*t/b^2       *dimensionless time* 

  Thus:
  
    d(theta)/d(tau)=d^2(theta)/d((eta))^2 

  Initial & boundary conditions:
  
    IC:   theta=1 when tau=0
    BC1:  theta=0 when eta=-1 for tau>0
    BC2:  theta=0 when eta=+1 for tau>0

### The MATLAB code:
```
function Dimless_unsteadyT  
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
    
    figure;
    hold on
    
    iterationCount=0;
    tauPlot=[0.01 0.04 0.1 0.2 0.4 0.6 1.0];
    
    % the array ind stores the index locations of tauPlot values
    ind=zeros(1,length(tauPlot));
    for i=1:length(tauPlot)
        [~,ind(i)]=closest(tau,tauPlot(i)); % function that finds the index 
                                            % whose value closest matching 
                                            % in tauPlot
    end
    
    % loop that plots the solution for various tauPlot values, serves as a
    % comparison to the solutions of Bird et al. to show the accuracy of the
    % solution found using pdepe and the analytical one from the authors.
    for i=1:length(ind)
        k=ind(i);
        iterationCount=iterationCount+1;
        ph=plot(etaR,theta2(k,:));
        label(ph,sprintf('\\tau =%0.2f',tau(k)),'location','center','slope');
    end
    hold off
    
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

    %defines the ICs for the model defined in pdefun
    function theta0 = pdeic(eta)
        theta0 = 1;      
    end

    %defines the BCs for the model defined in pdefun
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

        fig = figure; % creates new (empty figure) window 
        sVal = 1;   % slider index, used to initialize value of indices 
                    % which represent the index value of the dependant 
                    % variable array, eg. c_g(i,j)=f(t,x), (i or j)=sindex 
                    % depending on the slider wishing to display.

        theta_tSlide = @(sVal)theta2(round(sVal),:);
        
        subplot(2,1,1); 

        hold on
        ploth1 = plot(etaR,theta_tSlide(sVal));
        hold off
        xlabel('\eta =x/b');
        ylabel('\theta =(T-T_0)/(T_1-T_0)');
        ylim([0,1]);
        xlim([0,1]);
        ax = gca; % setting gca (MATLAB keyword for current active
                  % axis) to ax for changing layout
        ax.Position(2) = ax.Position(2)-0.10; %bottom
        ax.Position(4) = ax.Position(4)+0.10; %height
        box on
        panel1 = uipanel('Parent',fig,'Units','normalized',...
            'Position',[0.5-0.2/2 0.0 0.2 0.1]);
        panel2 = uipanel('Parent',fig,'Units','normalized',...
            'Position',[0.0 0.1 1.0 0.1]);
        
        pos_f=ax.Position; %figure position
        pos_p2=panel2.Position; %panel position
        
        slider1 = uicontrol(panel2,'Style','slider',...
            'Units','Normalized',...
            'value',sVal,...
            'Position',[pos_f(1) 0.2 pos_f(3) 0.6],...
            'min',1,'max',length(tau));

        txt1 = uicontrol(panel1,'Style','text',...
            'Units','Normalized',...
            'FontName','symbol',...
            'FontSize',12,...
            'Position',[0, 0, 0.5, 1],...
            'String','t');
       
        txt2 = uicontrol(panel1,'Style','text',...
            'Units','Normalized',...
            'FontSize',12,...
            'Position',[0.5, 0, 0.5, 1],...
            'String',sprintf('= %0.1f',tau(round(sVal))));

        hLstn1 = addlistener(slider1,'ContinuousValueChange',@updateplot1);

        function updateplot1(~,~)
            sVal = get(slider1,'value');
            set(ploth1,'YData',theta_tSlide(sVal));
            set(txt2,'String',sprintf('= %0.1f',...
                tau(round(sVal))));
        end
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
```

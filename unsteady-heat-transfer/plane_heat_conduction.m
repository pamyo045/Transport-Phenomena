function plane_heat_conduction       
    % parameter values
    tend  = 1000;   % s; time
    alfa  = 3e-3;   % Pa.s; thermal diffusivity
    width = 6;      % m; slab total width
    tau   = alfa*tend/width^2;
    % calcs solution array for temp profile using plotValFcn
    plotVal = plotValFcn('tend',tend); 
    YMax = max(plotVal.Y); % find maximum value to set limit of y-axis
    XMax = max(plotVal.X);
    
    createGUI; % calls createGUI to built interactive plot and sliders
    %----------------------------------------------------------------------
    function fun = plotValFcn(varStr,value)
        eval(sprintf('%s = %d;',varStr,value));
        fun = transient_HT(tend,alfa,width);
    end
    %----------------------------------------------------------------------
    function createGUI
        hFig1 = figure('ResizeFcn',@figPosition,'Visible','off');
        sVal1 = tend;
        sVal2 = alfa;
        sVal3 = width;
        sVal4 = width;
        
        sVal1Name = 'tend';
        sVal2Name = 'alfa';
        sVal3Name = 'width';
        sVal4Name = 'width';
        
        stxt1 = 'TIME = %.1f s : ';
        stxt2 = 'THERMAL DIFFUSIVITY = %.1e : ';
        stxt3 = 'WIDTH = %.2f m : ';
        
        sRangeMin = 0.01;
        sRangeMax = 5;
        
        ax1    = axes('Parent',hFig1,'Units','pixel');
        hPlot1 = plot(ax1,plotVal.X,plotVal.Y);
        title('PLANE SLAB TRANSIENT HEAT CONDUCTION');
        xlabel('\eta');
        ylabel('\theta(\tau,\eta)');
        ylim([0,YMax]);
        %xlim([0,delta]);
        plotXRange = hPlot1.XData(end) - hPlot1.XData(1);
        plotYRange = hPlot1.YData(end) - hPlot1.YData(1);
        XMid   = hPlot1.XData(1) + plotXRange/2;
        ptxt1 = text(XMid,0.9*YMax,'$\frac{d\theta}{d\tau} = \frac{d^2\theta}{d\eta^2}$','interpreter','latex');
        ptxt1.FontSize = 20;
        ptxt1.HorizontalAlignment = 'center';
        ptxt2 = text(XMid,0.8*YMax,sprintf('\\tau = %.3g',tau));
        ptxt2.HorizontalAlignment = 'center';
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
%         slider4 = uicontrol('Parent',hFig1,'Style','slider',...
%                             'value',sVal4,...
%                             'min',sRangeMin*sVal4,'max',sRangeMax*sVal4);
%         txt4    = uicontrol('Style','text','String',...
%                             sprintf(stxt4,sVal4));   
%         hLstn4  = addlistener(slider4,'ContinuousValueChange',@updateplot4);
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
            
%             slider4.Position    = slider3.Position;
%             slider4.Position(2) = slider3.Position(2) + slider3.Position(4);
%             txt4.Position       = txt3.Position;
%             txt4.Position(2)    = slider4.Position(2);
%             txt4.HorizontalAlignment = 'right';
            
            shiftUp   = slider3.Position(2) + slider3.Position(4);

            ax1pos(2) = shiftUp;
            ax1pos(4) = figHeight - shiftUp;
            ax1pos(3) = figWidth;
            ax1pos(1) = 0;

            set(ax1,'OuterPosition',ax1pos);
            
        end
        %------------------------------------------------------------------
        function update_disp
            tau = alfa*tend/width^2;
            ptxt2.String = sprintf('\\tau = %.3g',tau);
        end
        %------------------------------------------------------------------
        function updateplot1(~,~)
            sVal1 = get(slider1,'value');
            plotVal = plotValFcn(sVal1Name,sVal1);
            set(hPlot1,'YData',plotVal.Y);
            set(txt1,'String',sprintf(stxt1,sVal1));
            update_disp;
        end
        %------------------------------------------------------------------
        function updateplot2(~,~)
            sVal2 = get(slider2,'value');
            plotVal = plotValFcn(sVal2Name,sVal2);
            set(hPlot1,'YData',plotVal.Y);
            set(txt2,'String',sprintf(stxt2,sVal2));
            update_disp;
        end
        %------------------------------------------------------------------
        function updateplot3(~,~)
            sVal3 = get(slider3,'value');
            plotVal = plotValFcn(sVal3Name,sVal3);
            set(hPlot1,'YData',plotVal.Y);
            set(txt3,'String',sprintf(stxt3,sVal3));
            update_disp;
        end
        %------------------------------------------------------------------
        function updateplot4(~,~)
            sVal4 = get(slider4,'value');
            plotVal = plotValFcn(sVal4Name,sVal4);
            set(hPlot1,'YData',plotVal.Y);
            set(txt4,'String',sprintf(stxt4,sVal4));
            update_disp;
        end
        %------------------------------------------------------------------
        set(hFig1,'Visible','on');
    end
    %----------------------------------------------------------------------
end

function out = transient_HT(tend,alfa,width)
    % temporal domain
    nt   = 2;
    tau  = [0 tend*alfa/width^2];
    dtau = abs(tau(2) - tau(1));

    % spatial domain
    nz   = 200;
    eta  = linspace(-1,1,nz);
    deta = abs(eta(2) - eta(1));

    r    = dtau/deta^2;

    % initial conditions
    theta0 = 1; % @tau = 0

    % boundary conditions
    g1  = 0;    % @eta = -1;
    g2  = 0;    % @eta = +1;

    % initialize
    theta  = zeros(nt,nz);
    A = zeros(nz,1);
    B = zeros(nz,1);
    C = zeros(nz,1);
    D = zeros(nz,1);

    for j = 1:nt
        if(j == 1)
            theta(j,:)  = theta0;
            thetaOld    = theta0*ones(nz,1);
        else
            thetaOld   = theta(j-1,:)';
        end

        for i = 1:nz
            if(i == 1)      % lower boundary            
                A(i) =  0;
                B(i) =  1;
                C(i) =  0;
                D(i) =  g1;         
            elseif(i == nz) % upper boundary            
                A(i) =  0;
                B(i) =  1;
                C(i) =  0;
                D(i) =  g2;            
            else            % interior nodes  
                A(i) = -r;
                B(i) =  1 + 2*r;
                C(i) = -r;
                D(i) =  thetaOld(i);
            end  
        end 

        thetaNew   = TDMA(nz,A,B,C,D);
        theta(j,:) = (g1 - thetaNew)./(g1 - theta0);
    end  
    theta = 1 - theta(end,:);

    out.Y = theta;
    out.X = eta;
end

function E = TDMA(L,A,B,C,D)
    E     = zeros(L,1);
    BETA  = zeros(L,1);
    GAMMA = zeros(L,1);
    
    BETA(1)  = B(1);
    GAMMA(1) = D(1)/BETA(1);
    
 
   for I=2:L
        BETA(I)  = B(I) - A(I)*C(I - 1)/BETA(I - 1);
        GAMMA(I) = (D(I) - A(I)*GAMMA(I - 1))/BETA(I);
    end

    E(L) = GAMMA(L);
    for K = 1:L-1
        I    = L - K;
        E(I) = GAMMA(I) - C(I)*E(I + 1)/BETA(I);
    end
end

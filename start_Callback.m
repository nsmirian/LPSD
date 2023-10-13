properties (Access = public)
        arduinoObj % Description
        message
    end
    
    properties (Access = public)
        T % Description
        
    end
    
    properties (Access = public)
        
        region1
%         
        v_thick1 % Description
        v_thick2 % Description
        h_thick1 % Description
        h_thick2 % Description
        
        v_thickness_1
        v_thickness_2
        h_thickness_1
        h_thickness_2
        v_or_h_array
        v_or_h
        amplitude_array
        exp_counter
        
    end
    
  
    
    properties (Access = public)
        amplitude % Description
    end
    
    methods (Access = private)
        
      
   function mycallback(app,src,event)
   display(event.SelectedOption);
   end
   
    end
    function startupFcn(app)
    
     clc
            % Read experiment data from a CSV file
            
            
            % Plot patch on uiaxes
            hold(app.UIAxes, 'on')
            
           % region1 = patch(app.UIAxes,[-10 10 10 -10],[-5 -5 -4.4 -4.4],'r','FaceAlpha',1,...
              %'LineWidth',0.01,'LineStyle','-','tag','region1');
            
            load_folder = "C:\Users\student\Desktop\FRANCK\thesis\excel_data\";
            load_name   = "excel_data.xlsx";
            load_addr   = load_folder + load_name;
            
            app.T = readtable(load_addr,'NumHeaderLines',1);
            
            app.exp_counter = 1;
            app.v_thickness_1   = app.T.Var1;
            app.v_thickness_2   = app.T.Var2;
            app.h_thickness_1   = app.T.Var3;
            app.h_thickness_2   = app.T.Var4;
            app.amplitude_array = app.T.Var5;
            app.v_or_h_array    = app.T.Var6;
            
           app.v_thick1 = app.v_thickness_1(app.exp_counter);
           app.v_thick2 = app.v_thickness_2(app.exp_counter);
            
           app.h_thick1 = app.h_thickness_1(app.exp_counter);
           app.h_thick2 = app.h_thickness_2(app.exp_counter);
            
           app.v_or_h = app.v_or_h_array(app.exp_counter);
            
                      
            %Vertical line
            if app.v_or_h == 0
              app.region1 = patch(app.UIAxes, ...
                                [app.v_thick1 app.v_thick2 app.v_thick2 app.v_thick1], ...
                                [-10 -10 10 10],'r', ...
                                'FaceAlpha',1,...
                                'LineWidth',0.01, ...
                                'LineStyle','-','tag','region1'); 
            
            %Horizontal line
            elseif app.v_or_h == 1                
                app.region1 = patch(app.UIAxes,[-10 10 10 -10], ...
                                 [app.h_thick1 app.h_thick1 app.h_thick2 app.h_thick2], ...
                                 'r','FaceAlpha',1,...
                                 'LineWidth',0.01, ...
                                 'LineStyle','-','tag','region1');
            end
           
            % Define pointer behavior for patch
            pm.enterFcn = @(~,~) cursorPositionFeedback(app, app.region1, 'in');
            pm.exitFcn  = @(~,~) cursorPositionFeedback(app, app.region1, 'out');
            pm.traverseFcn = [];
            iptSetPointerBehavior(app.region1, pm)
            % Enable pointer manager for app
            iptPointerManager(app.UIFigure,'enable');
            
            % Create the Arduino serial object
            app.arduinoObj = serialport("COM3", 9600);
            configureTerminator(app.arduinoObj,"CR/LF");
            %flush(app.arduinoObj);
            % 
            for i=1:8 
                app.message = readline(app.arduinoObj);
                disp(app.message)
            end
function cursorPositionFeedback(app, hobj, inout)
% When inout is 'in', change hobj facecolor to green and update textbox.
% When inout is 'out' change hobj facecolor to red, and clear textbox.
% Check tag property of hobj to identify the object.
switch lower(inout)
    case 'in'
        facecolor = 'g';
        txt = 'Motor(s) vibrating';
        pointer = 'fleur';
        writeline(app.arduinoObj, "4&MOTOR_1_2&0!");
%         message = readline(app.arduinoObj);
%         disp(message)
    case 'out'
        facecolor = 'r';
        txt = 'No';
        pointer = 'crosshair';
        writeline(app.arduinoObj, "0&NO_MOTOR&0!");          %% THIS IS THE LINE THAT İS SUPPOSED TO STOP İT
%         message = readline(app.arduinoObj);
%         disp(message)
        
end
hobj.FaceColor = facecolor;
app.TextAreaEditField.Value = txt;
set(app.UIFigure, 'Pointer', pointer)
end  
   % Determine if mouse is within uiaxes
            cp = app.UIFigure.CurrentPoint;
            isInAxes = cp(1) >= app.UIAxes.Position(1) && ...
                cp(1) <= sum(app.UIAxes.Position([1,3])) && ...
                cp(2) >= app.UIAxes.Position(2) && ...
                cp(2) <= sum(app.UIAxes.Position([2,4]));
            if isInAxes
                set(app.CurrentPositionEditField, 'Value',...
                    sprintf('%.2f,  %.2f', app.UIAxes.CurrentPoint(1,1:2)))
            else
                set(app.CurrentPositionEditField, 'Value', '')
            end
            
    function NEXTButton_2Pushed(app, event)
            
            
            uiconfirm(app.UIFigure,'Are You sure?','Confirm Close',...
            'CloseFcn',@(src,event)mycallback(app,src,event));
            
            app.exp_counter = app.exp_counter + 1;
            app.v_or_h = app.v_or_h_array(app.exp_counter);
            
    
            if ishandle(app.region1)
               delete(app.region1);
           end
            
            %Vertical line
            if app.v_or_h == 0
                app.region1 = patch(app.UIAxes,...
                     [app.v_thick1 app.v_thick2 app.v_thick2 app.v_thick1],...
                     [-10 -10 10 10],'r',...
                     'FaceAlpha',1,...
                     'LineWidth',0.01,...
                     'LineStyle','-','tag','region1'); 
           
            %Horizontal line
            elseif app.v_or_h == 1                
                app.region1 = patch(app.UIAxes,[-10 10 10 -10],...
                     [app.h_thick1 app.h_thick1 app.h_thick2 app.h_thick2],...
                     'r','FaceAlpha',1,...
                     'LineWidth',0.01,...
                     'LineStyle','-','tag','region1');
            end
            
              
            % Define pointer behavior for patch
            pm.enterFcn = @(~,~) cursorPositionFeedback(app, app.region1, 'in');
            pm.exitFcn  = @(~,~) cursorPositionFeedback(app, app.region1, 'out');
            pm.traverseFcn = [];
            iptSetPointerBehavior(app.region1, pm);
            % Enable pointer manager for app
            iptPointerManager(app.UIFigure,'enable');
            
            % Create the Arduino serial object
            %app.arduinoObj = serialport("COM6", 9600);
            %configureTerminator(app.arduinoObj,"CR/LF");
            %flush(app.arduinoObj);
            % 
            for i=1:8 
                app.message = readline(app.arduinoObj);
                disp(app.message)
            end
            
            function cursorPositionFeedback(app, hobj, inout)
% When inout is 'in', change hobj facecolor to green and update textbox.
% When inout is 'out' change hobj facecolor to red, and clear textbox.
% Check tag property of hobj to identify the object.
switch lower(inout)
    case 'in'
        facecolor = 'g';
        txt = 'Motor(s) vibrating';
        pointer = 'fleur';
        writeline(app.arduinoObj, "4&MOTOR_1_2&0!");
%         message = readline(app.arduinoObj);
%         disp(message)
    case 'out'
        facecolor = 'r';
        txt = 'No';
        pointer = 'crosshair';
        writeline(app.arduinoObj, "0&NO_MOTOR&0!");             %% THIS IS THE LINE THAT İS SUPPOSED TO STOP İT
%         message = readline(app.arduinoObj);
%         disp(message)
        
end
hobj.FaceColor = facecolor;
app.TextAreaEditField.Value = txt;
set(app.UIFigure, 'Pointer', pointer)
end  
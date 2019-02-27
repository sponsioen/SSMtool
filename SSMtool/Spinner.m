function Spinner(Figure,Message,Action,Position)

persistent jObj


switch lower(Action)
    
    case 'start'
        
        if ~isempty(jObj)
            jObj.stop;
            jObj.setBusyText(Message);
            [~, hContainer] = javacomponent(jObj.getComponent, Position, Figure);
            delete(hContainer);
        end
        
        try
            
            % R2010a and newer
            iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
            iconsSizeEnums = javaMethod('values',iconsClassName);
            SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
            jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32,Message);
            jObj.getComponent.setBackground(java.awt.Color(0.94, 0.94, 0.94));
            jObj.setPaintsWhenStopped(0);
            
        catch
            
            % R2009b and earlier
            redColor   = java.awt.Color(1,0,0);
            blackColor = java.awt.Color(0,0,0);
            jObj = com.mathworks.widgets.BusyAffordance(redColor, blackColor);
            
        end
        
        jObj.setPaintsWhenStopped(true);  % default = false
        jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
        javacomponent(jObj.getComponent, Position, Figure);
        jObj.start;
        drawnow();
        
    case 'update'
        
        jObj.setBusyText(Message);
        
    case 'stop'
        
        if isempty(jObj)
            return
        end
        
        jObj.stop;
        jObj.setBusyText(Message);
        pause(0.5);
        [~, hContainer] = javacomponent(jObj.getComponent, Position, Figure);
        delete(hContainer);
        
    otherwise
        error('Unknown spinner action : "%s".',Action);
        
end

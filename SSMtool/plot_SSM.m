function [output] = plot_SSM(handles)
nm = load('cs.mat');
SSM_function = @(z1,z2)eval(strcat('SSM_function_',nm.folder_id,'(z1,z2)'));
n_dof = numel(handles.sys.spv);
p = handles.plot.select;
modal =  handles.sys.modal;
T = cell(1,n_dof);

if modal 
    for i = 1:n_dof
       T{i} = sprintf('q_%d',i);
    end
else 
    for i = 1:n_dof
       T{i} = char(handles.sys.spv(i));
    end
end
    
handles.sys.Toutput = T;  
set(handles.popup_x_cor,'String',{T{:}});
set(handles.popup_x_cor,'Value',handles.plot.select(1));
set(handles.popup_y_cor,'String',{T{:}});
set(handles.popup_y_cor,'Value',handles.plot.select(2));
set(handles.popup_z_cor,'String',{T{:}});
set(handles.popup_z_cor,'Value',handles.plot.select(3));    
xmin = handles.plot.xmin;
xmax = handles.plot.xmax;
ymin = handles.plot.ymin; 
ymax = handles.plot.ymax;
zmin = handles.plot.zmin; 
zmax = handles.plot.zmax;   
r_max = handles.plot.r_max;    
rho = linspace(0,r_max,200); 
theta = linspace(0,2*pi,200); 
[Rho, Theta] = meshgrid(rho,theta);
S = T;
    
if handles.complex_cor
    [S{:}] = SSM_function(Rho.*exp(1i.*Theta),Rho.*exp(-1i.*Theta));
else
    [S{:}] = SSM_function(Rho.*cos(Theta),Rho.*sin(Theta));
end
  
X = real(S{p(1)});
Y = real(S{p(2)});
Z = real(S{p(3)});

if Z == 0
    Z = zeros(100);
end

if Y == 0
    Y = zeros(100);
end

if X == 0 
    X = zeros(100);
end

rho_mesh = linspace(0.1,r_max,5);
  
if ~isempty(handles.fig)
    if isvalid(handles.fig)
        clf(handles.fig)
        close(handles.fig)
    end
end
  
handles.fig = figure;
hold on
h = surfl(X,Y,Z);
handles.plot.h = h;
     
for i = 1:5 
    S_mesh = T;
    if handles.complex_cor
        [S_mesh{:}] = SSM_function(rho_mesh(i).*exp(1i.*theta),rho_mesh(i).*exp(-1i.*theta));
    else
        [S_mesh{:}] = SSM_function(rho_mesh(i).*cos(theta),rho_mesh(i).*sin(theta));
    end
    
    X_mesh = real(S_mesh{p(1)});
    Y_mesh = real(S_mesh{p(2)});
    Z_mesh = real(S_mesh{p(3)});

    if Z_mesh == 0
        Z_mesh = zeros(1,100);
    end

    if Y_mesh == 0
        Y_mesh = zeros(1,100);
    end

    if X_mesh == 0 
        X_mesh = zeros(1,100);
    end

    line(X_mesh,Y_mesh,Z_mesh,'Color',[180/255,180/255,180/255],'LineWidth',1.5)

end
 
grid on
handles.axis_handle = handles.fig.CurrentAxes;
set(handles.axis_handle,'FontSize',12)
handles.legend_handle = {'SSM'};
legend(handles.axis_handle,handles.legend_handle);
xlabel(T{p(1)})
ylabel(T{p(2)})
zlabel(T{p(3)})
axis auto
axis square
box on
h.EdgeColor= 'r';
shading interp
h.FaceColor=[255/255 165/255,0];
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.6; 
h.DiffuseStrength = 0.85;
h.SpecularStrength = 0.9;
h.SpecularExponent = 35;
h.BackFaceLighting = 'unlit';
view(-60,20);
camlight(0,90)
camlight(0,-90)
set(handles.axis_handle,'FontUnits','normalized');
output = handles;
    
end

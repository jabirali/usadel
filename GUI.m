function GUI()
control = figure;
control_p = uipanel('Title','Superconductor/Ferromagnet Bilayer','FontSize',16,'FontName','Arial','BackgroundColor',[0.6,0.6,1]);

                
%control_s_diffusion_text = text(control_s, 'Diffusion constant', [0.1 0.8 0.8 0.1]);
%control_s_diffusion_edit = edit(control_s, 1.0, [0.1 0.6 0.8 0.1]);

control_c   = panel(control_p,    'Entire system', [12 210 400-18 60]);
temperature = textedit(control_c, 'Temperature', 1.0,  [  12 20 180 20]);
debyeenergy = textedit(control_c, 'Debye energy', 1.0, [ 280 20 180 20]);

control_c   = panel(control_p,    'Superconductor', [12 140 400-18 60]);
% Diffusion constant
% Material length
% Superconductor strength (N₀λ)
% Interface parameter
%control_s   = panel(control_p, 'Superconductor', [0.02 1/3 0.96 1/3]);

% USE EVAL FOR GETTING NUMS!

% Exchange strength
% Exchange angle
% SOC strength
% SOC angle
% Diffusion constant
% Material length
% Interface parameter


% control_s = panel(control_p, 'Superconductor', [0.02 0.20 0.46 0.40]);
% control_s_diffusion = textedit(control_s, 'Diffusion constant', 1.0, [ 12 200 200 20]);
% control_s_diffusion = textedit(control_s, 'Material constant', 1.0, [ 12 150 200 20]);
% 
% control_f = panel(control_p, 'Ferromagnet', [0.52 0.2 0.46 0.75]);
% control_s_diffusion = textedit(control_f, 'Diffusion constant', 1.0, [ 12 200 200 20]);
% control_s_diffusion = textedit(control_f, 'Debye energy', 1.0, [ 12 175 200 20]);
% control_s_diffusion = textedit(control_f, 'SOC angle', 1.0, [ 12 150 200 20]);
% control_s_diffusion = textedit(control_f, 'Temperature', 1.0, [ 12 125 200 20]);

                
control_bu = uicontrol('Parent',control_p,'String','Update model', 'Position',[ 12 18 110 36]);
control_be = uicontrol('Parent',control_p,'String','Show solution', 'Position',[130 18 110 36]);
control_be = uicontrol('Parent',control_p,'String','Show gap', 'Position',[250 18 110 36]);
control_bx = uicontrol('Parent',control_p,'String','Show DoS', 'Position',[370 18 110 36]);
end

function handle = panel(parent, string, position)
	handle = uipanel('Parent',parent,                ...
                     'Title',string,                 ...
                     'FontSize',14,                  ...
                     'FontName','Arial',             ...
                     'BackgroundColor',[0.8,1,0.8],  ...
                     'Units', 'points',              ...
                     'Position', position);
end

function handle = text(parent, string, position)
    handle = uicontrol('Parent', parent, 'Style', 'text', 'String', string, 'Position', position, 'BackgroundColor', [0.8,1,0.8], 'HorizontalAlignment', 'Left');
end

function handle = edit(parent, string, position)
    handle = uicontrol('Parent', parent, 'Style', 'edit', 'String', string, 'Position', position, 'BackgroundColor', [1,1,1]);
end

function handle = textedit(parent, name, value, position)
    text(parent, name, [position(1) position(2) position(3)-60 position(4)]);
    handle = edit(parent, value, [position(1)+position(3)-55 position(2) 55 position(4)]);
end
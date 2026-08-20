model = 'Mdl_BianchiOL/WECS';   % replace with your model name
load_system(model)

% Find all Inport and Outport blocks
inports  = find_system(model, 'BlockType', 'Inport');
outports = find_system(model, 'BlockType', 'Outport');

% Combine them
ports = find_system(model); %[inports; outports];

% Set font size
for idx = 1: length(ports)
    set_param(ports{idx}, 'FontSize', '14');
end

%You can also set other font properties similarly:
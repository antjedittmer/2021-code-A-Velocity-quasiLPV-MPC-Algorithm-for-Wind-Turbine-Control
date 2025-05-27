% Define the structure fields for the WECS bus
fields = {
    'N', 1;
    'Ng', 1;
    'mh', 1;
    'mbl', 1;
    'mb', 1;
    'mn', 1;
    'mtower', 1;
    'mr', 1;
    'mt', 1;
    'mtb', 1;
    'H', 1;
    'Jg', 1;
    'Jr', 1;
    'JrL', 1;
    'Js', 1;
    'wnb', 1;
    'wnt', 1;
    'wntsw', 1;
    'zetat', 1;
    'zetab', 1;
    'Kt', 1;
    'Bt', 1;
    'Ktsw', 1;
    'Btsw', 1;
    'Kb', 1;
    'Bb', 1;
    'Ks', 1;
    'Bs', 1;
    'Rr', 1;
    'rb', 1;
    'etag', 1;
};

% Create bus elements
elems = Simulink.BusElement.empty;
for i = 1:size(fields, 1)
    elem = Simulink.BusElement;
    elem.Name = fields{i,1};
    elem.Dimensions = fields{i,2};
    elems(end+1) = elem; %#ok<SAGROW>
end

% Create the Bus object
wecsBus = Simulink.Bus;
wecsBus.Elements = elems;
wecsBus.Description = 'Bus for Wind Energy Conversion System (WECS)';

% Load or create the data dictionary
dictName = 'DD_Mdl1.sldd';
if ~isfile(dictName)
    Simulink.data.dictionary.create(dictName);
end
dictObj = Simulink.data.dictionary.open(dictName);
section = getSection(dictObj, 'Design Data');

% Add or update the bus object in the dictionary
entryName = 'wecsBus';
try
    entry = getEntry(section, entryName);
    setValue(entry, wecsBus);  % Update existing entry
catch
    addEntry(section, entryName, wecsBus);  % Add new entry
end

% Save and close the dictionary
saveChanges(dictObj);
close(dictObj);

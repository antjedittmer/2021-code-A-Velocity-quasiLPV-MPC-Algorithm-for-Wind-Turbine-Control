function restoreBaseVars(values)
% restoreBaseVars assigns the fields of the struct values to the base
% workspace (e.g. the old values returned by setFilteredWindQLMPC).

names = fieldnames(values);
for idx = 1:numel(names)
    assignin('base', names{idx}, values.(names{idx}));
end
end

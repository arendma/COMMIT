function generateCOMMGENInputFiles(model, outDir)

%% reactions.tsv
eq = regexprep(printRxnFormula(model, 'printFlag', 0), '\[c', '@MNXC1');
eq = regexprep(eq, '\[e', '@MNXC2');
eq = regexprep(eq, '\]', '');
eq = regexprep(eq, '\+\ MNX', '+ 1\ MNX');
eq = regexprep(eq, '[>-]\ MNX', '> 1\ MNX');
eq = regexprep(eq, '^MNX', '1\ MNX');
eq = regexprep(eq, '->', '-->');
eq = regexprep(eq, '<-', '<--');
eq = regexprep(eq, '<=>', '<==>');

reactions = cell2table([...
    strtok(model.rxns, '_'),...
    eq,...
    strtok(model.rxns, '_'),...
    cellfun(@(x)strcat('Protein', num2str(x)),  num2cell(1:numel(model.rxns)),...
    'UniformOutput', false)',...
    strtok(model.rxns, '_'),...
    model.EC,...
    cell(numel(model.rxns), 1),...
    cell(numel(model.rxns), 1),...
    num2cell(model.lb),...
    num2cell(model.ub),...
    model.rxnNames,...
    repmat({'--'}, numel(model.rxns), 1),...
    num2cell(zeros(numel(model.rxns), 1)),...
    num2cell(zeros(numel(model.rxns), 1))...
    ]);

writetable(reactions, fullfile(outDir, 'reactions.csv'),...
    'WriteVariableNames', false,...
    'Delimiter', '\t');
unix(['mv ',...
    fullfile(outDir, 'reactions.csv '),...
    fullfile(outDir, 'reactions.tsv')]);

%% chemicals.tsv
if isfield(model, 'metCharges')
    charges = num2cell(model.metCharges);
else
    charges = num2cell(zeros(numel(model.mets), 1));
end

species = cell2table([...
    strtok(model.mets, '['),...
    regexprep(model.metNames, '_[a-z0]+$', ''),...
    strtok(model.mets, '['),...
    model.metFormulas,...
    num2cell(zeros(numel(model.mets), 1)),...
    charges,...
    repmat({'--'}, numel(model.mets), 1),...
    ]);

writetable(species, fullfile(outDir, 'chemicals.csv'),...
    'WriteVariableNames', false,...
    'Delimiter', '\t');
unix(['mv ',...
    fullfile(outDir, 'chemicals.csv '),...
    fullfile(outDir, 'chemicals.tsv')]);

%% compartments.tsv

compartments = cell2table([...
    {'MNXC1'; 'MNXC2'},...
    {'Cytosol'; 'Extracellular'},...
    {'c'; 'e'}...
    ]);

writetable(compartments, fullfile(outDir, 'compartments.csv'),...
    'WriteVariableNames', false,...
    'Delimiter', '\t');
unix(['mv ',...
    fullfile(outDir, 'compartments.csv '),...
    fullfile(outDir, 'compartments.tsv')]);

%% model.tsv

model = cell2table([...
    cellstr(model.id);...
    cellstr(model.description);...
    num2cell(min(model.lb));...
    num2cell(max(model.ub));...
    {'2019/02/13'};...
    cellstr(date);...
    cellstr(model.id);...
    cellstr(model.id);...
    {'text/xml'}],...
    'RowNames',...
    {'ID', 'Description', 'LB', 'UB', 'MNXref Version', 'Processed Date',...
    'Source ID', 'Orig filename', 'MIME type'});
    
writetable(model, fullfile(outDir, 'model.csv'),...
    'WriteVariableNames', false,...
    'WriteRowNames', true, 'Delimiter', '\t');
unix(['mv ',...
    fullfile(outDir, 'model.csv '),...
    fullfile(outDir, 'model.tsv')]);

end



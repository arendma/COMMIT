function [trList, fzList] = translateIDs(idList, id_type, dbTable, source, target,  verbose, fztranslate)
%% trList = translateIDs(idList, id_type, dbTable, source, target, verbose)
% Translate metabolite IDs from one namespace to another (ModelSEED, KEGG,
% MetaCyc, BiGG, Rhea, EC, ChEBI).
% Input:
%       cellstr idList:             array containing the ids to be translated
%       char id_type:               either 'rxn' or 'met'
%       table dbTable:              (optional) table that contains the 
%                                   translation of identifiers for the given 
%                                   type (if empty attempts to load table
%                                   from 'COMMIT/data/tables/MNXref')
%       char source:                source namespace ('ModelSEED',
%                                   'KEGG', 'MetaCyc', 'BiGG', 'MNXref',
%                                   (for rxns: 'Rhea', 'EC')
%                                   (for mets: 'ChEBI', 'NAMES')
%       char target:                target namespace
%       logical verbose (optional): if true, print warnings and
%                                   progress statements (default: true)
%       logical fztranslate(optional):        if true, for unmatched IDs a second 
%                                   cell of ids in the
%                                   source database with the shortes
%                                   levensthein distance to model IDs is
%                                   returned
% Output:
%       cell trList:                array containing the translated IDs
%       cell fzList:              if fztranslate==true a cell of closest
%                                   fuzzy
%                                   


if nargin < 6 || ~islogical(verbose)
    verbose = true;
end
if nargin < 7 
    fztranslate = false;
end

if isempty(dbTable)
    dbDir = 'data/tables/MNXref';
    if isequal(id_type, 'met')
        dbTable = readtable(fullfile(dbDir, 'MNXref-met-translation-table.csv'),...
            'Format', '%s\t%s\t%s\t%s\t%s\t%s\t%s');
    else
        dbTable = readtable(fullfile(dbDir, 'MNXref-rxn-translation-table.csv'),...
            'Format', '%s\t%s\t%s\t%s\t%s\t%s\t%s');
    end
end

if ~iscellstr(idList)
    if ischar(idList)
        idList = {idList};
    elseif iscell(idList)
        idList = cellfun(@(x)char(x), idList, 'UniformOutput', false);
    else
        error('The input idList is either not given or not of type cellstr')
    end
end

if ~ischar(source) || ~ischar(target) || any(~ismember(id_type, ['rxn' 'met']))
    error('The source and/or target namespace or type definition is incorrect')
end

% Available namespaces
rxnSources = {'MNXref', 'KEGG', 'BiGG', 'MetaCyc', 'ModelSEED', 'Rhea', 'EC'};
metSources = {'MNXref', 'KEGG', 'BiGG', 'MetaCyc', 'ModelSEED', 'ChEBI', 'NAMES'};

%% Find matched using the MNXref database

% get the respective column that contains the desired names spaces
if isequal(id_type, 'met')
    sourceID = find(contains(metSources, source));
    targetID = find(contains(metSources, target));
else
    sourceID = find(contains(rxnSources, source));
    targetID = find(contains(rxnSources, target));
end

if isempty(sourceID) || isempty(targetID)
    error('The requested namespace is not available')
end

% initialize translated IDs list
trList = repmat({''}, numel(idList), 1);
% initialize fuzy translation IDs list
if fztranslate 
    fzList = repmat({''}, numel(idList), 1);
elseif nargout>1
    fzList=[];
end

% Filter source  and target IDs by empty keys and keys that are definitely not
% contained in the given list
sourceIDs = dbTable.(source);
targetIDs = dbTable.(target);
sourceIDs = sourceIDs(contains(dbTable.(source), idList));
targetIDs = targetIDs(contains(dbTable.(source), idList));
clear dbTable;

% pad source IDs and queries with '|'
sourceIDs = strcat('|', sourceIDs);
sourceIDs = strcat(sourceIDs, '|');
idList = strcat(idList, '|');
idList = strcat('|', idList);

% match IDs
match_idx = cellfun(@(x)find(contains(sourceIDs, x)), idList, 'UniformOutput', false);
empty_idx = cellfun('isempty', match_idx);
trList(~empty_idx) = targetIDs([match_idx{:}]);


if verbose
    fprintf('\nTranslated %3.2f%%\n\n', 100*sum(~empty_idx)/numel(idList));
end

%create fuzzy string matched list
if fztranslate
    fprintf('Fuzzy string matching for ids not foung in %s', source)
    %create matrix of levensthein distances
    fzmatch=zeros(length(sourceIDs), sum(empty_idx));
    empty_int_idx=find(empty_idx);
    for i=1:length(empty_int_idx)
        query=strrep(idList{empty_int_idx(i)}, '|', '');
        for j=1:length(sourceIDs)
            reference=strrep(sourceIDs{j}, '|', '');
            min_dist=fzsearch(reference,query, 0, 1);
            min_dist=min_dist{1};
            fzmatch(j,i)=min_dist(1);
        end
        min_dist=min(fzmatch(:,i));
        fzList{empty_int_idx(i)}=strjoin(sourceIDs(fzmatch(:,i)==min_dist),';');
        if mod(i, 100)==0
            fprintf("Processed %u ids\n", i)
        end
    end
end

% make it a column vector
trList = reshape(strtrim(trList), numel(trList), 1);

trList = strrep(trList, ';', '|');

end

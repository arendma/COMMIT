function model_conversion_iCre1355
%Script to convert the iCre1355 model MetaCyc format
rng('default')
verbose=true;
translationDB=loadTranslationDB;
pcyc_smartabpath='data/tables/smartTabs_Plantcyc';

mod=readCbModel('data/models/iCre1355/iCre1355_mixo_upd.xml');
mod=creategrRulesField(mod);
%import KEGG ids 
rxn_inf=readtable('data/models/iCre1355/tpj13059-sup-0011-tables5.xlsx');
%in 180 cases ids are not identical but the ordering of the table is
%identical to the ordering of rxns, make sure ordering overlaps 
if sum(~cellfun(@strcmp, rxn_inf.ReactionID, mod.rxns))<=180
    mod.rxnKEGGID=rxn_inf.KEGGID;
else
    error("The iCre1355 model rxn IDs are likely not ordered identical to the supplemental table rxn IDs")
end
met_inf=readtable('data/models/iCre1355/tpj13059-sup-0012-tables6.xlsx');
%same for the metabolite table 
if sum(~cellfun(@strcmp, met_inf.Abbreviation, mod.mets))<=836
    mod.metKEGGID=met_inf.KEGGCmpdID;
else
    error("The iCre1355 model met IDs are likely not ordered identical to the supplemental table met IDs")
end

%translate BIGG ids using COMMIT
[conv_mod, umRxns, umMet]=translateModel(mod, 'BiGG', 'MetaCyc', translationDB, 0, true, true);

%translate KEGG IDs
rxn_from_KEGG=repmat({''}, length(mod.rxns),1);
rxn_from_KEGG(~cellfun(@isempty, mod.rxnKEGGID))=translateIDs(mod.rxnKEGGID(~cellfun(@isempty, mod.rxnKEGGID)), 'rxn', translationDB.rxnTab, 'KEGG', 'MetaCyc', true, false);
met_from_KEGG=repmat({''}, length(mod.mets),1);
met_comps = regexp(mod.mets, '\[\w\d?\]', 'match');
met_comps = vertcat(met_comps{:});
met_from_KEGG(~cellfun(@isempty, mod.metKEGGID))=translateIDs(mod.metKEGGID(~cellfun(@isempty, mod.metKEGGID)), 'met', translationDB.metTab, 'KEGG', 'MetaCyc', true, false);
met_from_KEGG(~cellfun(@isempty, met_from_KEGG))=strcat(met_from_KEGG(~cellfun(@isempty, met_from_KEGG)), met_comps(~cellfun(@isempty,met_from_KEGG)));
%create a conversion table 
is_tr_rxns=ones(length(mod.rxns),1); %logical indicating if sucessfull translated
is_tr_mets=ones(length(mod.mets),1);
fz_match_rxns=cell(length(mod.rxns),1); %if not translated save fuzzy match
fz_match_mets=cell(length(mod.mets),1);
tr_rxns=conv_mod.rxns;%Array only keeping translated ids
tr_mets=conv_mod.mets;
tr_rxns_KEGG=rxn_from_KEGG; 
tr_mets_KEGG=met_from_KEGG;
for i=1:length(mod.rxns)
    if strcmp(conv_mod.rxns(i), mod.rxns(i))
        is_tr_rxns(i)=0;
        fz_match_rxns(i)=umRxns(ismember(umRxns(:,1), conv_mod.rxns(i)),2);
        tr_rxns{i}='';
    end
    if strcmp(tr_rxns_KEGG(i), mod.rxnKEGGID(i))
        tr_rxns_KEGG{i}='';
    end
end
for i=1:length(mod.mets)
    if strcmp(conv_mod.mets(i), mod.mets(i))
        is_tr_mets(i)=0;
        fz_match_mets(i)=umMet(ismember(umMet(:,1), conv_mod.mets(i)),2);
        tr_mets{i}='';
    end
    if strcmp(tr_mets_KEGG(i), mod.metKEGGID(i))
        tr_mets_KEGG{i}='';
    end
end

conv_rxn_out=table(mod.rxns, mod.grRules, tr_rxns, tr_rxns_KEGG, is_tr_rxns, fz_match_rxns, 'VariableNames', {'Orig_rxnID', 'GPR', 'COMMITtr_rxnID', 'COMMITtr_rxnID_fromKEGG' ,'is_tr', 'COMMITfz_rxnID'});
conv_mets_out=table(mod.mets, tr_mets, tr_mets_KEGG, is_tr_mets, fz_match_mets, 'VariableNames', {'Orig_metID', 'COMMITtr_metID', 'COMMITtr_rxnID_fromKEGG', 'is_tr', 'COMMITfz_metID'});
%save fuzzy matches
writetable(conv_rxn_out, 'data/models/iCre1355/iCre1355_mixo_rxn_conv.xlsx')
writetable(conv_mets_out, 'data/models/iCre1355/iCre1355_mixo_met_conv.xlsx')

%match reactions using model borgifier
%Generate a MetaCyc supermodel using RAVEN functions
McycMod=getModelFromMetaCyc([],true, true, true);
%removing reactions without non-zero entries in stochiometric matrix
%since calculation of network features fails otherwise
McycMod=removeRxns(McycMod, McycMod.rxns(sum(abs(McycMod.S),1)==0));

borg_mod=mod;
%using found metacyc IDs for metabolites to ease fitting
for i=1:length(borg_mod.mets)
    if ~isempty(conv_mets_out.COMMITtr_metID{i}) && ~contains(conv_mets_out.COMMITtr_metID{i}, '|')
        borg_mod.mets(i)=conv_mets_out.COMMITtr_metID(i);
    elseif ~isempty(conv_mets_out.COMMITtr_rxnID_fromKEGG{i}) && ~contains(conv_mets_out.COMMITtr_rxnID_fromKEGG{i}, '|')
        borg_mod.mets(i)=conv_mets_out.COMMITtr_rxnID_fromKEGG(i);
    end
end
%try fixing based on metabolites
rxns_fromMCYC_met=metbased_match(borg_mod, McycMod);
if isequal(conv_rxn_out.Orig_rxnID, rxns_fromMCYC_met.Orig_rxnID)
    conv_rxn_out=[conv_rxn_out rxns_fromMCYC_met(:,2:3)];
else
    error('Output of metbased_match does not overlap with the ordering for conv_rxn_out table. Have the reactions in borg_mod been modified?')
end
writetable(conv_rxn_out, 'data/models/iCre1355/iCre1355_mixo_rxn_conv.xlsx')

%fixing duplicates that can not be autonumbered in conv_mod
borg_mod.rxns{1857}='MALATE-DEHYDROGENASE-NADP+-RXN_2';
borg_mod=verifyModelBorg(borg_mod, 'Verbose', true);
McycMod=verifyModelBorg(McycMod, 'Verbose', true);


%build Tmodel
Tmod=buildTmodel(McycMod);

[borg_mod, Tmod, score, Stats] = compareCbModels(borg_mod, Tmod, 'Verbose', true);
%[bestMatch, bestMatchIdx]=max(Stats.scoreTotal, [], 1);
%save best match of commit untranslated reactions
% match metabolites

%create a rxn match vector from commit to train the borgify classifier
borg_rxnList=ones(length(borg_mod.rxns),1)*-1;

for i=1:length(borg_rxnList)
    %get entry in conv_table
    conv_idx=strcmpi(borg_mod.rxns{i}, conv_rxn_out.Orig_rxnID); 
        if sum(conv_idx)~=1 %check if borg cmodel ids is mappable to original model
            if sum(conv_idx)==0 && verbose
                warning([borg_mod.rxns{i} ' is not mappable back to original iCre1355 model. Skipping...'])
            elseif sum(conv_idx)>1 && verbose
                warning([borg_mod.rxns{i} ' maps to more than one ID in the original iCre1355 model. Skipping...'])
            end
        else
            %check if there is a commit tranlated ID available
            if ~isempty(conv_rxn_out.COMMITtr_rxnID{conv_idx}) %if ID from BIGG is available take this one
                query=conv_rxn_out.COMMITtr_rxnID(conv_idx);
            elseif ~isempty(conv_rxn_out.COMMITtr_rxnID_fromKEGG{conv_idx})
                query=conv_rxn_out.COMMITtr_rxnID_fromKEGG(conv_idx);
            elseif ~isempty(conv_rxn_out.met_fw_match{conv_idx})
                query=conv_rxn_out.met_fw_match(conv_idx); %only use forward matches from metabolites for training
            end
            if ~isempty(query)
                %find commit translated metacyc id in metacyc model
                query=strsplit(query{1}, '|');
                query=strrep(query,'-', '_');
                crxn_idx=[];
                j=1;
                while isempty(crxn_idx) && j<=length(query)
                    crxn_idx=find(strcmpi(query{j}, Tmod.rxns));
                    j=j+1;
                end
                if length(crxn_idx)==1
                    borg_rxnList(i)=crxn_idx;
                elseif verbose
                    if isempty(crxn_idx)
                        warning([strjoin(query, '|'), ' (', num2str(i), ') not mappable to metacyc template model.Skipping...'])
                        
                    elseif length(crxn_idx)>1
                        warning([strjoin(query, '|') ' maps to more than one reaction in  metacyc template model.Skipping...'])
                    end
                end
                query='';
            end
        end
end

%train classifier
global SCORE TMODEL CMODEL
SCORE=score;
TMODEL=Tmod;
CMODEL=borg_mod;

%Create a training rxnList with only 80% of the assignments
test_set=randsample(find(borg_rxnList~=-1), ceil(sum(borg_rxnList~=-1)*0.2));
train_rxnList=borg_rxnList;
train_rxnList(test_set)=-1;

% read int metacyc reaction lists 
pcyc_tabfiles={dir([pcyc_smartabpath, '/PCyC*.txt']).name};
pcyc_rxnIDs={};
for nm=pcyc_tabfiles
    tmp_ids=readtable(fullfile(pcyc_smartabpath, nm{:}));
    tmp_ids=tmp_ids.Reaction;
    pcyc_rxnIDs=union(pcyc_rxnIDs, tmp_ids);
end
iCre1355_mcyc_rxns=union(conv_rxn_out.COMMITtr_rxnID(~cellfun(@isempty, conv_rxn_out.COMMITtr_rxnID)), ...
    conv_rxn_out.COMMITtr_rxnID_fromKEGG(~cellfun(@isempty,conv_rxn_out.COMMITtr_rxnID_fromKEGG) & cellfun(@isempty, conv_rxn_out.COMMITtr_rxnID)));
%seperate ambigous mappped ids
iCre1355_mcyc_rxns=cellfun(@strsplit, iCre1355_mcyc_rxns, repmat({'|'}, length(iCre1355_mcyc_rxns), 1), 'UniformOutput', false);
iCre1355_mcyc_rxns=unique(horzcat(iCre1355_mcyc_rxns{:})');
fprintf('%3.0f Percent of metacyc reactions mapped from iCre1355 are in plantcyc and chlorophyta cyc union\n', sum(ismember(iCre1355_mcyc_rxns, pcyc_rxnIDs))/length(iCre1355_mcyc_rxns)*100)

%Declare meta net idx outside of plantcyc as negatives
neg_idx=[];
for i=1:length(Tmod.rxns)
    if ~any(strcmpi(Tmod.rxns{i}, strrep(pcyc_rxnIDs, '-', '_')))
        neg_idx=[neg_idx, i];
    end
end
train_RFStat=optimalScores(train_rxnList, 'RF');
train_linStats=optimalScores(train_rxnList, 'linear');

RFacc=sum(train_RFStat.bestMatchIndex(test_set)==borg_rxnList(test_set))/length(test_set);
linacc=sum(train_linStats.bestMatchIndex(test_set)==borg_rxnList(test_set))/length(test_set);
avgacc=sum(Stats.bestMatchIndex(test_set)==borg_rxnList(test_set))/length(test_set);
disp(['The accuracy of machine learning in the test set', num2str(length(test_set)), ' (20%) is ', num2str(linacc), ' (linear reg) ', ...
    num2str(RFacc), '(RandomForest) vs ', num2str(avgacc), 'for the average score without training'])
plot_borgperf(borg_rxnList, train_linStats, test_set, 'figures/borgifier/linear_perf_withmet')
plot_borgperf_pcycneg(borg_rxnList, train_linStats, test_set, 'figures/borgifier/linear_perf_withmet', neg_idx)

        
%choose the better performing model and train on all data
[~, mod_choice]=max([linacc, RFacc]);
learn_mod={'linear', 'RF'};
OptStats=optimalScores(borg_rxnList, learn_mod{mod_choice});




borg_bm_rxn=repmat({''}, length(mod.rxns),1); %best match rxn ID
borg_bm_rxn_score=nan(length(mod.rxns),1);%best match score
mcyc_newID=repmat({''}, length(mod.rxns),1); %liet of new ids automatically assigned
mcyc_newID_source=repmat({''}, length(mod.rxns),1);
%use conv_mod rxns for mapping to the conversion table, since
%con_rxn_out.COMMITtr_rxnId does not have all entries
for i=1:size(conv_rxn_out,1)
    if any(strcmpi(strrep(mod.rxns{i}, '-', '_'), borg_mod.rxns), 'all') %minus are replaced by underscores in borg models and all caps
        tmpidx=strcmpi(strrep(mod.rxns{i}, '-', '_'), borg_mod.rxns);
        borg_bm_rxn(i)=Tmod.rxns(OptStats.bestMatchIndex(tmpidx)); 
        borg_bm_rxn_score(i)=OptStats.bestMatch(tmpidx);
        if ~isempty(conv_rxn_out.COMMITtr_rxnID{i})
            mcyc_newID(i)=conv_rxn_out.COMMITtr_rxnID(i);
            mcyc_newID_source{i}='COMMIT_BIGG';
        elseif ~isempty(conv_rxn_out.COMMITtr_rxnID_fromKEGG{i})
            mcyc_newID(i)=conv_rxn_out.COMMITtr_rxnID_fromKEGG(i);
            mcyc_newID_source{i}='COMMIT_KEGG';
        elseif ~isempty(conv_rxn_out.met_fw_match{i})
            mcyc_newID(i)=conv_rxn_out.met_fw_match(i);
            mcyc_newID_source{i}='MET_FW';
        elseif ~isempty(conv_rxn_out.met_rev_match{i})
            mcyc_newID(i)=conv_rxn_out.met_rev_match(i);
            mcyc_newID_source{i}='MET_REV';
        elseif borg_bm_rxn_score(i)>0.65
            mcyc_newID(i)=strrep(upper(borg_bm_rxn(i)),'_', '-');
            mcyc_newID_source{i}='BORGIFIER';
        end
    end
end

%report
%Save a table with conversion suggestions
tmp_tab=table(borg_bm_rxn, borg_bm_rxn_score, mcyc_newID, mcyc_newID_source, 'VariableNames', {'Borg_best_match', 'Norm_Borg_score', 'Mcyc_ID', 'Mcyc_ID_source'});
conv_rxn_out=[conv_rxn_out, tmp_tab];
writetable(conv_rxn_out, 'data/models/iCre1355/iCre1355_mixo_rxn_conv.xlsx')


%sum(logical(~cellfun(@isempty,regexp(conv_rxn_out.Orig_rxnID, '^EX|^DM|Biomass','ignorecase')) + ~cellfun(@isempty, conv_rxn_out.COMMITtr_rxnID) + ~cellfun(@isempty, conv_rxn_out.COMMITtr_rxnID_fromKEGG)))/size(conv_rxn_out,1)

%use borigifier to map metabolites - without machine learning

%generate a model with original metabolit ids to come up with suggestions
%for the missing ones
origmet_borg_mod=mod;
origmet_borg_mod.rxns{1857}='MALATE-DEHYDROGENASE-NADP+-RXN_2';
origmet_borg_mod=verifyModelBorg(origmet_borg_mod, 'Verbose', true);

metScores=compareAllMets(origmet_borg_mod, Tmod);
[metBestMatch, metBestMatchIdx] = max(metScores, [], 2);
%save results
save('borgout.mat', 'borg_mod', 'Tmod', 'score', 'Stats','OptStats', 'train_rxnList', 'metBestMatch', 'metBestMatchIdx')

borg_bm_met=repmat({''}, length(mod.mets),1); %best match met ID
borg_bm_met_score=nan(length(mod.mets),1); %best match score
%use conv_mod mets for mapping to the conversion table, since
%con_rxn_out.COMMITtr_metId does not have all entries
for i=1:size(conv_mets_out,1)
    if any(strcmpi(strrep(mod.mets{i}, '-', '_'), origmet_borg_mod.mets), 'all') %minus are replaced by underscores in borg models and all caps
        tmpidx=strcmpi(strrep(mod.mets{i}, '-', '_'), origmet_borg_mod.mets);
        if metBestMatch(tmpidx)>0 %only assign matches for score higher 0 
            borg_bm_met(i)=Tmod.mets(metBestMatchIdx(tmpidx)); 
            borg_bm_met_score(i)=metBestMatch(tmpidx);
        end
    end
end

%Save a table with conversion suggestions
tmp_tab=table(borg_bm_met, borg_bm_met_score, 'VariableNames', {'Borg_best_match', 'Norm_Borg_score'});
conv_mets_out=[conv_mets_out, tmp_tab];
writetable(conv_mets_out, 'data/models/iCre1355/iCre1355_mixo_met_conv.xlsx')

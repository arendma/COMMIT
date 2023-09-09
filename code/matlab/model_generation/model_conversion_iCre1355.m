function model_conversion_iCre1355
%Script to convert the iCre1355 model MetaCyc format
verbose=true;
translationDB=loadTranslationDB;
mod=readCbModel('data/models/iCre1355/iCre1355_mixo_upd.xml');
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
met_from_KEGG(~cellfun(@isempty, mod.metKEGGID))=translateIDs(mod.metKEGGID(~cellfun(@isempty, mod.metKEGGID)), 'met', translationDB.metTab, 'KEGG', 'MetaCyc', true, false);

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

conv_rxn_out=table(mod.rxns, tr_rxns, tr_rxns_KEGG, is_tr_rxns, fz_match_rxns, 'VariableNames', {'Orig_rxnID', 'COMMITtr_rxnID', 'COMMITtr_rxnID_fromKEGG' ,'is_tr', 'COMMITfz_rxnID'});
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
%fixing duplicates that can not be autonumbered in conv_mod
borg_mod=mod;
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

train_RFStat=optimalScores(train_rxnList, 'RF');
train_linStats=optimalScores(train_rxnList, 'linear');

RFacc=sum(train_RFStat.bestMatchIndex(test_set)==borg_rxnList(test_set))/length(test_set);
linacc=sum(train_linStats.bestMatchIndex(test_set)==borg_rxnList(test_set))/length(test_set);
avgacc=sum(Stats.bestMatchIndex(test_set)==borg_rxnList(test_set))/length(test_set);
disp(['The accuracy of machine learning in the test set', num2str(length(test_set)), ' (20%) is ', num2str(linacc), ' (linear reg) ', ...
    num2str(RFacc), '(RandomForest) vs ', num2str(avgacc), 'for the average score without training'])
        
%choose the better performing model and train on all data
[~, mod_choice]=max([linacc, RFacc]);
learn_mod={'linear', 'RF'};
OptStats=optimalScores(borg_rxnList, learn_mod{mod_choice});


metScores=compareAllMets(borg_mod, Tmod);
[metBestMatch, metBestMatchIdx] = max(metScores, [], 2);
%save results
save('borgout.mat', 'borg_mod', 'Tmod', 'score', 'Stats','OptStats', 'train_rxnList', 'metBestMatch', 'metBestMatchIdx')

borg_bm_rxn=repmat({''}, length(mod.rxns),1); %best match rxn ID
borg_bm_rxn_score=nan(length(mod.rxns),1); %best match score
%use conv_mod rxns for mapping to the conversion table, since
%con_rxn_out.COMMITtr_rxnId does not have all entries
for i=1:size(conv_rxn_out,1)
    if any(strcmpi(strrep(conv_mod.rxns{i}, '-', '_'), borg_mod.rxns), 'all') %minus are replaced by underscores in borg models and all caps
        tmpidx=strcmpi(strrep(conv_mod.rxns{i}, '-', '_'), borg_mod.rxns);
        borg_bm_rxn(i)=Tmod.rxns(OptStats.bestMatchIndex(tmpidx)); 
        borg_bm_rxn_score(i)=OptStats.bestMatch(tmpidx);
    end
end

%Save a table with conversion suggestions
tmp_tab=table(borg_bm_rxn, borg_bm_rxn_score, 'VariableNames', {'Borg_best_match', 'Norm_Borg_score'});
conv_rxn_out=[conv_rxn_out, tmp_tab];
writetable(conv_rxn_out, 'data/models/iCre1355/iCre1355_mixo_rxn_conv.xlsx')


borg_bm_met=repmat({''}, length(mod.mets),1); %best match met ID
borg_bm_met_score=nan(length(mod.mets),1); %best match score
%use conv_mod mets for mapping to the conversion table, since
%con_rxn_out.COMMITtr_metId does not have all entries
for i=1:size(conv_mets_out,1)
    if any(strcmpi(strrep(conv_mod.mets{i}, '-', '_'), borg_mod.mets), 'all') %minus are replaced by underscores in borg models and all caps
        tmpidx=strcmpi(strrep(conv_mod.mets{i}, '-', '_'), borg_mod.mets);
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

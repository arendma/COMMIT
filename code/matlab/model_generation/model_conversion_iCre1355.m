function model_conversion_iCre1355
%Script to convert the iCre1355 model MetaCyc format
translationDB=loadTranslationDB;
mod=readCbModel('data/models/iCre1355/iCre1355_mixo_upd.xml');

%translate BIGG ids using COMMIT
[conv_mod, umRxns, umMet]=translateModel(mod, 'BiGG', 'MetaCyc', translationDB, 0, true, true);

%create a conversion table 
is_tr_rxns=ones(length(mod.rxns),1); %logical indicating if sucessfull translated
is_tr_mets=ones(length(mod.mets),1);
fz_match_rxns=cell(length(mod.rxns),1); %if not translated save fuzzy match
fz_match_mets=cell(length(mod.mets),1);
tr_rxns=conv_mod.rxns;%Array only keeping translated ids
tr_mets=conv_mod.mets;
for i=1:length(mod.rxns)
    if strcmp(conv_mod.rxns(i), mod.rxns(i))
        is_tr_rxns(i)=0;
        fz_match_rxns(i)=umRxns(ismember(umRxns(:,1), conv_mod.rxns(i)),2);
        tr_rxns{i}='';
    end
end
for i=1:length(mod.mets)
    if strcmp(conv_mod.mets(i), mod.mets(i))
        is_tr_mets(i)=0;
        fz_match_mets(i)=umMet(ismember(umMet(:,1), conv_mod.mets(i)),2);
        tr_mets{i}='';
    end
end

conv_rxn_out=table(mod.rxns, tr_rxns, is_tr_rxns, fz_match_rxns, 'VariableNames', {'Orig_rxnID', 'COMMITtr_rxnID', 'is_tr', 'COMMITfz_rxnID'});
conv_mets_out=table(mod.mets, tr_mets, is_tr_mets, fz_match_mets, 'VariableNames', {'Orig_metID', 'COMMITtr_metID', 'is_tr', 'COMMITfz_metID'});
%save fuzzy matches
writetable(conv_rxn_out, 'data/models/iCre1355/iCre1355_mixo_rxn_conv.xlsx')
writetable(conv_mets_out, 'data/models/iCre1355/iCre1355_mixo_met_conv.xlsx')

%match reactions using model borgifier
%Generate a MetaCyc supermodel using RAVEN functions
McycMod=getModelFromMetaCyc([],false,false,false);
%removing reactions without non-zero entries in stochiometric matrix
%since calculation of network features fails otherwise
McycMod=removeRxns(McycMod, McycMod.rxns(sum(abs(McycMod.S),1)==0));
%fixing duplicates that can not be autonumbered in conv_mod
borg_conv_mod=conv_mod;
borg_conv_mod.rxns{1857}='MALATE-DEHYDROGENASE-NADP+-RXN_2';
borg_conv_mod=verifyModelBorg(borg_conv_mod, 'Verbose', true);
McycMod=verifyModelBorg(McycMod, 'Verbose', true);


%build Tmodel
Tmod=buildTmodel(McycMod);

[borg_conv_mod, Tmod, score, Stats] = compareCbModels(borg_conv_mod, Tmod, 'Verbose', true);
%[bestMatch, bestMatchIdx]=max(Stats.scoreTotal, [], 1);
%save best match of commit untranslated reactions
borg_bm_rxn=repmat({''}, length(mod.rxns),1); %best match rxn ID
borg_bm_rxn_score=nan(length(mod.rxns),1); %best match score
%use conv_mod rxns for mapping to the conversion table, since
%con_rxn_out.COMMITtr_rxnId does not have all entries
for i=1:size(conv_rxn_out,1)
    if any(strcmpi(strrep(conv_mod.rxns{i}, '-', '_'), borg_conv_mod.rxns), 'all') %minus are replaced by underscores in borg models and all caps
        tmpidx=strcmpi(strrep(conv_mod.rxns{i}, '-', '_'), borg_conv_mod.rxns);
        borg_bm_rxn(i)=Tmod.rxns(Stats.bestMatchIndex(tmpidx)); 
        borg_bm_rxn_score(i)=Stats.bestMatch(tmpidx);
    end
end

%Save a table with conversion suggestions
tmp_tab=table(borg_bm_rxn, borg_bm_rxn_score, 'VariableNames', {'Borg_best_match', 'Norm_Borg_score'});
conv_rxn_out=[conv_rxn_out, tmp_tab];
writetable(conv_rxn_out, 'data/models/iCre1355/iCre1355_mixo_rxn_conv.xlsx')
% match metabolites
metScores=compareAllMets(borg_conv_mod, Tmod);
[metBestMatch, metBestMatchIdx] = max(metScores, [], 2);

borg_bm_met=repmat({''}, length(mod.mets),1); %best match met ID
borg_bm_met_score=nan(length(mod.mets),1); %best match score
%use conv_mod mets for mapping to the conversion table, since
%con_rxn_out.COMMITtr_metId does not have all entries
for i=1:size(conv_mets_out,1)
    if any(strcmpi(strrep(conv_mod.mets{i}, '-', '_'), borg_conv_mod.mets), 'all') %minus are replaced by underscores in borg models and all caps
        tmpidx=strcmpi(strrep(conv_mod.mets{i}, '-', '_'), borg_conv_mod.mets);
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

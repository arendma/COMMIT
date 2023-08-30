function model_conversion_iCre1355
%Script to convert the iCre1355 model MetaCyc format
translationDB=loadTranslationDB;
mod=readCbModel('data/models/iCre1355/iCre1355_mixo_upd.xml');

%translate BIGG ids using COMMIT
[conv_mod, umRxns, umMet]=translateModel(mod, 'BiGG', 'MetaCyc', translationDB, 0, true, true);
%save fuzzy matches
writecell(umRxns, 'data/models/iCre1355/iCre1355_mixo_unmatched_rxns.xlsx')
writecell(umMet, 'data/models/iCre1355/iCre1355_mixo_unmatched_mets.xlsx')

%match reactions using model borgifier
%Generate a MetaCyc supermodel using RAVEN functions
McycMod=getModelFromMetaCyc([],false,false,false);
%removing reactions without non-zero entries in stochiometric matrix
%since calculation of network features fails otherwise
McycMod=removeRxns(McycMod, McycMod.rxns(sum(abs(McycMod.S),1)==0));
%fixing duplicates that can not be autonumbered in conv_mod
conv_mod.rxns{1857}='MALATE-DEHYDROGENASE-NADP+-RXN_2';
conv_mod=verifyModelBorg(conv_mod, 'Verbose', true);
McycMod=verifyModelBorg(McycMod, 'Verbose', true);


%build Tmodel
Tmod=buildTmodel(conv_mod);

[McycMod, Tmod, score, Stats] = compareCbModels(McycMod, Tmod, 'Verbose', true);

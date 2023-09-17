function res = metbased_match(new_mod, template_mod)
%remove location info from template mod
map_mod=new_mod;
map_mod.mets=regexprep(map_mod.mets, '\[\w\d?\]$','');
%remove metabolites not in the list but keep reactions by summing over
%removed rows
i=1;
while i<=length(map_mod.mets)
    ident_mets=find(strcmp(map_mod.mets{i}, map_mod.mets));
    if length(ident_mets)>1
        map_mod.S(ident_mets(1),:)=sum(map_mod.S(ident_mets,:));
        map_mod=rm_met(map_mod, ident_mets(2:end));
    end
    i=i+1;
end
%remove empty reactions
%remove metabolites not in the templatemod and remove reactions that
%contain these
map_mod=removeMetabolites(map_mod, map_mod.mets(~ismember(map_mod.mets, template_mod.mets)), true, 'inclusive');
nm_comp_rxns=map_mod.rxns(sum(abs(map_mod.S),1)>0);
nm_comp_S=map_mod.S(:,sum(abs(map_mod.S),1)>0);
template_mod=removeMetabolites(template_mod, template_mod.mets(~ismember(template_mod.mets, map_mod.mets)));
[~, match_idx]=ismember(map_mod.mets,template_mod.mets );
tm_comp_S=template_mod.S(match_idx(match_idx>0),:);
tm_comp_rxns=template_mod.rxns(sum(abs(tm_comp_S),1)>0);
tm_comp_S=tm_comp_S(:, sum(abs(tm_comp_S),1)>0);
cor_mat=corr(nm_comp_S,full(tm_comp_S));

%generate result table
fw=cell(length(new_mod.rxns),1);
rev=cell(length(new_mod.rxns),1);
for i=1:length(fw)
    if ismember(new_mod.rxns(i), nm_comp_rxns)
        fw{i}=strjoin(tm_comp_rxns(cor_mat(ismember(nm_comp_rxns, new_mod.rxns(i)),:)==1), '|');
        rev{i}=strjoin(tm_comp_rxns(cor_mat(ismember(nm_comp_rxns, new_mod.rxns(i)),:)==-1), '|');
    end
end
res=table(new_mod.rxns, fw, rev, 'VariableNames',{'Orig_rxnID','met_fw_match', 'met_rev_match'});
end
function mod = rm_met(mod, met_idx)
    %Function to remove metabolites from a model and the linked fields
    fn=fieldnames(mod);
      mod.S(met_idx,:)=[];
        for f=[fn(contains(fn, 'met'))', 'b', 'csense']
            mod.(f{1})(met_idx)=[];
        end
end
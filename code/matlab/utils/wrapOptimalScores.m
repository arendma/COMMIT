function Stats = wrapOptimalScores(rxnList, optimizer,CMODEL, TMODEL, SCORE)
%Wrapper function for the optimalScores function from COBRA borgifier tool
%to avoid creating global variables in calling function/script
%INPUT:
%  rxnList:       reactions` list, an integer vector of length
%                     CMODEL.rxns, that contains the matched index of the 
%                     TMODEL for each reaction, -1 if the reaction has unknown matches
%                     0 if the reaction is declared new and not contained in
%                     the template model
%      optimizer:     Either 'svm', 'RF', 'linear', or 'exp'; the latter two are custom
%                     functions contained within the functions `optWeightLin` and
%                     optWeightExp.
%      CMODEL:        global input
%      TMODEL:        global input
%      SCORE:         global input
global SCORE CMODEL TMODEL
Stats=optimalScores(rxnList, optimizer)
end
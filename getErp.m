function [erp, sigCell] = getErp(EEG, tlength, Fs)
% This function takes the original EEG signal, time lengths of each trial 
% and sampling frequency as input, and output the event related potentail 
% and the orginal signal in cell format
nLength = round(tlength/Fs);
nRepeat = length(EEG(1,:))/nLength;
erp = zeros(1,nLength);
for i = 1 : nRepeat
    sigCell{i} = EEG(:,(i-1)*nLength+1:i*nLength);
    erp = erp + EEG(:,(i-1)*nLength+1:i*nLength);
end
erp = erp / nRepeat;
end
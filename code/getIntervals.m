function [T0, T1, T2] = getIntervals(recName,annotator,fs,sigLen)
    T0=[];
    T1=[];
    T2=[];
    [annot, ~, ~, ~, ~, cmt] = rdann(recName, annotator);
    cmt = string(cmt);
    for i=1:size(annot,1)
        splt=split(cmt(i), ' ');
        dur = str2double(splt(3))*fs;
        intEnd = annot(i) + dur-1;
        if (i<length(annot))
            intEnd = min(annot(i+1)-1, intEnd);
        else
            intEnd = min(sigLen, intEnd);
        end
        if (splt(1)=="T0")
            T0 = [T0; [annot(i) intEnd]];
        elseif (splt(1)=="T1")
            T1 = [T1; annot(i) intEnd];
        elseif (splt(1)=="T2")
            T2 = [T2; annot(i) intEnd];
        end
    end
end
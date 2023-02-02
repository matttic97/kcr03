function project03(recNr, subject, dir)
    % Pridobi bloke signalov
    record = getRecord(recNr);
    [T1, t1s, T2, t2s, fs] = getSignalBlocks(record, subject, dir, 4);

    % Povprečje začetnih blokov
    interval = 3;
    t1sm = getLearnMean(t1s, T1, interval);
    t2sm = getLearnMean(t2s, T2, interval);

    % DEFINICIJA KONSTANT ZA FILTRACIJO
    F = [0 8 8 13 13 fs/2] / (fs/2);
    A = [0 0 1 1 0 0];
    n=75; %sprobaj različne vrednsoti
    b = firls(n, F, A);

    % Metoda skupnih prostorskih vzorcev
    [W] = f_CSP(t1sm, t2sm);

    % IZRACUN ZNACILK
    lvt1 = calcFeatures(W, b, t1s);
    lvt2 = calcFeatures(W, b, t2s);

    figure;
    scatter(lvt1(:, 1), lvt1(:, 2));
    hold on;
    scatter(lvt2(:, 1), lvt2(:, 2));

    % ZAPIS ZNACILK V DATOTEKE
    featVFile = strcat(subject, 'featureVectors.txt');
    classFile = strcat(subject, 'referenceClass.txt');
    
    fvf = fopen(featVFile, 'wt');
    rcf = fopen(classFile, 'wt');
    for i=1:size(lvt1,1)
        fprintf(fvf, '%.8f %.8f\n', lvt1(i,1), lvt1(i,2));
        fprintf(rcf, 'T1\n');
    end
    for i=1:size(lvt2,1)
        fprintf(fvf, '%.8f %.8f\n', lvt2(i,1), lvt2(i,2));
        fprintf(rcf, 'T2\n');
    end
    fclose(fvf);
    fclose(rcf);

end

function [record] = getRecord(recNr)
    record=[];
    for i=0:2
        record = [record, string(num2str(recNr+(i*4),'%02d'))]
    end
end

function [T1, t1s, T2, t2s, fs] = getSignalBlocks(record, subject, dir, blockTime)
    subject=string(subject);
    t1s = {};
    t2s = {};
    
    dirName = "";
    nrCd = 0;
    for i=1:size(record,2)
        if (nargin<3) % 
            recName = strcat('/',dir,'/',subject,'/',subject,'R',record(i),'.edf');
        elseif (strcmp(dir,"")==1)
            recName = strcat(subject,'R', record(i), '.edf');
        else
            dirName=strcat(dir,'/',subject,'/');
            nrCd=count(dirName,'/');
            cd (dirName);
            recName=strcat(subject, 'R', record(i), '.edf');
        end
        recName=convertStringsToChars(recName);
        disp(recName)
        [sig, fs, ~] = rdsamp(recName, 1:64);
        size(sig)
        [~, T1, T2] = getIntervals(recName, 'event', fs, size(sig,1));
        sig = sig.';

        interval = fs * blockTime;
        [~,n] = size(sig);        
        for j=1:size(T1, 1)
          if (T1(j,1) + interval) > n
              t1s{end+1}=sig(:, T1(j,1):n);
          else
              t1s{end+1}=sig(:, T1(j,1):T1(j,1) + interval);
          end
        end    
        for j=1:size(T2,1)
          if (T2(j,1) + interval) > n
              t2s{end+1}=sig(:, T2(j,1):n);
          else
              t2s{end+1}=sig(:, T2(j,1):T2(j,1) + interval);
          end
        end      

        for i=1:nrCd
            cd('..');
        end
    end
end

function [ts_mean] = getLearnMean(ts, T, interval)
    ts_mean = cell2mat(ts(1));
    for j = 2: interval
      ts_mean = ts_mean + cell2mat(ts(j));
    end
    ts_mean = ts_mean / size(T,2);
end

function [features] = calcFeatures(W, b, ts)
    features = [];
    for i=1:size(ts,2)
        tmp = W*cell2mat(ts(i));
        tmp = [tmp(1,:).' tmp(64,:).'].';
        tmp = filter(b,1,tmp);
        features = [features; log(var(tmp(1,:))) log(var(tmp(2,:)))];
    end
end
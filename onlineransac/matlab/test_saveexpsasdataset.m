% get all scenarios data to save them

clear;
close all;

folderdataset = './datamanagementplan/';

expcat = ExperimentCatalog(1);
numexps = length(expcat);

for f = 1:numexps
    fprintf("Preparing experiment #%d...\n",f);
    [~,~,data] = ExperimentGet(expcat,expcat{f}.class,expcat{f}.index,1,Inf,0,NaN,0);
    
    datafile = sprintf("%styrell_roundtrip%03d.csv",folderdataset,f);

    [fileID,message] = fopen(datafile,'w');
    if fileID < 0
        error(message);
    end
    fclose(fileID);
    writematrix(data,datafile);

    metadatafile = sprintf("%styrell_roundtrip%03d_metadata.txt",folderdataset,f);
    [fileID,message] = fopen(metadatafile,'w');
    if fileID < 0
        error(message);
    end
    fclose(fileID);
    txt = {};
    txt{end+1} = convertStringsToChars(sprintf("Experiment ID: %03d",f));
    txt{end+1} = convertStringsToChars(sprintf("Name: %s",expcat{f}.name));
    txt{end+1} = convertStringsToChars(sprintf("Class: %s",expcat{f}.class));
    txt{end+1} = convertStringsToChars(sprintf("Time units: ms"));
    txt{end+1} = convertStringsToChars(sprintf("Density (bytes): %d",expcat{f}.dens));
    txt{end+1} = convertStringsToChars(sprintf("Network: %s",expcat{f}.netw));
    txt{end+1} = convertStringsToChars(sprintf("Geoloc: %s",expcat{f}.dist));
    txt{end+1} = convertStringsToChars(sprintf("Software: %s (robot), %s (station)",expcat{f}.serversw,expcat{f}.clientsw));
    writelines(txt,metadatafile);
    
end
function snirfObj = rearrange_channels(snirfObj)
intensity = snirfObj.GetData();
% dod = hmrR_Intensity2OD(snirfObj.data);
file = '2020-12-25_006_probe_correctSD.snirf';
snirf = SnirfClass(['data',filesep,file]);
measurementlist = snirf.data.measurementList;
new_dataTimeSeries = zeros(size(intensity.dataTimeSeries,1),length(measurementlist));
for i = 1:1:length(measurementlist)
    channel = measurementlist(i);
    sourceId = channel.sourceIndex;
    detectorId = channel.detectorIndex;
    waveId = channel.wavelengthIndex;
    % now we have source detector and wave, we should find the same in dod
    d_measurementlist = intensity.measurementList;
    for j = 1:1:length(d_measurementlist)
        d_channel = d_measurementlist(j);
        d_sourceId = d_channel.sourceIndex;
        d_detectorId = d_channel.detectorIndex;
        d_waveId = d_channel.wavelengthIndex;
        if sourceId == d_sourceId && detectorId == d_detectorId && waveId == d_waveId
            data = intensity.dataTimeSeries(:,j);
            break
        end
    end
    new_dataTimeSeries(:,i) = data;
end
snirfObj.data.SetDataTimeSeries(new_dataTimeSeries)
snirfObj.data.measurementList = measurementlist;
snirfObj.probe = snirf.probe;
end
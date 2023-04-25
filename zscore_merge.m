%z-socring across conditions
%data in cell format
%EEG data with format of timepoints x channels
function dOut = zscore_merge(dataIn)
tempOut = [];
[nsubj,ncon] = size(dataIn);
temp = dataIn{1,1};
ndim = ndims(temp);
dmerge = cat(ndim+1,dataIn{:});
%dmerge should be 3D matrix
for ele = 1:size(temp,2)
    tempd = squeeze(dmerge(:,ele,:));
    tempmerge = tempd(:);
    mm = nanmean(tempmerge,1);
    sd = std(tempmerge,1);
    for nc = 1:ncon
        tempOut(:,ele,nc) =  (dmerge(:,ele,nc)-mm)/sd;
    end
end
for nc = 1:ncon
    dOut{nc} = squeeze(tempOut(:,:,nc));
end
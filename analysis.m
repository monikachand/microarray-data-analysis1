gseData = getgeodata('GSE5847', 'ToFile', 'GSE5847.txt')
gseData = geoseriesread('GSE5847.txt');
get(gseData.Data)
gseData.Data(1:10, 1:10)
gseData.Header
gseData.Header.Series
gseData.Header.Samples
gseData.Header.Samples.data_processing(1)
sampleSources = unique(gseData.Header.Samples.source_name_ch1);
sampleSources{:}
gseData.Header.Samples.characteristics_ch1(:,1)
sampleGrp = gseData.Header.Samples.characteristics_ch1(1,:);
gplData = geosoftread('GPL96.txt');
gplData.ColumnNames
gplProbesetIDs = gplData.Data(:, strcmp(gplData.ColumnNames, 'ID'));
geneSymbols = gplData.Data(:, strcmp(gplData.ColumnNames, 'Gene Symbol'));
gseData.Data = rownames(gseData.Data, ':', geneSymbols);
gseData.Data(1:10, 1:10)
stromaIdx = strcmpi(sampleSources{1}, gseData.Header.Samples.source_name_ch1);
nStroma = sum(stromaIdx)
stromaData = gseData.Data(:, stromaIdx);
stromaGrp = sampleGrp(stromaIdx);
nStromaIBC = sum(strcmp('diagnosis: IBC',stromaGrp))
nStromaNonIBC = sum(strcmp('diagnosis: non-IBC',stromaGrp))
stromaData = colnames(stromaData,':',stromaGrp);
[mask, stromaData] = genevarfilter(stromaData);
stromaData.NRows
rng default
IBC = stromaData(:, 'diagnosis: IBC');
NonIBC= stromaData(:, 'diagnosis: non-IBC');
[pvalues, tscores] = mattest(IBC,NonIBC,...
                           'Showhist', true', 'showplot', true, 'permute', 1000);
cutoff = 0.01;
sum(pvalues < cutoff)  
figure;
[pFDR, qvalues] = mafdr(pvalues, 'showplot', true)
sum(qvalues < cutoff)
pvaluesBH = mafdr(pvalues, 'BHFDR', true);
sum(pvaluesBH < cutoff)
testResults = [tscores pvalues pFDR qvalues pvaluesBH];
testResults = colnames(testResults, 5, {'FDR_BH'});
testResults = sortrows(testResults, 2);
testResults(1:20, :)
diffStruct = mavolcanoplot(NonIBC,IBC, pvalues)
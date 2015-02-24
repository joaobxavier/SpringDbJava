function [mutationMatrix, mutationList] = getMutationMatrix(db, geneList, referenceStrain, useMuscle)

% parse inputs
% if no genome_id is provided as reference strain then use the consensus
% seuqence
if nargin == 2
    referenceStrain = 'consensus';
    useMuscle = false;
end

genome_id = db.getPhenotypeColumn('genome_id');


% % fetch sequences from database
% for i = 1:length(geneList)
%     seqsQuery = cell(db.getSequencesOfGene(geneList(i)));
%     for j = 1:length(seqsQuery)
%         aaseq = nt2aa(seqsQuery{j});
%         seqs(j).Sequence = aaseq(1:end-1); % must translate
%         seqs(j).Header = num2str(genome_id(j));
%     end
%     geneListSeqs(i).seqs = seqs;
% end

% get all sequences in a single query
seqsQuery = cell(db.getSequencesOfGenes(geneList));
for i = 1:size(seqsQuery, 2) % iterate over gene list
    for j = 1:size(seqsQuery, 1) % iterate over geenoe list
        aaseq = nt2aa(seqsQuery{j, i});
        seqs(j).Sequence = aaseq(1:end-1); % must translate
        seqs(j).Header = num2str(genome_id(j));
    end
    geneListSeqs(i).seqs = seqs;
end
    




% if a single gene is provided rather than a  list of genes then
% put gene in a list with a single element, becuase the rest of the
% algorithm requires a list
if ischar(geneList)
    geneList = {geneList};
end


% get mutations for each gene
mutationMatrix = [];
mutationList = {};


mutationData = [];
for i = 1:length(geneList)
    [mm, ml] = getMutationMatrixAux(geneListSeqs(i).seqs, geneList{i}, referenceStrain, useMuscle, genome_id);
    mutationData(i).mutationMatrix = mm;
    mutationData(i).mutationList = ml;
end

for i = 1:length(geneList)
    mutationMatrix = [mutationMatrix, mutationData(i).mutationMatrix];
    mutationList = [mutationList, mutationData(i).mutationList];
end


% consolidate co-occuring mutations
clusters = clusterdata(mutationMatrix', 'cutoff', eps, 'distance', 'hamming');
[~, iMutation] = unique(clusters);
mutationMatrix = mutationMatrix(:, iMutation);
mutationListClustered = {};
for i = 1:length(iMutation)
    mutationListClustered{end+1} =...
        sprintf('%s;', mutationList{clusters == i});
end
mutationList = mutationListClustered;



function [mutationMatrix, mutationList] = getMutationMatrixAux(seqs, gene, referenceStrain, useMuscle, genome_id)

if useMuscle
    fastawrite('~/tmp/test.fasta', seqs);
    [status,cmdout] = system('muscle -in ~/tmp/test.fasta -msf -out ~/tmp/test.msf -maxiters 1 -diags1 -sv -distance1 kbit20_3');
    %[status,cmdout] = system('muscle -in ~/tmp/test.fasta -msf -out ~/tmp/test.msf -maxiters 1 -diags1 -sv -distance1 kbit20_3');
    SeqsMultiAligned = multialignread('~/tmp/test.msf');
    
    [~, iSort] = sort(cellfun(@str2num,{SeqsMultiAligned.Header}));
    SeqsMultiAligned = SeqsMultiAligned(iSort);
    system('rm ~/tmp/test.*');
else
    SeqsMultiAligned = multialign(seqs);
end
% find the mutations


if strcmp(referenceStrain, 'consensus')
    % replaced this code by my own calculation of consensus
    % because Matlab' seqconsensus was giving inconsistent results for wspC
    %consensus = seqconsensus(SeqsMultiAligned);
    
    % get the consensus sequence
    consensus = zeros(1, length(SeqsMultiAligned(1).Sequence));
    for i = 1:length(consensus)
        seqs = char(1, length(SeqsMultiAligned));
        for j = 1:length(SeqsMultiAligned)
            seqs(j) = SeqsMultiAligned(j).Sequence(i);
        end
        consensus(i) = mode(double(seqs));
    end
    consensus = char(consensus);
    reference = consensus;
else
    % use the reference strain provided as input
    reference = SeqsMultiAligned(genome_id == referenceStrain).Sequence;
end

% make list of mutations
mutationList = {};
for i = 1:length(SeqsMultiAligned)
    mutationsInStrain = find(SeqsMultiAligned(i).Sequence ~= reference);
    for j = 1:length(mutationsInStrain)
        r = mutationsInStrain(j);
        m =...
            sprintf('%s %s%d%s', gene, reference(r), r, SeqsMultiAligned(i).Sequence(r));
        if ~ismember(m, mutationList)
            mutationList{end+1} = m;
        end
    end
end
% make the mutation matrix
mutationMatrix = zeros(length(SeqsMultiAligned), length(mutationList));
for i = 1:length(SeqsMultiAligned)
    mutationsInStrain = find(SeqsMultiAligned(i).Sequence ~= reference);
    for j = 1:length(mutationsInStrain)
        r = mutationsInStrain(j);
        m =...
            sprintf('%s %s%d%s', gene, reference(r), r, SeqsMultiAligned(i).Sequence(r));
        if ismember(m, mutationList)
            mutationMatrix(i, strcmp(m,mutationList)) = 1;
        end
    end
end





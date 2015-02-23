function [mutationMatrix, mutationList] = getMutationMatrix(db, geneList, referenceStrain, useMuscle)

% parse inputs
% if no genome_id is provided as reference strain then use the consensus
% seuqence
if nargin == 2
    referenceStrain = 'consensus';
    useMuscle = false;
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
for i = 1:length(geneList)
    [mm, ml] = getMutationMatrixAux(db, geneList{i}, referenceStrain, useMuscle);
    mutationMatrix = [mutationMatrix, mm];
    mutationList = [mutationList, ml];
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



function [mutationMatrix, mutationList] = getMutationMatrixAux(db, gene, referenceStrain, useMuscle)

seqsQuery = cell(db.getSequencesOfGene(gene));
genome_id = db.getPhenotypeColumn('genome_id');

for i = 1:length(seqsQuery)
    aaseq = nt2aa(seqsQuery{i});
    seqs(i).Sequence = aaseq(1:end-1); % must translate
    seqs(i).Header = num2str(genome_id(i));
end


if useMuscle
    fastawrite('~/tmp/test.fasta', seqs);
    [status,cmdout] = system('muscle -in ~/tmp/test.fasta -msf -out ~/tmp/test.msf');
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
    g = db.getPhenotypeColumn('genome_id');
    reference = SeqsMultiAligned(g == referenceStrain).Sequence;
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





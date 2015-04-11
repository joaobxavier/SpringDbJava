function [SeqsMultiAligned] = makeAlignment(db, gene)

genome_id = db.getPhenotypeColumn('genome_id');


if isnumeric(gene)
    seqsQuery = cell(db.getSequencesOfGenesFromId(gene));
else
    seqsQuery = cell(db.getSequencesOfGene(gene));
end

% get all sequences in a single query (great speed improvemnt)
for i = 1:size(seqsQuery, 1) % iterate over genome list
    aaseq = nt2aa(seqsQuery{i, 1});
    seqs(i).Sequence = aaseq(1:end-1); % must translate
    seqs(i).Header = num2str(genome_id(i));
end

fastawrite('~/tmp/test.fasta', seqs);
[status,cmdout] = system('muscle -in ~/tmp/test.fasta -msf -out ~/tmp/test.msf -maxiters 1 -diags1 -sv -distance1 kbit20_3');
SeqsMultiAligned = multialignread('~/tmp/test.msf');
[~, iSort] = sort(cellfun(@str2num,{SeqsMultiAligned.Header}));
SeqsMultiAligned = SeqsMultiAligned(iSort);
system('rm ~/tmp/test.*');

% replace the character used for gaps
for i = 1:length(SeqsMultiAligned)
    SeqsMultiAligned(i).Sequence =...
        strrep(SeqsMultiAligned(i).Sequence, '.', '-');
end










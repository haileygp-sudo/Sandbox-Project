function codons = sequenceCodons(dnaSequence)

    dnaSequence = upper(dnaSequence); % makes all letters uppercase
    seqLength = length(dnaSequence);

    numberCodons = floor(seqLength / 3); % number of complete codons
    codons = cell(1, numberCodons);

    for i = 1:numberCodons

        % start and end of each codon
        startIdx = (i-1)*3 + 1;
        endIdx = startIdx + 2;
        codons{i} = dnaSequence(startIdx:endIdx);
    end
end

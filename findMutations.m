% find mutations
function mutations = findMutations(seq1, seq2)

    lend = min(length(seq1), length(seq2));
    mutationIndex = [];
    healthyBase = [];
    patientBase = [];

    for i = 1:len
        if seq1(i) ~= seq2(i)
            mutationIndex(end + 1) = i;
            healthyBase(end + 1) = seq1(i);
            patientBase(end + 1) = seq2(i);
        end
    end
 
    mutations.index = mutationIndex;
    mutations.healthy = healthyBase;
    mutations.patient = patientBase;
end 